#pragma once
#include <cplib/common.hpp>

namespace cplib {
// 登録点を圧縮した2次元Fenwick。構築O(N log N)、更新・矩形和O(log² N)。
template <class K, class T = Int> class CompressedFenwick2D {
    std::vector<K> xs, ys, pointYs;
    std::vector<Int> offsets, pointOffsets;
    std::vector<T> data, pointValues;

    static Int lowerIndex(const std::vector<K> &c, Int start, Int end, const K &y) {
        return std::lower_bound(c.begin() + start, c.begin() + end, y) - (c.begin() + start);
    }

    Int yIndex(Int node, const K &y) const {
        return lowerIndex(ys, offsets[node], offsets[node + 1], y);
    }

    Int pointIndex(Int xi, const K &y) const {
        Int start = pointOffsets[xi], end = pointOffsets[xi + 1],
            i = start + lowerIndex(pointYs, start, end, y);
        return i < end && pointYs[i] == y ? i : -1;
    }

    void addImpl(Int xi, const K &y, const T &delta) {
        for (Int node = xi + 1; node <= Int(xs.size()); node += node & -node) {
            Int start = offsets[node], size = offsets[node + 1] - start;
            for (Int i = yIndex(node, y); i < size; i |= i + 1)
                data[start + i] += delta;
        }
    }

    T innerSum(Int node, Int l, Int r) const {
        Int start = offsets[node];
        T result{}, left{};
        while (r > l) {
            result += data[start + r - 1];
            r &= r - 1;
        }
        while (l > r) {
            left += data[start + l - 1];
            l &= l - 1;
        }
        return result - left;
    }

    template <class P> void build(std::span<const P> points) {
        constexpr bool weighted = std::tuple_size_v<P> == 3;
        std::vector<P> ps(points.begin(), points.end());
        auto cmp = [](const P &a, const P &b) {
            if (std::get<0>(a) < std::get<0>(b))
                return true;
            if (std::get<0>(b) < std::get<0>(a))
                return false;
            return std::get<1>(a) < std::get<1>(b);
        };
        std::sort(ps.begin(), ps.end(), cmp);
        Int count = 0;
        for (Int i = 0; i < Int(ps.size()); ++i) {
            if (count == 0 || std::get<0>(ps[count - 1]) != std::get<0>(ps[i]) ||
                std::get<1>(ps[count - 1]) != std::get<1>(ps[i]))
                ps[count++] = ps[i];
            else if constexpr (weighted)
                std::get<2>(ps[count - 1]) += std::get<2>(ps[i]);
        }
        ps.resize(count);
        pointYs.resize(count);
        pointValues.resize(count);
        std::vector<std::pair<K, Int>> ranked(count);
        for (Int i = 0; i < count; ++i) {
            const auto &p = ps[i];
            if (xs.empty() || xs.back() != std::get<0>(p)) {
                xs.push_back(std::get<0>(p));
                pointOffsets.push_back(i);
            }
            pointYs[i] = std::get<1>(p);
            if constexpr (weighted)
                pointValues[i] = std::get<2>(p);
            ranked[i] = {std::get<1>(p), Int(xs.size())};
        }
        pointOffsets.push_back(count);
        Int n = xs.size();
        std::sort(ranked.begin(), ranked.end(),
                  [](const auto &a, const auto &b) { return a.first < b.first; });
        offsets.assign(n + 2, 0);
        std::vector<Int> cursor(n + 1, -1);
        // 内側の座標の総数を数え、平坦な配列を一度だけ確保する。
        for (Int i = 0; i < count; ++i) {
            auto [y, x] = ranked[i];
            for (Int node = x; node <= n; node += node & -node) {
                if (cursor[node] != -1 && ranked[cursor[node]].first == y)
                    break;
                ++offsets[node + 1];
                cursor[node] = i;
            }
        }
        for (Int i = 1; i <= n + 1; ++i)
            offsets[i] += offsets[i - 1];
        ys.resize(offsets[n + 1]);
        data.resize(ys.size());
        std::fill(cursor.begin(), cursor.end(), 0);
        std::vector<Int> pointCursor;
        if constexpr (weighted)
            pointCursor.resize(n);
        for (auto [y, x] : ranked) {
            T value{};
            if constexpr (weighted) {
                Int xi = x - 1;
                value = pointValues[pointOffsets[xi] + pointCursor[xi]++];
            }
            for (Int node = x; node <= n; node += node & -node) {
                Int start = offsets[node];
                if (cursor[node] > 0 && ys[start + cursor[node] - 1] == y) {
                    if constexpr (weighted)
                        data[start + cursor[node] - 1] += value;
                    else
                        break;
                } else {
                    Int index = start + cursor[node];
                    ys[index] = y;
                    if constexpr (weighted)
                        data[index] += value;
                    ++cursor[node];
                }
            }
        }
        if constexpr (weighted)
            for (Int node = 1; node <= n; ++node) {
                Int start = offsets[node], size = offsets[node + 1] - start;
                for (Int i = 0; i < size; ++i) {
                    Int parent = i | (i + 1);
                    if (parent < size)
                        data[start + parent] += data[start + i];
                }
            }
    }

public:
    // 更新点を事前登録し、全点を零で初期化します。重複点は一つにまとめます。
    // Kには一貫した<と==、TにはT{}を零とする可換加法群の演算が必要です。
    explicit CompressedFenwick2D(std::span<const std::pair<K, K>> p) {
        build(p);
    }

    // (x,y,w)の列から一括構築します。同じ座標の重みは加算し、重みが零の点も登録します。
    explicit CompressedFenwick2D(std::span<const std::tuple<K, K, T>> p) {
        build(p);
    }

    Int len() const {
        return pointYs.size();
    }

    Int size() const {
        return len();
    }

    // 登録点(x,y)にdeltaをO(log² N)で加算します。未登録点は更新できません。
    void add(const K &x, const K &y, const T &delta) {
        Int xi = std::lower_bound(xs.begin(), xs.end(), x) - xs.begin();
        assert(xi < Int(xs.size()) && xs[xi] == x);
        Int pi = pointIndex(xi, y);
        assert(pi >= 0);
        pointValues[pi] += delta;
        addImpl(xi, y, delta);
    }

    // x<xUpperかつy<yUpperを満たす登録点の和をO(log² N)で返します。
    T prefix(const K &xUpper, const K &yUpper) const {
        T result{};
        Int node = std::lower_bound(xs.begin(), xs.end(), xUpper) - xs.begin();
        while (node > 0) {
            result += innerSum(node, 0, yIndex(node, yUpper));
            node &= node - 1;
        }
        return result;
    }

    // xl<=x<xrかつy<yUpperを満たす登録点の和をO(log² N)で返します。
    T getLess(const K &xl, const K &xr, const K &yUpper) const {
        assert(!(xr < xl));
        Int l = std::lower_bound(xs.begin(), xs.end(), xl) - xs.begin(),
            r = std::lower_bound(xs.begin(), xs.end(), xr) - xs.begin();
        T result{}, left{};
        while (r > l) {
            result += innerSum(r, 0, yIndex(r, yUpper));
            r &= r - 1;
        }
        while (l > r) {
            left += innerSum(l, 0, yIndex(l, yUpper));
            l &= l - 1;
        }
        return result - left;
    }

    // [xl,xr)×[yl,yr)の和をO(log² N)で返します。境界は未登録でも構いません。
    T get(const K &xl, const K &xr, const K &yl, const K &yr) const {
        assert(!(xr < xl) && !(yr < yl));
        if (!(yl < yr))
            return T{};
        Int l = std::lower_bound(xs.begin(), xs.end(), xl) - xs.begin(),
            r = std::lower_bound(xs.begin(), xs.end(), xr) - xs.begin();
        T result{}, left{};
        while (r > l) {
            result += innerSum(r, yIndex(r, yl), yIndex(r, yr));
            r &= r - 1;
        }
        while (l > r) {
            left += innerSum(l, yIndex(l, yl), yIndex(l, yr));
            l &= l - 1;
        }
        return result - left;
    }

    // 点取得のみO(log N)。未登録なら零。
    T get(const K &x, const K &y) const {
        Int xi = std::lower_bound(xs.begin(), xs.end(), x) - xs.begin();
        if (xi == Int(xs.size()) || xs[xi] != x)
            return T{};
        Int pi = pointIndex(xi, y);
        return pi >= 0 ? pointValues[pi] : T{};
    }

    void set(const K &x, const K &y, const T &value) {
        Int xi = std::lower_bound(xs.begin(), xs.end(), x) - xs.begin();
        assert(xi < Int(xs.size()) && xs[xi] == x);
        Int pi = pointIndex(xi, y);
        assert(pi >= 0);
        T delta = value - pointValues[pi];
        pointValues[pi] = value;
        addImpl(xi, y, delta);
    }

    T operator[](std::pair<K, K> p) const {
        return get(p.first, p.second);
    }

    struct Reference {
        CompressedFenwick2D *tree;
        std::pair<K, K> p;

        operator T() const {
            return tree->get(p.first, p.second);
        }

        Reference &operator=(const T &value) {
            tree->set(p.first, p.second, value);
            return *this;
        }

        Reference &operator=(const Reference &r) {
            return *this = T(r);
        }
    };

    Reference operator[](std::pair<K, K> p) {
        return {this, p};
    }

    T get_all() const {
        T result{};
        Int node = xs.size();
        while (node > 0) {
            result += innerSum(node, 0, offsets[node + 1] - offsets[node]);
            node &= node - 1;
        }
        return result;
    }
};

template <class K, class T = Int> auto initCompressedFenwick2D(std::span<const std::pair<K, K>> p) {
    return CompressedFenwick2D<K, T>(p);
}

template <class K, class T = Int>
auto initCompressedFenwick2D(const std::vector<std::pair<K, K>> &p) {
    return CompressedFenwick2D<K, T>(p);
}

template <class K, class T = Int, std::size_t N>
auto initCompressedFenwick2D(const std::array<std::pair<K, K>, N> &p) {
    return CompressedFenwick2D<K, T>(p);
}

template <class K, class T> auto initCompressedFenwick2D(std::span<const std::tuple<K, K, T>> p) {
    return CompressedFenwick2D<K, T>(p);
}

template <class K, class T>
auto initCompressedFenwick2D(const std::vector<std::tuple<K, K, T>> &p) {
    return CompressedFenwick2D<K, T>(p);
}

template <class K, class T> Int len(const CompressedFenwick2D<K, T> &t) {
    return t.len();
}

template <class K, class T>
void add(CompressedFenwick2D<K, T> &t, const K &x, const K &y, const T &d) {
    t.add(x, y, d);
}

template <class K, class T>
auto prefix(const CompressedFenwick2D<K, T> &t, const K &x, const K &y) {
    return t.prefix(x, y);
}

template <class K, class T>
auto getLess(const CompressedFenwick2D<K, T> &t, const K &l, const K &r, const K &y) {
    return t.getLess(l, r, y);
}

template <class K, class T>
auto get(const CompressedFenwick2D<K, T> &t, const K &xl, const K &xr, const K &yl, const K &yr) {
    return t.get(xl, xr, yl, yr);
}

template <class K, class T> auto get_all(const CompressedFenwick2D<K, T> &t) {
    return t.get_all();
}
}
