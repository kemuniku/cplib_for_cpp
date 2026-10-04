#pragma once
#include <cplib/common.hpp>

namespace cplib {
// 可換モノイド用。内側の葉数を2冪に丸めず、座標M個・値2M個を保持する。
template <class K, class T> class CompressedSegmentTree2D {
    std::vector<K> xs, ys;
    std::vector<Int> offsets;
    std::vector<T> data;
    Int base = 1, count = 0;
    T identity;
    std::function<T(T, T)> merge;
    using Self = CompressedSegmentTree2D;
    std::function<void(Self &, K, K, T)> updateImpl;
    std::function<T(const Self &, K, K, K, K)> rangeImpl;

    Int yIndex(Int node, const K &y) const {
        Int start = offsets[node];
        return std::lower_bound(ys.begin() + start, ys.begin() + offsets[node + 1], y) -
               (ys.begin() + start);
    }

    T valueAt(Int node, const K &y) const {
        Int start = offsets[node], size = offsets[node + 1] - start, i = yIndex(node, y);
        return i < size && ys[start + i] == y ? data[2 * start + size + i] : identity;
    }

public:
    // 登録点、または初期値付き登録点からO(N log N)で一括構築する。
    // 更新座標を事前登録します。構築後の座標追加はできません。Kには一貫した<と==が必要です。
    // 座標だけを与えた場合は各点をeで初期化し、初期値付きの場合は同じ座標の値をopでマージします。
    // opには可換モノイドの演算、eには単位元を指定してください。
    template <class P, class Op>
    CompressedSegmentTree2D(std::span<const P> points, Op op, T e) : identity(e), merge(op) {
        constexpr bool weighted = std::tuple_size_v<P> == 3;
        std::vector<P> ps(points.begin(), points.end());
        std::sort(ps.begin(), ps.end(), [](const P &a, const P &b) {
            if (std::get<0>(a) < std::get<0>(b))
                return true;
            if (std::get<0>(b) < std::get<0>(a))
                return false;
            return std::get<1>(a) < std::get<1>(b);
        });
        for (Int i = 0; i < Int(ps.size()); ++i) {
            if (count == 0 || std::get<0>(ps[count - 1]) != std::get<0>(ps[i]) ||
                std::get<1>(ps[count - 1]) != std::get<1>(ps[i]))
                ps[count++] = ps[i];
            else if constexpr (weighted)
                std::get<2>(ps[count - 1]) = merge(std::get<2>(ps[count - 1]), std::get<2>(ps[i]));
        }
        ps.resize(count);
        std::vector<std::pair<K, Int>> ranked(count);
        for (Int i = 0; i < count; ++i) {
            const auto &p = ps[i];
            if (xs.empty() || xs.back() != std::get<0>(p))
                xs.push_back(std::get<0>(p));
            ranked[i] = {std::get<1>(p), Int(xs.size()) - 1};
        }
        while (base < Int(xs.size()))
            base *= 2;
        std::sort(ranked.begin(), ranked.end(),
                  [](const auto &a, const auto &b) { return a.first < b.first; });
        Int nodes = base * 2;
        offsets.assign(nodes + 1, 0);
        std::vector<Int> cursor(nodes, -1);
        for (Int i = 0; i < count; ++i) {
            auto [y, x] = ranked[i];
            for (Int node = base + x; node > 0; node >>= 1) {
                if (cursor[node] != -1 && ranked[cursor[node]].first == y)
                    break;
                ++offsets[node + 1];
                cursor[node] = i;
            }
        }
        for (Int i = 1; i <= nodes; ++i)
            offsets[i] += offsets[i - 1];
        ys.resize(offsets[nodes]);
        data.assign(2 * ys.size(), e);
        std::fill(cursor.begin(), cursor.end(), 0);
        for (auto [y, x] : ranked)
            for (Int node = base + x; node > 0; node >>= 1) {
                Int start = offsets[node];
                if (cursor[node] > 0 && ys[start + cursor[node] - 1] == y)
                    break;
                ys[start + cursor[node]++] = y;
            }
        if constexpr (weighted) {
            Int xi = 0, yi = 0;
            for (const auto &p : ps) {
                while (xs[xi] < std::get<0>(p)) {
                    ++xi;
                    yi = 0;
                }
                Int node = base + xi, start = offsets[node], size = offsets[node + 1] - start;
                data[2 * start + size + yi++] = std::get<2>(p);
            }
            for (Int node = base * 2 - 1; node > 0; --node) {
                Int start = offsets[node], size = offsets[node + 1] - start;
                if (node < base) {
                    Int ls = offsets[node * 2], rs = offsets[node * 2 + 1], ln = rs - ls,
                        rn = offsets[node * 2 + 2] - rs, li = 0, ri = 0;
                    for (Int i = 0; i < size; ++i) {
                        while (li < ln && ys[ls + li] < ys[start + i])
                            ++li;
                        while (ri < rn && ys[rs + ri] < ys[start + i])
                            ++ri;
                        T lv = e, rv = e;
                        if (li < ln && ys[ls + li] == ys[start + i])
                            lv = data[2 * ls + ln + li];
                        if (ri < rn && ys[rs + ri] == ys[start + i])
                            rv = data[2 * rs + rn + ri];
                        data[2 * start + size + i] = merge(lv, rv);
                    }
                }
                for (Int i = size - 1; i > 0; --i)
                    data[2 * start + i] =
                        merge(data[2 * start + 2 * i], data[2 * start + 2 * i + 1]);
            }
        }
    }

    template <class Op> void updateBody(const K &x, const K &y, T value, Op op) {
        Int xi = std::lower_bound(xs.begin(), xs.end(), x) - xs.begin();
        assert(xi < Int(xs.size()) && xs[xi] == x);
        Int node = base + xi, yi = yIndex(node, y);
        assert(yi < offsets[node + 1] - offsets[node] && ys[offsets[node] + yi] == y);
        T current = value;
        while (true) {
            Int start = offsets[node], pos = offsets[node + 1] - start + yi;
            data[2 * start + pos] = current;
            while (pos > 1) {
                pos >>= 1;
                data[2 * start + pos] =
                    op(data[2 * start + 2 * pos], data[2 * start + 2 * pos + 1]);
            }
            if (node == 1)
                break;
            current = op(current, valueAt(node ^ 1, y));
            node >>= 1;
            yi = yIndex(node, y);
        }
    }

    template <class Op>
    T rangeBody(const K &xl, const K &xr, const K &yl, const K &yr, Op op) const {
        assert(!(xr < xl) && !(yr < yl));
        T acc = identity;
        if (xl < xr && yl < yr) {
            Int l = std::lower_bound(xs.begin(), xs.end(), xl) - xs.begin() + base,
                r = std::lower_bound(xs.begin(), xs.end(), xr) - xs.begin() + base;
            auto consume = [&](Int node) {
                Int start = offsets[node], size = offsets[node + 1] - start,
                    a = yIndex(node, yl) + size, b = yIndex(node, yr) + size;
                while (a < b) {
                    if (a & 1)
                        acc = op(acc, data[2 * start + a++]);
                    if (b & 1)
                        acc = op(acc, data[2 * start + --b]);
                    a >>= 1;
                    b >>= 1;
                }
            };
            while (l < r) {
                if (l & 1)
                    consume(l++);
                if (r & 1)
                    consume(--r);
                l >>= 1;
                r >>= 1;
            }
        }
        return acc;
    }

    // 点更新・矩形取得はO(log² N)、点取得はO(log N)。
    // 登録点(x,y)をO(log² N)で上書きします。
    void update(const K &x, const K &y, T value) {
        if (updateImpl)
            updateImpl(*this, x, y, value);
        else
            updateBody(x, y, value, merge);
    }

    // [xl,xr)×[yl,yr)の積をO(log² N)で返します。境界は未登録でも構いません。
    T get(const K &xl, const K &xr, const K &yl, const K &yr) const {
        return rangeImpl ? rangeImpl(*this, xl, xr, yl, yr) : rangeBody(xl, xr, yl, yr, merge);
    }

    // 点の値をO(log N)で返します。未登録なら単位元です。
    T get(const K &x, const K &y) const {
        Int xi = std::lower_bound(xs.begin(), xs.end(), x) - xs.begin();
        return xi < Int(xs.size()) && xs[xi] == x ? valueAt(base + xi, y) : identity;
    }

    T operator[](std::pair<K, K> p) const {
        return get(p.first, p.second);
    }

    struct Reference {
        Self *tree;
        std::pair<K, K> p;

        operator T() const {
            return tree->get(p.first, p.second);
        }

        Reference &operator=(T value) {
            tree->update(p.first, p.second, value);
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
        return count == 0 ? identity : data[2 * offsets[1] + 1];
    }

    Int len() const {
        return count;
    }

    Int size() const {
        return len();
    }

    template <class Op> void specialize(Op op) {
        updateImpl = [op](Self &t, K x, K y, T value) { t.updateBody(x, y, value, op); };
        rangeImpl = [op](const Self &t, K xl, K xr, K yl, K yr) {
            return t.rangeBody(xl, xr, yl, yr, op);
        };
    }
};

template <class P, class Op, class T>
auto initCompressedSegmentTree2D(std::span<const P> points, Op op, T e) {
    using K = std::tuple_element_t<0, P>;
    return CompressedSegmentTree2D<K, T>(points, op, e);
}

template <class P, class Op, class T>
auto initCompressedSegmentTree2D(const std::vector<P> &points, Op op, T e) {
    return initCompressedSegmentTree2D(std::span<const P>(points), op, e);
}

template <class Points, class Op, class T>
auto newCompressedSeg2DWith(const Points &points, Op op, T e) {
    auto result = initCompressedSegmentTree2D(points, op, e);
    result.specialize(op);
    return result;
}

template <class K, class T> Int len(const CompressedSegmentTree2D<K, T> &s) {
    return s.len();
}

template <class K, class T> T get_all(const CompressedSegmentTree2D<K, T> &s) {
    return s.get_all();
}

template <class K, class T> T get(const CompressedSegmentTree2D<K, T> &s, K xl, K xr, K yl, K yr) {
    return s.get(xl, xr, yl, yr);
}

template <class K, class T> void update(CompressedSegmentTree2D<K, T> &s, K x, K y, T v) {
    s.update(x, y, v);
}
}
