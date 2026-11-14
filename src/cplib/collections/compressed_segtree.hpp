#pragma once
#include <cplib/collections/segtree.hpp>
#include <cplib/collections/compressed_coordinates_internal.hpp>

namespace cplib {
template <class K, class T> class CompressedSegmentTree {
    std::vector<K> coords;
    SegmentTree<T> tree;
    T identity;
    std::function<T(T, T)> merge;
    std::vector<Int> indexSlots;
    using Self = CompressedSegmentTree;
    std::function<void(Self &, Int, T)> updateImpl;
    std::function<T(const Self &, K, K, bool)> rangeImpl;

    Int coordinateIndex(const K &x) const {
        return detail::findCompressedCoordinate<K>(coords, indexSlots, x);
    }

public:
    // 入力座標数Nに対しO(N log(N+1))。演算・比較・initialはO(1)を想定する。
    // 座標をソート・重複除去し、登録点をinitial(x)（省略時はe）で初期化します。
    // 構築後の座標追加はできません。Kには一貫した<と==が必要です。
    template <class Op, class Initial = std::nullptr_t>
    CompressedSegmentTree(std::span<const K> input, Op op, T e, Initial initial = nullptr)
        : coords(detail::uniqueCompressedCoordinates(input)), tree(Int(coords.size()), op, e),
          identity(e), merge(op), indexSlots(detail::initCompressedCoordinateIndex<K>(coords)) {
        if constexpr (!std::is_same_v<Initial, std::nullptr_t>) {
            std::vector<T> values;
            values.reserve(coords.size());
            for (const K &x : coords)
                values.push_back(initial(x));
            tree = initSegmentTree(values, op, e);
        }
    }

    // 座標の二分探索と部分木積の取得を同じ走査で行う。O(log N)。
    template <bool Inclusive, class Op> T rangeProductBody(const K &l, const K &r, Op op) const {
        auto beforeEnd = [&](const K &x) {
            if constexpr (Inclusive)
                return !(r < x);
            else
                return x < r;
        };
        Int n = coords.size();
        if (n == 0)
            return identity;
        if (coords[n - 1] < l || !beforeEnd(coords[0]))
            return identity;
        if (!(coords[0] < l) && beforeEnd(coords[n - 1]))
            return tree.arr[1];
        Int node = 1, a = 0, b = tree.arr.size() / 2;
        while (b - a > 1) {
            Int m = (a + b) >> 1;
            if (m >= n || !beforeEnd(coords[m])) {
                node <<= 1;
                b = m;
            } else if (coords[m] < l) {
                node = (node << 1) | 1;
                a = m;
            } else {
                T left = identity, right = identity;
                Int ln = node << 1, la = a, lb = m;
                while (lb - la > 1) {
                    Int mid = (la + lb) >> 1;
                    if (coords[mid] < l) {
                        ln = (ln << 1) | 1;
                        la = mid;
                    } else {
                        left = op(tree.arr[(ln << 1) | 1], left);
                        ln <<= 1;
                        lb = mid;
                    }
                }
                if (!(coords[la] < l))
                    left = op(tree.arr[ln], left);
                Int rn = (node << 1) | 1, ra = m, rb = b;
                while (rb - ra > 1) {
                    Int mid = (ra + rb) >> 1;
                    if (mid >= n || !beforeEnd(coords[mid])) {
                        rn <<= 1;
                        rb = mid;
                    } else {
                        right = op(right, tree.arr[rn << 1]);
                        rn = (rn << 1) | 1;
                        ra = mid;
                    }
                }
                if (beforeEnd(coords[ra]))
                    right = op(right, tree.arr[rn]);
                return op(left, right);
            }
        }
        return !(coords[a] < l) && beforeEnd(coords[a]) ? tree.arr[node] : identity;
    }

    // 登録済みの座標xの値をO(log N)で上書きします。
    void update(const K &x, const T &value) {
        Int i = coordinateIndex(x);
        assert(i >= 0);
        if (updateImpl)
            updateImpl(*this, i, value);
        else
            tree.update(i, value);
    }

    // 座標xの値をO(log N)で返します。未登録の座標では単位元を返します。
    T point(const K &x) const {
        Int i = coordinateIndex(x);
        return i >= 0 ? tree[i] : identity;
    }

    // 座標xの値を返す。未登録なら単位元。O(log(N+1))。
    T operator[](const K &x) const {
        return point(x);
    }

    struct Reference {
        Self *tree;
        K key;

        // 現在値を読み取る。O(log(N+1))。
        operator T() const {
            return tree->point(key);
        }

        // 値を代入する。O(log(N+1))。
        Reference &operator=(const T &x) {
            tree->update(key, x);
            return *this;
        }

        // 別の参照の値を読み取り、この位置へ代入する。O(log(N+1))。
        Reference &operator=(const Reference &r) {
            return *this = T(r);
        }
    };

    // 座標xへの読み書き用参照を返す。取得O(1)、読み書きO(log(N+1))。
    Reference operator[](const K &x) {
        return {this, x};
    }

    // 半開区間[l,r)内の登録点の積を座標順にO(log N)で返します。両端は未登録でも構いません。
    T get(const K &l, const K &r) const {
        assert(!(r < l));
        return rangeImpl ? rangeImpl(*this, l, r, false) : rangeProductBody<false>(l, r, merge);
    }

    // 閉区間内の登録点の積をO(log N)で返します。逆順の区間では単位元を返します。
    T get(ClosedSlice<K, K> s) const {
        if (s.b < s.a)
            return identity;
        return rangeImpl ? rangeImpl(*this, s.a, s.b, true)
                         : rangeProductBody<true>(s.a, s.b, merge);
    }

    // 閉区間内の登録点の積を返す。O(log(N+1))。
    T operator[](ClosedSlice<K, K> s) const {
        return get(s);
    }

    // 登録点の全体の積を返す。O(1)。
    T get_all() const {
        return tree.get_all();
    }

    // 重複除去後の座標数を返す。O(1)。
    Int size() const {
        return coords.size();
    }

    // 圧縮順の全値を文字列で返す。出力長LとしてO(N+L)。
    std::string str() const {
        return tree.str();
    }

    // 内部のセグメント木配列への参照を返す。O(1)。
    std::vector<T> &compressedData() {
        return tree.arr;
    }

    // 内部のセグメント木配列への読み取り専用参照を返す。O(1)。
    const std::vector<T> &compressedData() const {
        return tree.arr;
    }

    // 重複除去済み座標列への参照を返す。O(1)。
    const std::vector<K> &compressedCoordinates() const {
        return coords;
    }

    // 単位元を返す。値のコピーがO(1)ならO(1)。
    T compressedIdentity() const {
        return identity;
    }

    // 更新と区間取得の実装を差し替える。関数オブジェクトのコピー費用を除きO(1)。
    template <class Update, class Range>
    void setCompressedOperations(Update updateOp, Range rangeOp) {
        updateImpl = updateOp;
        rangeImpl = rangeOp;
    }
};

// 座標数Nから圧縮セグメント木を構築する。比較・演算・初期化がO(1)ならO(N log(N+1))。
template <class K, class T, class Op, class Initial = std::nullptr_t>
auto initCompressedSegmentTree(std::span<const K> coords, Op op, T e, Initial initial = nullptr) {
    return CompressedSegmentTree<K, T>(coords, op, e, initial);
}

// 座標数Nから圧縮セグメント木を構築する。比較・演算・初期化がO(1)ならO(N log(N+1))。
template <class K, class T, class Op, class Initial = std::nullptr_t>
auto initCompressedSegmentTree(const std::vector<K> &coords, Op op, T e,
                               Initial initial = nullptr) {
    return CompressedSegmentTree<K, T>(coords, op, e, initial);
}

// 式で指定する版は、区間内の演算を直接呼ぶ専用処理を登録する。
// 座標を圧縮し、直接演算するセグメント木を構築する。O(N log(N+1))。
template <class Coords, class Op, class T, class Initial = std::nullptr_t>
auto newCompressedSegWith(const Coords &coords, Op op, T e, Initial initial = nullptr) {
    auto result = initCompressedSegmentTree(coords, op, e, initial);
    result.setCompressedOperations(
        [op](auto &self, Int i, T value) {
            auto &arr = self.compressedData();
            Int node = i + arr.size() / 2;
            arr[node] = value;
            while (node > 1) {
                node >>= 1;
                arr[node] = op(arr[2 * node], arr[2 * node + 1]);
            }
        },
        [op](const auto &self, auto l, auto r, bool inclusive) {
            return inclusive ? self.template rangeProductBody<true>(l, r, op)
                             : self.template rangeProductBody<false>(l, r, op);
        });
    return result;
}

// 指定区間の積を返す。O(log(N+1))。
template <class K, class T> auto get(const CompressedSegmentTree<K, T> &s, const K &l, const K &r) {
    return s.get(l, r);
}

// 指定区間の積を返す。O(log(N+1))。
template <class K, class T> auto get(const CompressedSegmentTree<K, T> &s, ClosedSlice<K, K> r) {
    return s.get(r);
}

// 全体の積を返す。O(1)。
template <class K, class T> auto get_all(const CompressedSegmentTree<K, T> &s) {
    return s.get_all();
}

// 登録座標の値を上書きする。O(log(N+1))。
template <class K, class T> void update(CompressedSegmentTree<K, T> &s, const K &p, const T &x) {
    s.update(p, x);
}

// 圧縮順の全値を文字列で返す。出力長LとしてO(N+L)。
template <class K, class T> auto to_string(const CompressedSegmentTree<K, T> &s) {
    return s.str();
}
}
