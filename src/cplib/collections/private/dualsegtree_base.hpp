#pragma once
#include <cplib/common.hpp>
#include <cplib/utils/backwards_index.hpp>
#include <sstream>

namespace cplib::detail::dualsegtree_base {
// 区間作用・一点取得O(log N)、構築O(N)。通常版とhasLazy版は別の特殊化。
template <class S, class F, class Mapping = std::function<S(F, S)>,
          class Composition = std::function<F(F, F)>, bool HasFlag = false>
class DualSegmentTree {
    std::vector<S> data;
    std::vector<F> lazy;
    std::vector<bool> hasLazy;
    [[no_unique_address]] Mapping mapping;
    [[no_unique_address]] Composition composition;
    F id;
    Int lastnode = 1, log = 0, length = 0;

    void allApply(Int p, const F &f) {
        if (p < lastnode) {
            if constexpr (HasFlag) {
                if (hasLazy[p])
                    lazy[p] = composition(f, lazy[p]);
                else {
                    lazy[p] = f;
                    hasLazy[p] = true;
                }
            } else
                lazy[p] = composition(f, lazy[p]);
        } else {
            Int i = p - lastnode;
            if (i < length)
                data[i] = mapping(f, data[i]);
        }
    }

    [[gnu::noinline]] void pushFlag(Int p) {
        if (hasLazy[p]) {
            F f = lazy[p];
            allApply(2 * p, f);
            allApply(2 * p + 1, f);
            hasLazy[p] = false;
        }
    }

    void push(Int p) {
        if constexpr (HasFlag)
            pushFlag(p);
        else {
            allApply(2 * p, lazy[p]);
            allApply(2 * p + 1, lazy[p]);
            lazy[p] = id;
        }
    }

    void allPush(Int p) {
        for (Int i = log; i >= 1; --i)
            push(p >> i);
    }

public:
    // 初期配列と作用・合成を指定して構築する。要素数をnとしてO(n+1)。
    DualSegmentTree(std::span<const S> v, Mapping map, Composition comp, F id_)
        : data(v.begin(), v.end()), mapping(map), composition(comp), id(id_), length(v.size()) {
        while (lastnode < length) {
            lastnode *= 2;
            ++log;
        }
        if constexpr (HasFlag) {
            lazy.resize(lastnode);
            hasLazy.resize(lastnode);
        } else
            lazy.assign(lastnode, id);
    }

    // 同じ初期値を持つ長さnの列から構築する。O(n+1)。
    DualSegmentTree(Int n, S initial, Mapping map, Composition comp, F id_)
        : DualSegmentTree(std::vector<S>(n, initial), map, comp, id_) {
        assert(n >= 0);
    }

    // O(log(n+1))回の作用・合成。
    // 半開区間[l,r)の各要素にfを作用させます。
    void apply(Int l, Int r, const F &f) {
        assert(0 <= l && l <= r && r <= length);
        if (l == r)
            return;
        l += lastnode;
        r += lastnode;
        for (Int i = log; i >= 1; --i) {
            if (((l >> i) << i) != l)
                push(l >> i);
            if (((r >> i) << i) != r)
                push((r - 1) >> i);
        }
        while (l < r) {
            if (l & 1)
                allApply(l++, f);
            if (r & 1)
                allApply(--r, f);
            l >>= 1;
            r >>= 1;
        }
    }

    // O(log(n+1))回の作用・合成。
    // 閉区間segmentの各要素にfを作用させます。
    template <class L, class R> void apply(ClosedSlice<L, R> s, const F &f) {
        apply(resolve_index(length, s.a), resolve_index(length, s.b) + 1, f);
    }

    // i番目の値を取得する。O(log(n+1))回の作用・合成。
    S get(Int i) {
        assert(0 <= i && i < length);
        allPush(i + lastnode);
        return data[i];
    }

    // i番目の値を設定する。O(log(n+1))回の作用・合成。
    void update(Int i, const S &value) {
        assert(0 <= i && i < length);
        allPush(i + lastnode);
        data[i] = value;
    }

    // 要素数を返す。O(1)。
    Int size() const {
        return length;
    }

    struct Reference {
        DualSegmentTree *tree;
        Int index;

        // 対象の値を取得する。O(log(n+1))回の作用・合成。
        operator S() const {
            return tree->get(index);
        }

        // 対象の値を設定する。O(log(n+1))回の作用・合成。
        Reference &operator=(const S &v) {
            tree->update(index, v);
            return *this;
        }

        // 他の位置の値を読み取り代入する。O(log(n+1))回の作用・合成。
        Reference &operator=(const Reference &v) {
            return *this = S(v);
        }
    };

    // i番目を読み書きする参照を返す。O(1)。
    Reference operator[](Int i) {
        assert(0 <= i && i < length);
        return {this, i};
    }

    // 末尾からi番目を読み書きする参照を返す。O(1)。
    Reference operator[](BackwardsIndex i) {
        return (*this)[length - i.value];
    }

    struct Iterator {
        DualSegmentTree *tree;
        Int index;

        // 指す位置の値を取得する。O(log(n+1))回の作用・合成。
        S operator*() const {
            return tree->get(index);
        }

        // 次の位置へ進む。O(1)。
        Iterator &operator++() {
            ++index;
            return *this;
        }

        // 位置が異なるか返す。O(1)。
        bool operator!=(const Iterator &b) const {
            return index != b.index;
        }
    };

    // 先頭のイテレータを返す。O(1)。
    Iterator begin() {
        return {this, 0};
    }

    // 末尾の次のイテレータを返す。O(1)。
    Iterator end() {
        return {this, length};
    }

    // 列挙用に自身への参照を返す。O(1)。
    DualSegmentTree &items() {
        return *this;
    }

    // 全要素を配列にして返す。O(n log(n+1))回の作用・合成。
    std::vector<S> toSeq() {
        std::vector<S> result;
        result.reserve(length);
        for (Int i = 0; i < length; ++i)
            result.push_back(get(i));
        return result;
    }

    // 全要素を空白区切りで文字列化する。O(n log(n+1))回の作用・合成と出力処理。
    std::string str() {
        std::ostringstream s;
        for (Int i = 0; i < length; ++i) {
            if (i)
                s << ' ';
            s << get(i);
        }
        return s.str();
    }
};

// 配列と作用・合成を指定して構築する。要素数をnとしてO(n+1)。
template <class S, class F, class Map, class Comp>
auto initDualSegmentTree(std::span<const S> v, Map map, Comp comp, F id) {
    return DualSegmentTree<S, F>(v, map, comp, id);
}

// 配列と作用・合成を指定して構築する。要素数をnとしてO(n+1)。
template <class S, class F, class Map, class Comp>
auto initDualSegmentTree(const std::vector<S> &v, Map map, Comp comp, F id) {
    return DualSegmentTree<S, F>(v, map, comp, id);
}

// 長さnの同じ初期値を持つ列から構築する。O(n+1)。
template <class S, class F, class Map, class Comp>
auto initDualSegmentTree(Int n, S initial, Map map, Comp comp, F id) {
    return DualSegmentTree<S, F>(n, initial, map, comp, id);
}

// 初期列と作用・合成を指定して構築する。要素数をnとしてO(n+1)。
template <class... Args> auto newDualSegWith(Args &&...args) {
    return initDualSegmentTree(std::forward<Args>(args)...);
}

#define CPLIB_DUAL_TEMPLATE template <class S, class F, class M, class C, bool B>
#define CPLIB_DUAL_TYPE DualSegmentTree<S, F, M, C, B>

// i番目の値を取得する。O(log(n+1))回の作用・合成。
CPLIB_DUAL_TEMPLATE auto get(CPLIB_DUAL_TYPE &s, Int i) {
    return s.get(i);
}

// i番目の値を設定する。O(log(n+1))回の作用・合成。
CPLIB_DUAL_TEMPLATE void update(CPLIB_DUAL_TYPE &s, Int i, const S &value) {
    s.update(i, value);
}

// 半開区間[l,r)にfを作用させる。O(log(n+1))回の作用・合成。
CPLIB_DUAL_TEMPLATE void apply(CPLIB_DUAL_TYPE &s, Int l, Int r, const F &f) {
    s.apply(l, r, f);
}

// 全要素を配列にして返す。O(n log(n+1))回の作用・合成。
CPLIB_DUAL_TEMPLATE auto toSeq(CPLIB_DUAL_TYPE &s) {
    return s.toSeq();
}

// 列挙用にsへの参照を返す。O(1)。
CPLIB_DUAL_TEMPLATE auto &items(CPLIB_DUAL_TYPE &s) {
    return s;
}

// 全要素を空白区切りで文字列化する。O(n log(n+1))回の作用・合成と出力処理。
CPLIB_DUAL_TEMPLATE auto to_string(CPLIB_DUAL_TYPE &s) {
    return s.str();
}

// 両端を含む区間にfを作用させる。O(log(n+1))回の作用・合成。
template <class S, class F, class M, class C, bool B, class L, class R>
void apply(CPLIB_DUAL_TYPE &s, ClosedSlice<L, R> r, const F &f) {
    s.apply(r, f);
}

#undef CPLIB_DUAL_TEMPLATE
#undef CPLIB_DUAL_TYPE
}
