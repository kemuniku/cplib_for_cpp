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

    DualSegmentTree(Int n, S initial, Mapping map, Composition comp, F id_)
        : DualSegmentTree(std::vector<S>(n, initial), map, comp, id_) {
        assert(n >= 0);
    }

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

    // 閉区間segmentの各要素にfを作用させます。
    template <class L, class R> void apply(ClosedSlice<L, R> s, const F &f) {
        apply(resolve_index(length, s.a), resolve_index(length, s.b) + 1, f);
    }

    S get(Int i) {
        assert(0 <= i && i < length);
        allPush(i + lastnode);
        return data[i];
    }

    void update(Int i, const S &value) {
        assert(0 <= i && i < length);
        allPush(i + lastnode);
        data[i] = value;
    }

    Int len() const {
        return length;
    }

    Int size() const {
        return len();
    }

    struct Reference {
        DualSegmentTree *tree;
        Int index;

        operator S() const {
            return tree->get(index);
        }

        Reference &operator=(const S &v) {
            tree->update(index, v);
            return *this;
        }

        Reference &operator=(const Reference &v) {
            return *this = S(v);
        }
    };

    Reference operator[](Int i) {
        assert(0 <= i && i < length);
        return {this, i};
    }

    Reference operator[](BackwardsIndex i) {
        return (*this)[length - i.value];
    }

    struct Iterator {
        DualSegmentTree *tree;
        Int index;

        S operator*() const {
            return tree->get(index);
        }

        Iterator &operator++() {
            ++index;
            return *this;
        }

        bool operator!=(const Iterator &b) const {
            return index != b.index;
        }
    };

    Iterator begin() {
        return {this, 0};
    }

    Iterator end() {
        return {this, length};
    }

    DualSegmentTree &items() {
        return *this;
    }

    std::vector<S> toSeq() {
        std::vector<S> result;
        result.reserve(length);
        for (Int i = 0; i < length; ++i)
            result.push_back(get(i));
        return result;
    }

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

template <class S, class F, class Map, class Comp>
auto initDualSegmentTree(std::span<const S> v, Map map, Comp comp, F id) {
    return DualSegmentTree<S, F>(v, map, comp, id);
}

template <class S, class F, class Map, class Comp>
auto initDualSegmentTree(const std::vector<S> &v, Map map, Comp comp, F id) {
    return DualSegmentTree<S, F>(v, map, comp, id);
}

template <class S, class F, class Map, class Comp>
auto initDualSegmentTree(Int n, S initial, Map map, Comp comp, F id) {
    return DualSegmentTree<S, F>(n, initial, map, comp, id);
}

template <class... Args> auto newDualSegWith(Args &&...args) {
    return initDualSegmentTree(std::forward<Args>(args)...);
}

#define CPLIB_DUAL_TEMPLATE template <class S, class F, class M, class C, bool B>
#define CPLIB_DUAL_TYPE DualSegmentTree<S, F, M, C, B>

CPLIB_DUAL_TEMPLATE Int len(const CPLIB_DUAL_TYPE &s) {
    return s.len();
}

CPLIB_DUAL_TEMPLATE auto get(CPLIB_DUAL_TYPE &s, Int i) {
    return s.get(i);
}

CPLIB_DUAL_TEMPLATE void update(CPLIB_DUAL_TYPE &s, Int i, const S &value) {
    s.update(i, value);
}

CPLIB_DUAL_TEMPLATE void apply(CPLIB_DUAL_TYPE &s, Int l, Int r, const F &f) {
    s.apply(l, r, f);
}

CPLIB_DUAL_TEMPLATE auto toSeq(CPLIB_DUAL_TYPE &s) {
    return s.toSeq();
}

CPLIB_DUAL_TEMPLATE auto &items(CPLIB_DUAL_TYPE &s) {
    return s;
}

CPLIB_DUAL_TEMPLATE auto to_string(CPLIB_DUAL_TYPE &s) {
    return s.str();
}

template <class S, class F, class M, class C, bool B, class L, class R>
void apply(CPLIB_DUAL_TYPE &s, ClosedSlice<L, R> r, const F &f) {
    s.apply(r, f);
}

#undef CPLIB_DUAL_TEMPLATE
#undef CPLIB_DUAL_TYPE
}
