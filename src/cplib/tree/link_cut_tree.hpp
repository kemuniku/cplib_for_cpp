#pragma once
#include <cplib/tree/private/link_cut_tree_base.hpp>

namespace cplib::detail {
template <class S> struct LinkCutTreeState : LinkCutTreeOperations<LinkCutTreeState<S>, S> {
    struct Node {
        Int left = 0, right = 0, parent = 0;
        bool rev = false;
        S value, prod, rprod, virtualValue, all;
    };

    std::vector<Node> nodes;
    std::vector<Int> stack;
    std::function<S(S, S)> merge;
    S defaultValue;
    std::function<S(S)> inverse;

    LinkCutTreeState(std::span<const S> values, std::function<S(S, S)> op, S e,
                     std::function<S(S)> inv)
        : merge(std::move(op)), defaultValue(e), inverse(std::move(inv)) {
        nodes.reserve(values.size() + 1);
        nodes.push_back({0, 0, 0, false, e, e, e, e, e});
        for (const auto &x : values)
            nodes.push_back({0, 0, 0, false, x, x, x, e, x});
    }

    void push(Int v) {
        if (!v || !nodes[v].rev)
            return;
        for (Int child : {nodes[v].left, nodes[v].right})
            this->toggle(child);
        nodes[v].rev = false;
    }

    void setParent(Int v, Int parent) {
        if (v)
            nodes[v].parent = parent;
    }

    void pull(Int v) {
        Int l = nodes[v].left, r = nodes[v].right;
        nodes[v].prod = merge(merge(nodes[l].prod, nodes[v].value), nodes[r].prod);
        nodes[v].rprod = merge(merge(nodes[r].rprod, nodes[v].value), nodes[l].rprod);
        if (inverse)
            nodes[v].all = merge(merge(nodes[l].all, nodes[r].all),
                                 merge(nodes[v].value, nodes[v].virtualValue));
    }

    void addVirtual(Int v, Int child) {
        if (inverse && child)
            nodes[v].virtualValue = merge(nodes[v].virtualValue, nodes[child].all);
    }

    void removeVirtual(Int v, Int child) {
        if (inverse && child)
            nodes[v].virtualValue = merge(nodes[v].virtualValue, inverse(nodes[child].all));
    }
};
}

namespace cplib {
// パスの非可換モノイド、または部分木・成分の可換群を管理するLink-Cut Tree。
template <class S>
class LinkCutTree : public detail::LinkCutTreeHandle<detail::LinkCutTreeState<S>, S> {
    using State = detail::LinkCutTreeState<S>;
    using Base = detail::LinkCutTreeHandle<State, S>;

public:
    // 頂点の値valuesを持つ、辺のない森を作る。時間・空間O(N)。
    // mergeとdefaultValueはモノイドをなすこと。パスは非可換でもよい。
    // 部分木・成分を集約する場合はinverseも渡し、可換群であることを利用者が保証する。
    // inverseを省略すると部分木情報を更新しない。各演算の計算量はmerge等がO(1)の場合。
    LinkCutTree(std::span<const S> values, std::function<S(S, S)> merge, S defaultValue,
                std::function<S(S)> inverse = {})
        : Base(std::make_shared<State>(values, std::move(merge), std::move(defaultValue),
                                       std::move(inverse))) {
    }

    // 全頂点の値がdefaultValueの、辺のない森を作る。時間・空間O(N)。
    LinkCutTree(Int n, std::function<S(S, S)> merge, S defaultValue,
                std::function<S(S)> inverse = {})
        : LinkCutTree(std::vector<S>(checked(n), defaultValue), std::move(merge), defaultValue,
                      std::move(inverse)) {
    }

private:
    static Int checked(Int n) {
        assert(n >= 0);
        return n;
    }
};

template <class V, class Op, class S, class Inv = std::function<S(S)>>
auto initLinkCutTree(V &&values, Op merge, S defaultValue, Inv inverse = {}) {
    return LinkCutTree<S>(std::forward<V>(values), std::move(merge), std::move(defaultValue),
                          std::move(inverse));
}

template <class... A> auto newLinkCutTreeWith(A &&...args) {
    return initLinkCutTree(std::forward<A>(args)...);
}
}
