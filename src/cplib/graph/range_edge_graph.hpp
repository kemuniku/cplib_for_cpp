#pragma once
#include <cplib/graph/graph.hpp>

namespace cplib {
template <class T = Int> class Range_Edge_Graph {
    struct Node {
        Int l, r, left = -1, right = -1, inId, outId;
    };

    T zero;
    Int root = -1;
    std::vector<Node> nodes;

    void validateRange(Int l, Int r) const {
        assert(0 <= l && l <= r && r <= n);
    }

    void validatePoint(Int v) const {
        assert(0 <= v && v < n);
    }

    Int build(Int l, Int r, Int &nextId) {
        Int idx = nodes.size();
        nodes.push_back({l, r, -1, -1, l, l});
        if (r - l == 1)
            return idx;
        Int inId = nextId++, outId = nextId++;
        nodes[idx].inId = inId;
        nodes[idx].outId = outId;
        Int m = (l + r) >> 1, left = build(l, m, nextId), right = build(m, r, nextId);
        nodes[idx].left = left;
        nodes[idx].right = right;
        G.add_edge(inId, nodes[left].inId, zero);
        G.add_edge(inId, nodes[right].inId, zero);
        G.add_edge(nodes[left].outId, outId, zero);
        G.add_edge(nodes[right].outId, outId, zero);
        return idx;
    }

    void connectRangeToVertex(Int idx, Int l, Int r, Int to, T cost) {
        auto node = nodes[idx];
        if (r <= node.l || node.r <= l)
            return;
        if (l <= node.l && node.r <= r) {
            G.add_edge(node.outId, to, cost);
            return;
        }
        connectRangeToVertex(node.left, l, r, to, cost);
        connectRangeToVertex(node.right, l, r, to, cost);
    }

    void connectVertexToRange(Int src, Int idx, Int l, Int r, T cost) {
        auto node = nodes[idx];
        if (r <= node.l || node.r <= l)
            return;
        if (l <= node.l && node.r <= r) {
            G.add_edge(src, node.inId, cost);
            return;
        }
        connectVertexToRange(src, node.left, l, r, cost);
        connectVertexToRange(src, node.right, l, r, cost);
    }

public:
    WeightedDirectedGraph<T> G;
    Int n, baseLen;

    // 元頂点をin/out木の葉として共有する。O(N)、頂点数3N-2。
    explicit Range_Edge_Graph(Int N, T zero = T{})
        : zero(zero), G(N == 0 ? 0 : 3 * N - 2), n(N), baseLen(N == 0 ? 0 : 3 * N - 2) {
        assert(N >= 0);
        if (N) {
            Int nextId = N;
            root = build(0, N, nextId);
            assert(nextId == baseLen);
        }
    }

    // 内部の重み付き有向グラフへの参照を返す。O(1)。
    auto &graph() {
        return G;
    }

    // 内部の重み付き有向グラフへの読み取り専用参照を返す。O(1)。
    const auto &graph() const {
        return G;
    }

    // 補助頂点を含む現在の頂点数を返す。O(1)。
    Int len() const {
        return G.len;
    }

    // 補助頂点を含む現在の頂点数を返す。O(1)。
    Int size() const {
        return len();
    }

    // 元の頂点数を返す。O(1)。
    Int original_len() const {
        return n;
    }

    // 範囲辺を追加する前の補助頂点を含む頂点数を返す。O(1)。
    Int base_len() const {
        return baseLen;
    }

    // 範囲同士は中継点を2個追加する。各範囲追加O(log N)、点同士O(1)。
    void add_range_to_range_edge(Int fl, Int fr, Int tl, Int tr, T cost) {
        validateRange(fl, fr);
        validateRange(tl, tr);
        if (fl == fr || tl == tr)
            return;
        Int from = G.len, to = G.len + 1;
        G.len += 2;
        G.edges.resize(G.len);
        connectRangeToVertex(root, fl, fr, from, zero);
        G.add_edge(from, to, cost);
        connectVertexToRange(to, root, tl, tr, zero);
    }

    // 点から半開区間への辺を追加する。償却O(log N)。
    void add_point_to_range_edge(Int src, Int tl, Int tr, T cost) {
        validatePoint(src);
        validateRange(tl, tr);
        if (tl != tr)
            connectVertexToRange(src, root, tl, tr, cost);
    }

    // 半開区間から点への辺を追加する。償却O(log N)。
    void add_range_to_point_edge(Int fl, Int fr, Int to, T cost) {
        validateRange(fl, fr);
        validatePoint(to);
        if (fl != fr)
            connectRangeToVertex(root, fl, fr, to, cost);
    }

    // 点から点への辺を追加する。償却O(1)。
    void add_point_to_point_edge(Int src, Int to, T cost) {
        validatePoint(src);
        validatePoint(to);
        G.add_edge(src, to, cost);
    }

    // 半開区間同士を結ぶ辺を追加する。償却O(log N)。
    void add_edge(Int fl, Int fr, Int tl, Int tr, T cost) {
        add_range_to_range_edge(fl, fr, tl, tr, cost);
    }

    // 点から点への辺を追加する。償却O(1)。
    void add_edge(Int src, Int to, T cost) {
        add_point_to_point_edge(src, to, cost);
    }
};

// n頂点から範囲辺グラフを構築する。O(n)。
template <class T = Int> auto initWeightedRangeGraph(Int n, T zero = T{}) {
    return Range_Edge_Graph<T>(n, zero);
}

template <class T> const auto &graph(const Range_Edge_Graph<T> &g) {
    return g.graph();
}

template <class T> auto &graph(Range_Edge_Graph<T> &g) {
    return g.graph();
}

// 補助頂点を含む現在の頂点数を返す。O(1)。
template <class T> Int len(const Range_Edge_Graph<T> &g) {
    return g.len();
}

// 元の頂点数を返す。O(1)。
template <class T> Int original_len(const Range_Edge_Graph<T> &g) {
    return g.original_len();
}

// 範囲辺追加前の補助頂点を含む頂点数を返す。O(1)。
template <class T> Int base_len(const Range_Edge_Graph<T> &g) {
    return g.base_len();
}

// 半開区間同士を結ぶ辺を追加する。償却O(log N)。
template <class T>
void add_range_to_range_edge(Range_Edge_Graph<T> &g, Int fl, Int fr, Int tl, Int tr, T cost) {
    g.add_range_to_range_edge(fl, fr, tl, tr, cost);
}

// 点から半開区間への辺を追加する。償却O(log N)。
template <class T>
void add_point_to_range_edge(Range_Edge_Graph<T> &g, Int src, Int tl, Int tr, T cost) {
    g.add_point_to_range_edge(src, tl, tr, cost);
}

// 半開区間から点への辺を追加する。償却O(log N)。
template <class T>
void add_range_to_point_edge(Range_Edge_Graph<T> &g, Int fl, Int fr, Int to, T cost) {
    g.add_range_to_point_edge(fl, fr, to, cost);
}

// 点から点への辺を追加する。償却O(1)。
template <class T> void add_point_to_point_edge(Range_Edge_Graph<T> &g, Int src, Int to, T cost) {
    g.add_point_to_point_edge(src, to, cost);
}
}
