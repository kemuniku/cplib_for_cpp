#pragma once
#include <cplib/common.hpp>
#include <unordered_map>

namespace cplib::legacy_graph {
// 旧版は辺番号を保持せず、無向辺の両方向を挿入時に格納する。
template <class T, bool Static, bool Directed, bool Weighted> struct BasicGraph {
    using value_type = T;
    using Edge = std::pair<std::int32_t, T>;
    static constexpr bool is_static = Static, is_directed = Directed, is_weighted = Weighted;
    Int len;
    std::vector<std::vector<Edge>> edges;
    std::vector<std::int32_t> src, dst, start;
    std::vector<T> cost;
    std::vector<Edge> elist;

    // n頂点のグラフを作る。動的版O(n)、静的版O(1)。
    explicit BasicGraph(Int n = 0, Int capacity = 0) : len{n} {
        assert(n >= 0);
        if constexpr (Static) {
            capacity *= Directed ? 1 : 2;
            src.reserve(capacity);
            dst.reserve(capacity);
            cost.reserve(capacity);
        } else
            edges.resize(n);
    }

    // 頂点数を返す。辺数とは異なる。O(1)。
    Int size() const {
        return len;
    }

    // 向きを指定して動的グラフに辺を追加する。償却O(1)。
    void add_edge_dynamic_impl(Int u, Int v, T w, bool directed) {
        edges[u].emplace_back(v, w);
        if (!directed)
            edges[v].emplace_back(u, w);
    }

    // 向きを指定して静的グラフに辺を追加する。償却O(1)。
    void add_edge_static_impl(Int u, Int v, T w, bool directed) {
        src.push_back(u);
        dst.push_back(v);
        cost.push_back(w);
        if (!directed) {
            src.push_back(v);
            dst.push_back(u);
            cost.push_back(w);
        }
    }

    // 重み付き辺を追加する。償却O(1)。静的版は再buildが必要。
    void add_edge(Int u, Int v, T w)
        requires Weighted
    {
        if constexpr (Static)
            add_edge_static_impl(u, v, w, Directed);
        else
            add_edge_dynamic_impl(u, v, w, Directed);
    }

    // 重みなし辺を追加する。償却O(1)。静的版は再buildが必要。
    void add_edge(Int u, Int v)
        requires(!Weighted)
    {
        if constexpr (Static)
            add_edge_static_impl(u, v, T(1), Directed);
        else
            add_edge_dynamic_impl(u, v, T(1), Directed);
    }

    // 静的グラフの隣接配列を構築する。O(V+E)。
    void build_impl()
        requires Static
    {
        start.assign(len + 1, 0);
        for (auto u : src)
            ++start[u];
        std::partial_sum(start.begin(), start.end(), start.begin());
        elist.resize(start.back());
        for (Int i = Int(src.size()) - 1; i >= 0; --i)
            elist[--start[src[i]]] = {dst[i], cost[i]};
    }

    // 静的グラフの隣接配列を構築する。O(V+E)。
    void build()
        requires Static
    {
        build_impl();
    }

    // 静的グラフが構築済みか検査する。O(1)。
    void static_graph_initialized_check() const
        requires Static
    {
        assert(!start.empty());
    }

    // 頂点xの隣接頂点と重みのビューを返す。取得O(1)、列挙O(deg(x))。
    std::span<const Edge> to_and_cost(Int x) const {
        if constexpr (Static) {
            static_graph_initialized_check();
            return std::span<const Edge>(elist).subspan(start[x], start[x + 1] - start[x]);
        } else
            return edges[x];
    }

    struct Adjacency {
        std::span<const Edge> edges;

        struct Iterator {
            std::span<const Edge> edges;
            std::size_t i;

            // 現在の隣接要素を返す。値のコピーがO(1)ならO(1)。
            auto operator*() const {
                const auto &e = edges[i];
                if constexpr (Weighted)
                    return std::pair<Int, T>{e.first, e.second};
                else
                    return Int(e.first);
            }

            // 次の隣接要素へ進む。O(1)。
            Iterator &operator++() {
                ++i;
                return *this;
            }

            // 二つの反復子が異なる位置を指すか返す。O(1)。
            bool operator!=(const Iterator &b) const {
                return i != b.i;
            }
        };

        // 隣接要素の列挙の開始位置を返す。O(1)。
        Iterator begin() const {
            return {edges, 0};
        }

        // 隣接要素の列挙の終了位置を返す。O(1)。
        Iterator end() const {
            return {edges, edges.size()};
        }
    };

    // 隣接頂点（重み付きでは頂点と重み）のビューを返す。取得O(1)、列挙O(deg(x))。
    auto operator[](Int x) const {
        return Adjacency{to_and_cost(x)};
    }
};
template <class G>
concept GraphTypes = requires {
    G::is_static;
    G::is_directed;
    G::is_weighted;
};
template <class G>
concept DirectedGraph = GraphTypes<G> && G::is_directed;
template <class G>
concept UnDirectedGraph = GraphTypes<G> && (!G::is_directed);
template <class G>
concept WeightedGraph = GraphTypes<G> && G::is_weighted;
template <class G>
concept UnWeightedGraph = GraphTypes<G> && (!G::is_weighted);
template <class G>
concept DynamicGraphTypes = GraphTypes<G> && (!G::is_static);
template <class G>
concept StaticGraphTypes = GraphTypes<G> && G::is_static;
template <class T = Int> using WeightedDirectedGraph = BasicGraph<T, false, true, true>;

// n頂点の動的有向重み付きグラフを作る。O(n)。
template <class T = Int> inline auto initWeightedDirectedGraph(Int n) {
    return WeightedDirectedGraph<T>(n);
}

using UnWeightedDirectedGraph = BasicGraph<Int, false, true, false>;

// n頂点の動的有向重みなしグラフを作る。O(n)。
inline auto initUnWeightedDirectedGraph(Int n) {
    return UnWeightedDirectedGraph(n);
}
template <class T = Int> using WeightedUnDirectedGraph = BasicGraph<T, false, false, true>;

// n頂点の動的無向重み付きグラフを作る。O(n)。
template <class T = Int> inline auto initWeightedUnDirectedGraph(Int n) {
    return WeightedUnDirectedGraph<T>(n);
}

using UnWeightedUnDirectedGraph = BasicGraph<Int, false, false, false>;

// n頂点の動的無向重みなしグラフを作る。O(n)。
inline auto initUnWeightedUnDirectedGraph(Int n) {
    return UnWeightedUnDirectedGraph(n);
}
template <class T = Int> using WeightedDirectedStaticGraph = BasicGraph<T, true, true, true>;

// n頂点の静的有向重み付きグラフを作る。O(1)。
template <class T = Int> inline auto initWeightedDirectedStaticGraph(Int n, Int capacity = 0) {
    return WeightedDirectedStaticGraph<T>(n, capacity);
}

using UnWeightedDirectedStaticGraph = BasicGraph<Int, true, true, false>;

// n頂点の静的有向重みなしグラフを作る。O(1)。
inline auto initUnWeightedDirectedStaticGraph(Int n, Int capacity = 0) {
    return UnWeightedDirectedStaticGraph(n, capacity);
}
template <class T = Int> using WeightedUnDirectedStaticGraph = BasicGraph<T, true, false, true>;

// n頂点の静的無向重み付きグラフを作る。O(1)。
template <class T = Int> inline auto initWeightedUnDirectedStaticGraph(Int n, Int capacity = 0) {
    return WeightedUnDirectedStaticGraph<T>(n, capacity);
}

using UnWeightedUnDirectedStaticGraph = BasicGraph<Int, true, false, false>;

// n頂点の静的無向重みなしグラフを作る。O(1)。
inline auto initUnWeightedUnDirectedStaticGraph(Int n, Int capacity = 0) {
    return UnWeightedUnDirectedStaticGraph(n, capacity);
}

// グラフに辺を追加する。償却O(1)、ラベル付きグラフでは期待償却O(1)。
template <GraphTypes G, class... A> void add_edge(G &g, A &&...a) {
    g.add_edge(std::forward<A>(a)...);
}

// 静的グラフの隣接配列を構築する。O(V+E)。
template <StaticGraphTypes G> void build(G &g) {
    g.build();
}

// 静的グラフの隣接配列を構築する。O(V+E)。
template <StaticGraphTypes G> void build_impl(G &g) {
    g.build_impl();
}

// 静的グラフが構築済みか検査する。O(1)。
template <StaticGraphTypes G> void static_graph_initialized_check(const G &g) {
    g.static_graph_initialized_check();
}

// 隣接頂点と重みのビューを返す。取得O(1)、列挙O(deg(x))。
template <GraphTypes G> auto to_and_cost(const G &g, Int x) {
    return g.to_and_cost(x);
}

// 向きを指定して動的グラフに辺を追加する。償却O(1)。
template <DynamicGraphTypes G>
void add_edge_dynamic_impl(G &g, Int u, Int v, typename G::value_type cost, bool directed) {
    g.add_edge_dynamic_impl(u, v, cost, directed);
}

// 向きを指定して静的グラフに辺を追加する。償却O(1)。
template <StaticGraphTypes G>
void add_edge_static_impl(G &g, Int u, Int v, typename G::value_type cost, bool directed) {
    g.add_edge_static_impl(u, v, cost, directed);
}

template <class K, class T, bool Directed, bool Weighted> struct TableGraph {
    std::unordered_map<K, Int> toi;
    std::vector<K> v;
    BasicGraph<T, false, Directed, Weighted> graph;

    // ラベル付きグラフを構築する。頂点数Vに対して期待O(V)、最悪O(V²)。
    explicit TableGraph(std::span<const K> input)
        : v(input.begin(), input.end()), graph(input.size()) {
        for (Int i = 0; i < Int(v.size()); ++i)
            toi[v[i]] = i;
    }

    // ラベルで指定した辺を追加する。期待償却O(1)。
    void add_edge(const K &a, const K &b)
        requires(!Weighted)
    {
        graph.add_edge(toi.at(a), toi.at(b));
    }

    // ラベルと重みを指定して辺を追加する。期待償却O(1)。
    void add_edge(const K &a, const K &b, T cost)
        requires Weighted
    {
        graph.add_edge(toi.at(a), toi.at(b), cost);
    }

    struct Adjacency {
        const TableGraph *owner;
        std::span<const typename decltype(graph)::Edge> edges;

        struct Iterator {
            const TableGraph *owner;
            std::span<const typename decltype(graph)::Edge> edges;
            std::size_t i;

            // 現在の隣接要素を返す。値のコピーがO(1)ならO(1)。
            auto operator*() const {
                const auto &e = edges[i];
                if constexpr (Weighted)
                    return std::pair<K, T>{owner->v[e.first], e.second};
                else
                    return owner->v[e.first];
            }

            // 次の隣接要素へ進む。O(1)。
            Iterator &operator++() {
                ++i;
                return *this;
            }

            // 二つの反復子が異なる位置を指すか返す。O(1)。
            bool operator!=(const Iterator &b) const {
                return i != b.i;
            }
        };

        // 隣接要素の列挙の開始位置を返す。O(1)。
        Iterator begin() const {
            return {owner, edges, 0};
        }

        // 隣接要素の列挙の終了位置を返す。O(1)。
        Iterator end() const {
            return {owner, edges, edges.size()};
        }
    };

    // ラベルkeyの隣接頂点のビューを返す。取得は期待O(1)、列挙O(deg(key))。
    auto operator[](const K &key) const {
        return Adjacency{this, graph.to_and_cost(toi.at(key))};
    }
};
template <class K, class T = Int> using WeightedDirectedTableGraph = TableGraph<K, T, true, true>;

// ラベル付き有向重み付きグラフを作る。頂点数Vに対して期待O(V)、最悪O(V²)。
template <class T = Int, class K> auto initWeightedDirectedTableGraph(std::span<const K> v) {
    return WeightedDirectedTableGraph<K, T>(v);
}

// ラベル付き有向重み付きグラフを作る。頂点数Vに対して期待O(V)、最悪O(V²)。
template <class T = Int, class K> auto initWeightedDirectedTableGraph(const std::vector<K> &v) {
    return WeightedDirectedTableGraph<K, T>(v);
}
template <class K> using UnWeightedDirectedTableGraph = TableGraph<K, Int, true, false>;

// ラベル付き有向重みなしグラフを作る。頂点数Vに対して期待O(V)、最悪O(V²)。
template <class K> auto initUnWeightedDirectedTableGraph(std::span<const K> v) {
    return UnWeightedDirectedTableGraph<K>(v);
}

// ラベル付き有向重みなしグラフを作る。頂点数Vに対して期待O(V)、最悪O(V²)。
template <class K> auto initUnWeightedDirectedTableGraph(const std::vector<K> &v) {
    return UnWeightedDirectedTableGraph<K>(v);
}
template <class K, class T = Int>
using WeightedUnDirectedTableGraph = TableGraph<K, T, false, true>;

// ラベル付き無向重み付きグラフを作る。頂点数Vに対して期待O(V)、最悪O(V²)。
template <class T = Int, class K> auto initWeightedUnDirectedTableGraph(std::span<const K> v) {
    return WeightedUnDirectedTableGraph<K, T>(v);
}

// ラベル付き無向重み付きグラフを作る。頂点数Vに対して期待O(V)、最悪O(V²)。
template <class T = Int, class K> auto initWeightedUnDirectedTableGraph(const std::vector<K> &v) {
    return WeightedUnDirectedTableGraph<K, T>(v);
}
template <class K> using UnWeightedUnDirectedTableGraph = TableGraph<K, Int, false, false>;

// ラベル付き無向重みなしグラフを作る。頂点数Vに対して期待O(V)、最悪O(V²)。
template <class K> auto initUnWeightedUnDirectedTableGraph(std::span<const K> v) {
    return UnWeightedUnDirectedTableGraph<K>(v);
}

// ラベル付き無向重みなしグラフを作る。頂点数Vに対して期待O(V)、最悪O(V²)。
template <class K> auto initUnWeightedUnDirectedTableGraph(const std::vector<K> &v) {
    return UnWeightedUnDirectedTableGraph<K>(v);
}

// グラフに辺を追加する。償却O(1)、ラベル付きグラフでは期待償却O(1)。
template <class K, class T, bool D, bool W, class... A>
void add_edge(TableGraph<K, T, D, W> &g, A &&...a) {
    g.add_edge(std::forward<A>(a)...);
}
}
