#pragma once
#include <cplib/common.hpp>
#include <unordered_map>

namespace cplib {
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

    explicit BasicGraph(Int n = 0, Int capacity = 0) : len(n) {
        assert(n >= 0);
        if constexpr (Static) {
            capacity *= Directed ? 1 : 2;
            src.reserve(capacity);
            dst.reserve(capacity);
            cost.reserve(capacity);
        } else
            edges.resize(n);
    }

    void add_edge_dynamic_impl(Int u, Int v, T w, bool directed) {
        edges[u].emplace_back(v, w);
        if (!directed)
            edges[v].emplace_back(u, w);
    }

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

    void add_edge(Int u, Int v, T w)
        requires Weighted
    {
        if constexpr (Static)
            add_edge_static_impl(u, v, w, Directed);
        else
            add_edge_dynamic_impl(u, v, w, Directed);
    }

    void add_edge(Int u, Int v)
        requires(!Weighted)
    {
        if constexpr (Static)
            add_edge_static_impl(u, v, T(1), Directed);
        else
            add_edge_dynamic_impl(u, v, T(1), Directed);
    }

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

    void build()
        requires Static
    {
        build_impl();
    }

    void static_graph_initialized_check() const
        requires Static
    {
        assert(!start.empty());
    }

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

            auto operator*() const {
                const auto &e = edges[i];
                if constexpr (Weighted)
                    return std::pair<Int, T>{e.first, e.second};
                else
                    return Int(e.first);
            }

            Iterator &operator++() {
                ++i;
                return *this;
            }

            bool operator!=(const Iterator &b) const {
                return i != b.i;
            }
        };

        Iterator begin() const {
            return {edges, 0};
        }

        Iterator end() const {
            return {edges, edges.size()};
        }
    };

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

template <class T = Int> inline auto initWeightedDirectedGraph(Int n) {
    return WeightedDirectedGraph<T>(n);
}

using UnWeightedDirectedGraph = BasicGraph<Int, false, true, false>;

inline auto initUnWeightedDirectedGraph(Int n) {
    return UnWeightedDirectedGraph(n);
}
template <class T = Int> using WeightedUnDirectedGraph = BasicGraph<T, false, false, true>;

template <class T = Int> inline auto initWeightedUnDirectedGraph(Int n) {
    return WeightedUnDirectedGraph<T>(n);
}

using UnWeightedUnDirectedGraph = BasicGraph<Int, false, false, false>;

inline auto initUnWeightedUnDirectedGraph(Int n) {
    return UnWeightedUnDirectedGraph(n);
}
template <class T = Int> using WeightedDirectedStaticGraph = BasicGraph<T, true, true, true>;

template <class T = Int> inline auto initWeightedDirectedStaticGraph(Int n, Int capacity = 0) {
    return WeightedDirectedStaticGraph<T>(n, capacity);
}

using UnWeightedDirectedStaticGraph = BasicGraph<Int, true, true, false>;

inline auto initUnWeightedDirectedStaticGraph(Int n, Int capacity = 0) {
    return UnWeightedDirectedStaticGraph(n, capacity);
}
template <class T = Int> using WeightedUnDirectedStaticGraph = BasicGraph<T, true, false, true>;

template <class T = Int> inline auto initWeightedUnDirectedStaticGraph(Int n, Int capacity = 0) {
    return WeightedUnDirectedStaticGraph<T>(n, capacity);
}

using UnWeightedUnDirectedStaticGraph = BasicGraph<Int, true, false, false>;

inline auto initUnWeightedUnDirectedStaticGraph(Int n, Int capacity = 0) {
    return UnWeightedUnDirectedStaticGraph(n, capacity);
}

template <GraphTypes G> Int len(const G &g) {
    return g.len;
}

template <GraphTypes G, class... A> void add_edge(G &g, A &&...a) {
    g.add_edge(std::forward<A>(a)...);
}

template <StaticGraphTypes G> void build(G &g) {
    g.build();
}

template <StaticGraphTypes G> void build_impl(G &g) {
    g.build_impl();
}

template <StaticGraphTypes G> void static_graph_initialized_check(const G &g) {
    g.static_graph_initialized_check();
}

template <GraphTypes G> auto to_and_cost(const G &g, Int x) {
    return g.to_and_cost(x);
}

template <DynamicGraphTypes G>
void add_edge_dynamic_impl(G &g, Int u, Int v, typename G::value_type cost, bool directed) {
    g.add_edge_dynamic_impl(u, v, cost, directed);
}

template <StaticGraphTypes G>
void add_edge_static_impl(G &g, Int u, Int v, typename G::value_type cost, bool directed) {
    g.add_edge_static_impl(u, v, cost, directed);
}

template <class K, class T, bool Directed, bool Weighted> struct TableGraph {
    std::unordered_map<K, Int> toi;
    std::vector<K> v;
    BasicGraph<T, false, Directed, Weighted> graph;

    explicit TableGraph(std::span<const K> input)
        : v(input.begin(), input.end()), graph(input.size()) {
        for (Int i = 0; i < Int(v.size()); ++i)
            toi[v[i]] = i;
    }

    void add_edge(const K &a, const K &b)
        requires(!Weighted)
    {
        graph.add_edge(toi.at(a), toi.at(b));
    }

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

            auto operator*() const {
                const auto &e = edges[i];
                if constexpr (Weighted)
                    return std::pair<K, T>{owner->v[e.first], e.second};
                else
                    return owner->v[e.first];
            }

            Iterator &operator++() {
                ++i;
                return *this;
            }

            bool operator!=(const Iterator &b) const {
                return i != b.i;
            }
        };

        Iterator begin() const {
            return {owner, edges, 0};
        }

        Iterator end() const {
            return {owner, edges, edges.size()};
        }
    };

    auto operator[](const K &key) const {
        return Adjacency{this, graph.to_and_cost(toi.at(key))};
    }
};
template <class K, class T = Int> using WeightedDirectedTableGraph = TableGraph<K, T, true, true>;

template <class T = Int, class K> auto initWeightedDirectedTableGraph(std::span<const K> v) {
    return WeightedDirectedTableGraph<K, T>(v);
}

template <class T = Int, class K> auto initWeightedDirectedTableGraph(const std::vector<K> &v) {
    return WeightedDirectedTableGraph<K, T>(v);
}
template <class K> using UnWeightedDirectedTableGraph = TableGraph<K, Int, true, false>;

template <class K> auto initUnWeightedDirectedTableGraph(std::span<const K> v) {
    return UnWeightedDirectedTableGraph<K>(v);
}

template <class K> auto initUnWeightedDirectedTableGraph(const std::vector<K> &v) {
    return UnWeightedDirectedTableGraph<K>(v);
}
template <class K, class T = Int>
using WeightedUnDirectedTableGraph = TableGraph<K, T, false, true>;

template <class T = Int, class K> auto initWeightedUnDirectedTableGraph(std::span<const K> v) {
    return WeightedUnDirectedTableGraph<K, T>(v);
}

template <class T = Int, class K> auto initWeightedUnDirectedTableGraph(const std::vector<K> &v) {
    return WeightedUnDirectedTableGraph<K, T>(v);
}
template <class K> using UnWeightedUnDirectedTableGraph = TableGraph<K, Int, false, false>;

template <class K> auto initUnWeightedUnDirectedTableGraph(std::span<const K> v) {
    return UnWeightedUnDirectedTableGraph<K>(v);
}

template <class K> auto initUnWeightedUnDirectedTableGraph(const std::vector<K> &v) {
    return UnWeightedUnDirectedTableGraph<K>(v);
}

template <class K, class T, bool D, bool W, class... A>
void add_edge(TableGraph<K, T, D, W> &g, A &&...a) {
    g.add_edge(std::forward<A>(a)...);
}
}
