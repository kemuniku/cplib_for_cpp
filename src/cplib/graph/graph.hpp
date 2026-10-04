#pragma once
#include <cplib/common.hpp>
#include <cplib/iterator_views.hpp>
#include <tuple>
#include <unordered_map>

namespace cplib {
template <class T> struct EdgeInfo {
    Int src, dst;
    T cost;
};

template <> struct EdgeInfo<void> {
    Int src, dst;
};

struct AdjacentEdge {
    std::int32_t dst, id;
};

template <class T> struct WeightedAdjacentEdge {
    std::int32_t dst, id;
    T cost;
};
template <class T>
using GraphAdjacent = std::conditional_t<std::is_void_v<T>, AdjacentEdge, WeightedAdjacentEdge<T>>;

// 重み付きグラフでは走査の高速化のため、重みを隣接配列にも保存する。公開配列の重みを直接変更しないこと。
template <class T, bool Static = false, bool Directed = true> class BasicGraph {
    int direction = 0;
    std::vector<bool> directed;

    // 通常は方向だけを保持し、混在時だけ辺ごとの方向を保存する。償却O(1)。
    void record_direction(bool dir) {
        int d = dir ? 1 : 2;
        if (!direction)
            direction = d;
        else if (direction != d) {
            if (direction != 3) {
                directed.assign(edge_info.size() - 1, direction == 1);
                direction = 3;
            }
            directed.push_back(dir);
        }
    }

public:
    using weight_type = T;
    using cost_type = std::conditional_t<std::is_void_v<T>, Int, T>;
    using adjacent_type = GraphAdjacent<T>;
    static constexpr bool is_static = Static, is_directed = Directed;
    Int len = 0;
    std::vector<std::vector<adjacent_type>> edges;
    std::vector<EdgeInfo<T>> edge_info;
    std::vector<adjacent_type> elist;
    std::vector<std::int32_t> start;

    // 頂点数nのグラフを構築する。動的版O(n)、静的版O(1)。
    explicit BasicGraph(Int n = 0, Int capacity = 0) : len(n) {
        edge_info.reserve(capacity);
        if constexpr (!Static)
            edges.resize(n);
    }

    // 辺を追加し、辺番号を返す。償却O(1)。静的版は再buildが必要。
    Int add_edge_impl(Int u, Int v, cost_type cost, bool dir) {
        Int id = edge_info.size();
        if constexpr (std::is_void_v<T>)
            edge_info.push_back({u, v});
        else
            edge_info.push_back({u, v, cost});
        if constexpr (Static) {
            record_direction(dir);
            start.clear();
        } else {
            if constexpr (std::is_void_v<T>) {
                edges[u].push_back({std::int32_t(v), std::int32_t(id)});
                if (!dir)
                    edges[v].push_back({std::int32_t(u), std::int32_t(id)});
            } else {
                edges[u].push_back({std::int32_t(v), std::int32_t(id), cost});
                if (!dir)
                    edges[v].push_back({std::int32_t(u), std::int32_t(id), cost});
            }
        }
        return id;
    }

    Int add_edge(Int u, Int v, cost_type cost = 1) {
        return add_edge_impl(u, v, cost, Directed);
    }

    // 隣接配列をCSRとして構築する。O(V+E)、追加順を保持する。
    void build()
        requires Static
    {
        start.assign(len + 1, 0);
        auto undirected = [&](Int id) {
            return direction == 2 || (direction == 3 && !directed[id]);
        };
        for (Int id = 0; id < Int(edge_info.size()); ++id) {
            auto e = edge_info[id];
            ++start[e.src + 1];
            if (undirected(id))
                ++start[e.dst + 1];
        }
        for (Int i = 0; i < len; ++i)
            start[i + 1] += start[i];
        elist.resize(start.back());
        auto cursor = start;
        for (Int id = 0; id < Int(edge_info.size()); ++id) {
            auto e = edge_info[id];
            if constexpr (std::is_void_v<T>) {
                elist[cursor[e.src]++] = {std::int32_t(e.dst), std::int32_t(id)};
                if (undirected(id))
                    elist[cursor[e.dst]++] = {std::int32_t(e.src), std::int32_t(id)};
            } else {
                elist[cursor[e.src]++] = {std::int32_t(e.dst), std::int32_t(id), e.cost};
                if (undirected(id))
                    elist[cursor[e.dst]++] = {std::int32_t(e.src), std::int32_t(id), e.cost};
            }
        }
    }

    void static_graph_initialized_check() const
        requires Static
    {
        assert(!start.empty());
    }

    // 辺と隣接要素の容量を予約する。既存要素を維持する。
    void reserve(Int capacity, std::span<const Int> degrees = {})
        requires(!Static)
    {
        assert(capacity >= 0);
        assert(degrees.empty() || Int(degrees.size()) == len);
        edge_info.reserve(capacity);
        for (std::size_t u = 0; u < degrees.size(); ++u) {
            assert(degrees[u] >= 0);
            edges[u].reserve(degrees[u]);
        }
    }

    // 辺数を返す。無向辺も一辺として数える。O(1)。
    Int edge_count() const {
        return edge_info.size();
    }

    // 辺番号から追加時の向きで、重みなしは (src, dst)、重みありは (src, dst, cost) を返す。O(1)。
    auto get_edge(Int id) const {
        auto e = edge_info[id];
        if constexpr (std::is_void_v<T>)
            return std::pair{e.src, e.dst};
        else
            return std::tuple{e.src, e.dst, e.cost};
    }

    // 隣接辺の参照範囲を返す。O(1)、走査O(deg(v))、追加領域O(1)。
    auto adjacency(Int v) const {
        if constexpr (Static) {
            static_graph_initialized_check();
            return std::span<const adjacent_type>(elist).subspan(start[v], start[v + 1] - start[v]);
        } else
            return std::span<const adjacent_type>(edges[v]);
    }

    // 隣接頂点と辺番号を追加順に列挙する。O(deg(v))。
    auto to_and_id(Int v) const {
        return detail::transform_range(adjacency(v),
                                       [](auto e) { return std::pair<Int, Int>{e.dst, e.id}; });
    }

    // 隣接頂点と重みを列挙する。重みなしは重み1を返す。O(deg(v))。
    auto to_and_cost(Int v) const {
        return detail::transform_range(adjacency(v), [](auto e) {
            if constexpr (std::is_void_v<T>)
                return std::pair<Int, Int>{e.dst, 1};
            else
                return std::pair<Int, T>{e.dst, e.cost};
        });
    }

    // 隣接頂点、重み、辺番号を列挙する。重みなしは重み1を返す。O(deg(v))。
    auto to_and_cost_and_id(Int v) const {
        return detail::transform_range(adjacency(v), [](auto e) {
            if constexpr (std::is_void_v<T>)
                return std::tuple<Int, Int, Int>{e.dst, 1, e.id};
            else
                return std::tuple<Int, T, Int>{e.dst, e.cost, e.id};
        });
    }

    auto operator[](Int v) const {
        if constexpr (std::is_void_v<T>)
            return detail::transform_range(adjacency(v), [](auto e) { return Int(e.dst); });
        else
            return to_and_cost(v);
    }
};
template <class T> using DynamicGraph = BasicGraph<T, false>;
template <class T> using StaticGraph = BasicGraph<T, true>;
template <class T = Int> using WeightedDirectedGraph = BasicGraph<T, false, true>;
template <class T = Int> using WeightedUnDirectedGraph = BasicGraph<T, false, false>;
using UnWeightedDirectedGraph = BasicGraph<void, false, true>;
using UnWeightedUnDirectedGraph = BasicGraph<void, false, false>;
template <class T = Int> using WeightedDirectedStaticGraph = BasicGraph<T, true, true>;
template <class T = Int> using WeightedUnDirectedStaticGraph = BasicGraph<T, true, false>;
using UnWeightedDirectedStaticGraph = BasicGraph<void, true, true>;
using UnWeightedUnDirectedStaticGraph = BasicGraph<void, true, false>;
template <class G>
concept GraphTypes = requires(const G &g) {
    typename G::cost_type;
    g.len;
    g.to_and_cost(Int(0));
};
template <class G>
concept DirectedGraph = GraphTypes<G> && G::is_directed;
template <class G>
concept UnDirectedGraph = GraphTypes<G> && !G::is_directed;
template <class G>
concept WeightedGraph = GraphTypes<G> && !std::is_void_v<typename G::weight_type>;
template <class G>
concept UnWeightedGraph = GraphTypes<G> && std::is_void_v<typename G::weight_type>;
template <class G>
concept DynamicGraphTypes = GraphTypes<G> && !G::is_static;
template <class G>
concept StaticGraphTypes = GraphTypes<G> && G::is_static;

template <class T = Int> auto initWeightedDirectedGraph(Int n, Int capacity = 0) {
    return WeightedDirectedGraph<T>(n, capacity);
}

template <class T = Int> auto initWeightedUnDirectedGraph(Int n, Int capacity = 0) {
    return WeightedUnDirectedGraph<T>(n, capacity);
}

inline auto initUnWeightedDirectedGraph(Int n, Int capacity = 0) {
    return UnWeightedDirectedGraph(n, capacity);
}

inline auto initUnWeightedUnDirectedGraph(Int n, Int capacity = 0) {
    return UnWeightedUnDirectedGraph(n, capacity);
}

template <class T = Int> auto initWeightedDirectedStaticGraph(Int n, Int capacity = 0) {
    return WeightedDirectedStaticGraph<T>(n, capacity);
}

template <class T = Int> auto initWeightedUnDirectedStaticGraph(Int n, Int capacity = 0) {
    return WeightedUnDirectedStaticGraph<T>(n, capacity);
}

inline auto initUnWeightedDirectedStaticGraph(Int n, Int capacity = 0) {
    return UnWeightedDirectedStaticGraph(n, capacity);
}

inline auto initUnWeightedUnDirectedStaticGraph(Int n, Int capacity = 0) {
    return UnWeightedUnDirectedStaticGraph(n, capacity);
}

template <GraphTypes G> Int len(const G &g) {
    return g.len;
}

template <GraphTypes G> Int edge_count(const G &g) {
    return g.edge_count();
}

template <GraphTypes G> auto get_edge(const G &g, Int id) {
    return g.get_edge(id);
}

template <StaticGraphTypes G> void build(G &g) {
    g.build();
}

template <StaticGraphTypes G> void build_impl(G &g) {
    g.build();
}

template <StaticGraphTypes G> void static_graph_initialized_check(const G &g) {
    g.static_graph_initialized_check();
}

template <DynamicGraphTypes G> void reserve(G &g, Int capacity, std::span<const Int> degrees = {}) {
    g.reserve(capacity, degrees);
}

template <DynamicGraphTypes G>
Int add_edge_dynamic_impl(G &g, Int u, Int v, typename G::cost_type cost, bool directed) {
    return g.add_edge_impl(u, v, cost, directed);
}

template <DynamicGraphTypes G>
    requires UnWeightedGraph<G>
Int add_edge_dynamic_impl(G &g, Int u, Int v, bool directed) {
    return g.add_edge_impl(u, v, 1, directed);
}

template <StaticGraphTypes G>
Int add_edge_static_impl(G &g, Int u, Int v, typename G::cost_type cost, bool directed) {
    return g.add_edge_impl(u, v, cost, directed);
}

template <StaticGraphTypes G>
    requires UnWeightedGraph<G>
Int add_edge_static_impl(G &g, Int u, Int v, bool directed) {
    return g.add_edge_impl(u, v, 1, directed);
}

template <class Label, class Weight = void, bool Directed = true> struct TableGraph {
    std::unordered_map<Label, Int> toi;
    std::vector<Label> v;
    BasicGraph<Weight, false, Directed> graph;

    explicit TableGraph(std::span<const Label> vertices)
        : v(vertices.begin(), vertices.end()), graph(vertices.size()) {
        for (Int i = 0; i < Int(v.size()); ++i)
            toi[v[i]] = i;
    }

    // ラベルで指定した辺を追加し、辺番号を返す。償却 O(1)。
    Int add_edge(const Label &u, const Label &w, typename decltype(graph)::cost_type cost = 1) {
        return graph.add_edge(toi.at(u), toi.at(w), cost);
    }

    auto operator[](const Label &x) const {
        if constexpr (std::is_void_v<Weight>)
            return detail::transform_range(graph[toi.at(x)], [this](Int i) { return v[i]; });
        else
            return detail::transform_range(
                graph[toi.at(x)], [this](auto e) { return std::pair{v[e.first], e.second}; });
    }

    // 隣接頂点のラベルと辺番号を列挙する。O(deg(x))。
    auto to_and_id(const Label &x) const {
        return detail::transform_range(graph.to_and_id(toi.at(x)),
                                       [this](auto e) { return std::pair{v[e.first], e.second}; });
    }

    // 隣接頂点のラベル、重み、辺番号を列挙する。O(deg(x))。
    auto to_and_cost_and_id(const Label &x) const {
        return detail::transform_range(graph.to_and_cost_and_id(toi.at(x)), [this](auto e) {
            return std::tuple{v[std::get<0>(e)], std::get<1>(e), std::get<2>(e)};
        });
    }
};
template <class T> using UnWeightedUnDirectedTableGraph = TableGraph<T, void, false>;
template <class T> using UnWeightedDirectedTableGraph = TableGraph<T, void, true>;
template <class T, class S = Int> using WeightedUnDirectedTableGraph = TableGraph<T, S, false>;
template <class T, class S = Int> using WeightedDirectedTableGraph = TableGraph<T, S, true>;

template <class Range> auto initUnWeightedUnDirectedTableGraph(const Range &v) {
    return UnWeightedUnDirectedTableGraph<typename Range::value_type>(v);
}

template <class Range> auto initUnWeightedDirectedTableGraph(const Range &v) {
    return UnWeightedDirectedTableGraph<typename Range::value_type>(v);
}

template <class S = Int, class Range> auto initWeightedUnDirectedTableGraph(const Range &v) {
    return WeightedUnDirectedTableGraph<typename Range::value_type, S>(v);
}

template <class S = Int, class Range> auto initWeightedDirectedTableGraph(const Range &v) {
    return WeightedDirectedTableGraph<typename Range::value_type, S>(v);
}

template <class G, class... Args> auto add_edge(G &g, Args &&...args) {
    return g.add_edge(std::forward<Args>(args)...);
}

template <class G, class V> auto to_and_id(const G &g, const V &v) {
    return g.to_and_id(v);
}

template <class G, class V> auto to_and_cost(const G &g, const V &v) {
    return g.to_and_cost(v);
}

template <class G, class V> auto to_and_cost_and_id(const G &g, const V &v) {
    return g.to_and_cost_and_id(v);
}
}
