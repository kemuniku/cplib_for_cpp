#include <cplib/graph/directed_mst.hpp>
#include <cplib/graph/dijkstra.hpp>
#include <cplib/graph/min_cost_b_flow.hpp>
#include <cplib/graph/remove_edge.hpp>
#include <cplib/graph/graph_visualization.hpp>
#include <cplib/tree/dsu_on_tree.hpp>
#include <cplib/geometry/closest_pair.hpp>
#include <cplib/geometry/convex_polygon.hpp>
#include <cplib/geometry/manhattan_mst.hpp>
#include <random>
#include <set>
#include <filesystem>

using namespace cplib;

void testDirectedMST(std::mt19937_64 &rng) {
    for (Int trial = 0; trial < 600; ++trial) {
        Int n = 1 + rng() % 5, root = rng() % n;
        WeightedDirectedGraph<Int> g(n);
        std::vector<std::vector<Int>> incoming(n);
        for (Int i = 0, m = rng() % 11; i < m; ++i) {
            Int u = rng() % n, v = rng() % n;
            Int id = g.add_edge(u, v, Int(rng() % 17) - 8);
            if (u != v && v != root)
                incoming[v].push_back(id);
        }
        std::optional<Int> expected;
        std::vector<Int> chosen(n, -1);
        auto search = [&](auto &&self, Int v, Int cost) -> void {
            if (v == n) {
                for (Int x = 0; x < n; ++x) {
                    Int current = x;
                    for (Int steps = 0; steps < n && current != root; ++steps)
                        current = g.edge_info[chosen[current]].src;
                    if (current != root)
                        return;
                }
                if (!expected || cost < *expected)
                    expected = cost;
                return;
            }
            if (v == root) {
                self(self, v + 1, cost);
            } else {
                for (Int id : incoming[v]) {
                    chosen[v] = id;
                    self(self, v + 1, cost + g.edge_info[id].cost);
                }
            }
        };
        search(search, 0, 0);
        auto actual = directedMST(g, root);
        assert(bool(actual) == bool(expected));
        if (actual) {
            assert(actual->cost == *expected && actual->inEdge[root] == -1);
            for (Int x = 0; x < n; ++x) {
                Int v = x;
                for (Int steps = 0; steps < n && v != root; ++steps) {
                    Int id = actual->inEdge[v];
                    assert(g.edge_info[id].dst == v);
                    v = g.edge_info[id].src;
                }
                assert(v == root);
            }
        }
    }
}

void testBFlow(std::mt19937_64 &rng) {
    for (Int trial = 0; trial < 600; ++trial) {
        Int n = 1 + rng() % 4, m = rng() % 7;
        MinCostBFlow graph(n);
        std::vector<MinCostBFlowEdge> edges;
        std::vector<Int> supply(n);
        for (Int v = 0; v + 1 < n; ++v) {
            supply[v] = Int(rng() % 7) - 3;
            supply.back() -= supply[v];
        }
        for (Int v = 0; v < n; ++v)
            graph.add_supply(v, supply[v]);
        for (Int i = 0; i < m; ++i) {
            Int u = rng() % n, v = rng() % n, lower = Int(rng() % 3) - 1;
            Int upper = lower + rng() % 3, cost = Int(rng() % 9) - 4;
            graph.add_edge(u, v, lower, upper, cost);
            edges.push_back({u, v, lower, upper, cost});
        }
        std::optional<Int> expected;
        std::vector<Int> balance(n);
        auto search = [&](auto &&self, Int i, Int cost) -> void {
            if (i == m) {
                if (balance == supply && (!expected || cost < *expected))
                    expected = cost;
                return;
            }
            auto e = edges[i];
            for (Int f = e.lower; f <= e.upper; ++f) {
                balance[e.src] += f;
                balance[e.dst] -= f;
                self(self, i + 1, cost + f * e.cost);
                balance[e.src] -= f;
                balance[e.dst] += f;
            }
        };
        search(search, 0, 0);
        auto actual = graph.solve();
        assert(actual.feasible == bool(expected));
        if (!actual.feasible)
            continue;
        assert(actual.cost == *expected);
        auto potential = graph.get_potential();
        std::fill(balance.begin(), balance.end(), 0);
        for (auto e : graph.get_edges()) {
            assert(e.lower <= e.flow && e.flow <= e.upper);
            balance[e.src] += e.flow;
            balance[e.dst] -= e.flow;
            auto reduced = __int128(e.cost) + potential[e.src] - potential[e.dst];
            if (e.flow < e.upper)
                assert(reduced >= 0);
            if (e.flow > e.lower)
                assert(reduced <= 0);
        }
        assert(balance == supply);
        assert(graph.solve().cost == actual.cost);
    }
}

void testGraphUtilities() {
    auto g = initUnWeightedUnDirectedStaticGraph(3);
    g.add_edge(0, 0);
    g.add_edge(0, 1);
    g.add_edge(1, 2);
    g.build();
    assert(g.degree(0) == 3 && degree(g, 1) == 2);
    Int id = 0;
    for (auto removed : remove_each_edge(g)) {
        assert(removed.edge_count() == 2);
        assert(removed.size() == 3);
        Int j = 0;
        for (Int i = 0; i < 3; ++i)
            if (i != id)
                assert(removed.get_edge(j++) == g.get_edge(i));
        ++id;
    }
    WeightedDirectedGraph<Int> weighted(4);
    weighted.add_edge(0, 1, 7);
    weighted.add_edge(1, 2, 5);
    weighted.add_edge(0, 2, 9);
    auto add = [](Int d, Int edge) { return std::max(d, edge); };
    auto distances = dijkstra(weighted, std::vector<Int>{0, 0}, add, Int(0), Int(100));
    assert((distances == std::vector<Int>{0, 7, 7, 100}));
    struct Cost {
        Int value;
    };
    auto custom = dijkstra(weighted, Int(0), [](Cost d, Int edge) { return Cost{d.value + edge}; },
        Cost{0}, Cost{100}, [](Cost a, Cost b) { return a.value < b.value ? -1 : a.value > b.value ? 1 : 0; });
    assert(custom[2].value == 9);
    std::vector<std::vector<Int>> tree{{1, 2}, {0, 3, 4}, {0}, {1}, {1}};
    Int active = 0;
    std::vector<Int> answers(5);
    auto addVertex = [&](Int v) { active += v + 1; };
    auto clearVertex = [&](Int v) { active -= v + 1; };
    dsuOnTree(tree, addVertex, [&](Int v) { answers[v] = active; }, clearVertex);
    assert((answers == std::vector<Int>{15, 11, 3, 4, 5}) && active == 0);
    std::vector<Int> automatic(5);
    dsuOnTree(tree, [&](RollbackLog &log, Int v) { log.write(active) += v + 1; },
        [&](Int v) { automatic[v] = active; });
    assert(automatic == answers && active == 0);
    bool caught = false;
    try {
        dsuOnTree(tree, addVertex, [&](Int) { throw 1; }, clearVertex);
    } catch (int) {
        caught = true;
    }
    assert(caught && active == 0);
    std::vector<UnWeightedDirectedGraph> graphs{UnWeightedDirectedGraph(0), UnWeightedDirectedGraph(2)};
    graphs[1].add_edge(0, 0);
    graphs[1].add_edge(0, 1);
    auto output = std::filesystem::temp_directory_path() / "cppcplib_nim_port_graphs.html";
    assert(visualize_graphs(graphs, [](const auto &) { return true; }, output.string(), 1) == 1);
    std::ifstream input(output);
    std::string html{std::istreambuf_iterator<char>(input), {}};
    assert(html.find("全件数は未集計") != std::string::npos);
    assert(visualize_graphs(graphs, [](const auto &) { return true; }, output.string()) == 2);
    std::ifstream complete(output);
    std::string allHtml{std::istreambuf_iterator<char>(complete), {}};
    assert(allHtml.find("nan") == std::string::npos);
    assert(allHtml.find("marker-end") != std::string::npos);
    std::filesystem::remove(output);
}

void testGeometry(std::mt19937_64 &rng) {
    for (Int trial = 0; trial < 400; ++trial) {
        std::vector<Point<Int>> points(rng() % 45);
        for (auto &p : points)
            p = {Int(rng() % 31) - 15, Int(rng() % 31) - 15};
        if (points.size() >= 2) {
            Int expected = std::numeric_limits<Int>::max();
            for (Int i = 0; i < Int(points.size()); ++i)
                for (Int j = i + 1; j < Int(points.size()); ++j)
                    expected = std::min(expected, norm(points[i] - points[j]));
            auto [i, j] = closest_pair(points);
            assert(0 <= i && i < j && j < Int(points.size()));
            assert(norm(points[i] - points[j]) == expected);
            assert(closest_pair_distance(points) == std::sqrt(double(expected)));
        }
        auto edges = manhattan_mst(points);
        assert(edges.size() == (points.empty() ? 0 : points.size() - 1));
        Int total = 0;
        for (auto [u, v] : edges)
            total += std::abs(points[u].x - points[v].x) + std::abs(points[u].y - points[v].y);
        Int expected = 0;
        std::vector<Int> distance(points.size(), std::numeric_limits<Int>::max());
        std::vector<bool> used(points.size());
        if (!points.empty())
            distance[0] = 0;
        for (Int step = 0; step < Int(points.size()); ++step) {
            Int v = -1;
            for (Int i = 0; i < Int(points.size()); ++i)
                if (!used[i] && (v < 0 || distance[i] < distance[v]))
                    v = i;
            expected += distance[v];
            used[v] = true;
            for (Int i = 0; i < Int(points.size()); ++i)
                distance[i] = std::min(distance[i], std::abs(points[v].x - points[i].x) + std::abs(points[v].y - points[i].y));
        }
        assert(total == expected);
        auto hull = convex_hull(points, false);
        auto convex = initConvexPolygon(hull);
        if (!points.empty()) {
            Int maxDistance = 0;
            for (auto a : points)
                for (auto b : points)
                    maxDistance = std::max(maxDistance, norm(a - b));
            assert(convex.diameter() == maxDistance);
        }
        for (Int q = 0; q < 30; ++q) {
            Point<Int> p{Int(rng() % 41) - 20, Int(rng() % 41) - 20};
            assert(convex.contains(p) == contains(hull, p, false));
            assert(convex.contains(p, true) == contains(hull, p, true));
            assert(convex.on_edge(p) == on_edge(hull, p));
        }
    }
    std::vector<Point<double>> points{{1e10, 1e10}, {1e10 + 0.01, 1e10}, {1e10 + 1, 1e10 + 1}};
    assert((closest_pair(points) == std::pair<Int, Int>{0, 1}));
    auto square = initPolygon<double>({{0, 0}, {2, 0}, {2, 2}, {0, 2}});
    auto clipped = convex_cut(square, initLine(Point<double>{1, -1}, Point<double>{1, 3}));
    assert(std::abs(area(clipped) - 2) < 1e-9);
    std::unordered_set<PointKey<double>> keys;
    keys.insert(toPointKey(Point<double>{0.0, -0.0}));
    assert(keys.contains(toPointKey(Point<double>{-0.0, 0.0})));
}
int main() {
    std::mt19937_64 rng(1234567);
    testDirectedMST(rng);
    testBFlow(rng);
    testGraphUtilities();
    testGeometry(rng);
}
