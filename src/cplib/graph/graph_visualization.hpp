#pragma once
#include <cplib/graph/graph.hpp>
#include <cmath>
#include <fstream>
#include <sstream>
#include <stdexcept>

namespace cplib::detail {
struct GraphDrawingPoint {
    double x, y;
};
template <GraphTypes G> auto graphDrawingPositions(const G &g) {
    Int n = g.size();
    std::vector<GraphDrawingPoint> positions(n), forces(n);
    if (n == 0)
        return positions;
    auto score = [&](const std::vector<GraphDrawingPoint> &sites) {
        double result = 0;
        auto cross = [](auto a, auto b, auto c) {
            return (b.x - a.x) * (c.y - a.y) - (b.y - a.y) * (c.x - a.x);
        };
        for (Int u = 0; u < n; ++u) {
            for (Int v = u + 1; v < n; ++v) {
                double shortfall = std::max(0.0, 44 - std::hypot(sites[u].x - sites[v].x, sites[u].y - sites[v].y));
                result += shortfall * shortfall * 100;
            }
        }
        for (Int i = 0; i < g.edge_count(); ++i) {
            const auto &e = g.edge_info[i];
            if (e.src == e.dst)
                continue;
            auto a = sites[e.src], b = sites[e.dst];
            double dx = b.x - a.x, dy = b.y - a.y;
            double lengthSquared = std::max(0.001, dx * dx + dy * dy);
            result += lengthSquared * 0.01;
            for (Int v = 0; v < n; ++v) {
                if (v == e.src || v == e.dst)
                    continue;
                auto point = sites[v];
                double t = std::clamp(((point.x - a.x) * dx + (point.y - a.y) * dy) / lengthSquared, 0.0, 1.0);
                double shortfall = std::max(0.0, 26 - std::hypot(point.x - a.x - t * dx, point.y - a.y - t * dy));
                result += shortfall * shortfall * 100;
            }
            for (Int j = i + 1; j < g.edge_count(); ++j) {
                const auto &f = g.edge_info[j];
                if (f.src == f.dst || e.src == f.src || e.src == f.dst || e.dst == f.src || e.dst == f.dst)
                    continue;
                auto c = sites[f.src], d = sites[f.dst];
                if (cross(a, b, c) * cross(a, b, d) < -0.001 && cross(c, d, a) * cross(c, d, b) < -0.001)
                    result += 1000000;
            }
        }
        return result;
    };
    double bestScore = std::numeric_limits<double>::infinity();
    std::vector<GraphDrawingPoint> best;
    for (Int attempt = 0; attempt < 4; ++attempt) {
        std::uint32_t state = 123456789 + std::uint32_t(attempt) * 10007;
        auto random = [&]() {
            state ^= state << 13;
            state ^= state >> 17;
            state ^= state << 5;
            return double(state) / std::numeric_limits<std::uint32_t>::max() * 200 - 100;
        };
        for (auto &p : positions)
            p = {random(), random()};
        double ideal = std::sqrt(40000.0 / n);
        for (Int iteration = 0; iteration < 300; ++iteration) {
            for (Int v = 0; v < n; ++v)
                forces[v] = {-positions[v].x * 2, -positions[v].y * 2};
            for (Int u = 0; u < n; ++u) {
                for (Int v = u + 1; v < n; ++v) {
                    double dx = positions[u].x - positions[v].x, dy = positions[u].y - positions[v].y;
                    double force = ideal * ideal / std::max(0.01, dx * dx + dy * dy);
                    forces[u].x += dx * force;
                    forces[u].y += dy * force;
                    forces[v].x -= dx * force;
                    forces[v].y -= dy * force;
                }
            }
            for (const auto &e : g.edge_info) {
                if (e.src == e.dst)
                    continue;
                double dx = positions[e.src].x - positions[e.dst].x;
                double dy = positions[e.src].y - positions[e.dst].y;
                double force = std::hypot(dx, dy) / ideal;
                forces[e.src].x -= dx * force;
                forces[e.src].y -= dy * force;
                forces[e.dst].x += dx * force;
                forces[e.dst].y += dy * force;
            }
            double temperature = 12.0 * (1 - iteration / 300.0) + 0.02;
            for (Int v = 0; v < n; ++v) {
                double length = std::hypot(forces[v].x, forces[v].y);
                if (length > 0) {
                    double step = std::min(length, temperature) / length;
                    positions[v].x += forces[v].x * step;
                    positions[v].y += forces[v].y * step;
                }
            }
        }
        double minX = positions[0].x, maxX = minX, minY = positions[0].y, maxY = minY;
        for (auto p : positions) {
            minX = std::min(minX, p.x);
            maxX = std::max(maxX, p.x);
            minY = std::min(minY, p.y);
            maxY = std::max(maxY, p.y);
        }
        double scale = 244 / std::max({1.0, maxX - minX, maxY - minY});
        for (auto &p : positions)
            p = {160 + (p.x - (minX + maxX) / 2) * scale, 160 + (p.y - (minY + maxY) / 2) * scale};
        double value = score(positions);
        if (value < bestScore) {
            bestScore = value;
            best = positions;
        }
    }
    return best;
}
template <GraphTypes G> void graphDrawingSvg(std::ostream &out, const G &g, Int indexed, Int serial) {
    auto positions = graphDrawingPositions(g);
    out << "<svg viewBox=\"0 0 320 320\" xmlns=\"http://www.w3.org/2000/svg\">\n"
        << "<defs><marker id=\"arrow" << serial
        << "\" viewBox=\"0 0 10 10\" refX=\"9\" refY=\"5\""
        << " markerWidth=\"6\" markerHeight=\"6\" orient=\"auto\">"
        << "<path d=\"M0 0 L10 5 L0 10 Z\" fill=\"#64748b\"/></marker></defs>\n";
    for (const auto &e : g.edge_info) {
        auto a = positions[e.src], b = positions[e.dst];
        if (e.src == e.dst) {
            out << "<path d=\"M " << a.x - 10 << ' ' << a.y - 12 << " C "
                << a.x - 50 << ' ' << a.y - 65 << ' ' << a.x + 50 << ' ' << a.y - 65 << ' '
                << a.x + 10 << ' ' << a.y - 12 << "\" fill=\"none\" stroke=\"#64748b\"";
        } else {
            double dx = b.x - a.x, dy = b.y - a.y, length = std::max(0.01, std::hypot(dx, dy));
            out << "<path d=\"M " << a.x + dx / length * 18 << ' ' << a.y + dy / length * 18
                << " Q " << (a.x + b.x) / 2 - dy * 0.15 << ' ' << (a.y + b.y) / 2 + dx * 0.15
                << ' ' << b.x - dx / length * 20 << ' ' << b.y - dy / length * 20
                << "\" fill=\"none\" stroke=\"#64748b\"";
        }
        if constexpr (G::is_directed)
            out << " marker-end=\"url(#arrow" << serial << ")\"";
        out << "/>\n";
    }
    for (Int v = 0; v < g.size(); ++v) {
        auto p = positions[v];
        out << "<circle cx=\"" << p.x << "\" cy=\"" << p.y << "\" r=\"17\" fill=\"#eff6ff\" stroke=\"#2563eb\"/>"
            << "<text x=\"" << p.x << "\" y=\"" << p.y << "\" text-anchor=\"middle\" dy=\".35em\">"
            << v + indexed << "</text>\n";
    }
    out << "</svg>\n";
}
}
namespace cplib {
// 条件を満たすグラフをHTMLへ保存して件数を返す。max_results=0は無制限。負値はinvalid_argument。
// Q個の列挙・判定費用に加え、一件の配置 O((V+E)²)、出力 O(V+E)、追加領域 O(V)。
template <class Range, class Predicate>
Int visualize_graphs(Range &&graphs, Predicate predicate, const std::string &output = "graphs.html",
                     Int max_results = 0, Int indexed = 0) {
    if (max_results < 0)
        throw std::invalid_argument("表示件数の上限は非負で指定してください");
    std::ofstream file(output);
    if (!file)
        throw std::runtime_error("HTMLファイルを開けません");
    file << "<!doctype html><html lang=\"ja\"><meta charset=\"utf-8\"><title>Graph enumeration</title>"
         << "<style>body{font-family:system-ui;background:#f8fafc}"
         << "main{display:grid;grid-template-columns:repeat(auto-fill,minmax(280px,1fr));gap:16px}"
         << "article{background:white;padding:16px}svg{width:100%;max-height:340px}</style>"
         << "<h1>条件を満たすグラフ</h1><main>\n";
    Int count = 0;
    bool limited = false;
    for (const auto &g : graphs) {
        if (!predicate(g))
            continue;
        ++count;
        file << "<article data-edge-count=\"" << g.edge_count() << "\"><h2>#" << count
             << " / " << g.size() << " 頂点 / " << g.edge_count() << " 辺</h2>\n";
        detail::graphDrawingSvg(file, g, indexed, count);
        file << "<details><summary>辺リスト</summary><pre>" << g.size() << ' ' << g.edge_count() << '\n';
        for (const auto &e : g.edge_info)
            file << e.src + indexed << ' ' << e.dst + indexed << '\n';
        file << "</pre></details></article>\n";
        if (max_results > 0 && count >= max_results) {
            limited = true;
            break;
        }
    }
    file << "</main><footer>表示件数: " << count;
    if (limited)
        file << " / 表示上限に達したため探索を終了しました（全件数は未集計）";
    file << "</footer></html>\n";
    if (!file)
        throw std::runtime_error("HTMLファイルの書き込みに失敗しました");
    return count;
}
}
