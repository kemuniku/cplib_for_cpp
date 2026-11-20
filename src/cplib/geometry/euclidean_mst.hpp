#pragma once
#include <cplib/geometry/delaunay_triangulation.hpp>
#include <cplib/collections/unionfind.hpp>

namespace cplib {
// 厳密なDelaunay三角形分割とKruskal法を使い、MSTを入力の頂点番号の組で返す。
// O(N log N)回の算術・比較、O(N)個の数値を格納。重複点は別頂点。
template <class T, class Id = void>
std::vector<std::pair<Int, Int>> euclidean_mst(std::span<const Point<T, Id>> points) {
    return detail::withDelaunayCoordinates(points, []<class Coord>(std::vector<detail::DelaunaySite<Coord>> sites) {
        auto mesh = detail::buildDelaunayMesh<false>(std::move(sites));
        auto &out = mesh.duplicates;
        std::vector<std::tuple<Coord, Int, Int>> candidates;
        candidates.reserve(mesh.topology.edges.size());
        for (auto [u, v] : mesh.topology.edges) {
            Coord dx = mesh.sites[u].x - mesh.sites[v].x, dy = mesh.sites[u].y - mesh.sites[v].y;
            candidates.emplace_back(dx * dx + dy * dy, u, v);
        }
        std::sort(candidates.begin(), candidates.end());
        UnionFind uf(mesh.sites.size());
        for (const auto &[weight, u, v] : candidates) {
            if (uf.issame(u, v))
                continue;
            uf.unite(u, v);
            auto [a, b] = std::minmax(mesh.sites[u].id, mesh.sites[v].id);
            out.emplace_back(a, b);
            if (uf.count == 1)
                break;
        }
        return std::move(out);
    });
}

// 座標の制約はspan版と同じ。
template <class T, class Id = void> auto euclidean_mst(const std::vector<Point<T, Id>> &p) {
    return euclidean_mst<T, Id>(std::span<const Point<T, Id>>(p));
}

// 座標の組で表すN点のMSTを返す。座標の制約はPoint版と同じ。
template <class T> auto euclidean_mst(std::span<const std::pair<T, T>> p) {
    std::vector<Point<T>> v;
    v.reserve(p.size());
    for (auto [x, y] : p)
        v.push_back({x, y});
    return euclidean_mst(v);
}

// 座標の制約はPoint版と同じ。
template <class T> auto euclidean_mst(const std::vector<std::pair<T, T>> &p) {
    return euclidean_mst<T>(std::span<const std::pair<T, T>>(p));
}
}
