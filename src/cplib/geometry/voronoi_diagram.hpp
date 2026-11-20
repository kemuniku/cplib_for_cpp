#pragma once
#include <cplib/geometry/delaunay_triangulation.hpp>
#include <cplib/math/fractions.hpp>

namespace cplib {

struct VoronoiDiagram {
    using Rational = Fraction<BigInt>;
    using P = Point<Rational>;

    struct Edge {
        std::array<Int, 2> sites;    // 元の入力添字。辺の進行方向から見て左、右。
        std::array<Int, 2> vertices; // 始点、終点。-1 はその向きの無限遠。
    };

    std::vector<P> sites;                  // 元の入力順。同一点も含む。
    std::vector<Int> representative;       // 同座標では最初の入力添字。
    std::vector<P> vertices;               // 有限頂点の厳密な座標。
    std::vector<Edge> edges;
    std::vector<std::vector<Int>> cell_edges;

    // 母点に接する辺の添字。順序は任意。同一母点には同じセルを返す。
    const std::vector<Int> &cell(Int site) const { return cell_edges[representative[site]]; }

    // 辺上の一点。有限端点があればその座標、直線の場合は母点の中点。
    P edge_anchor(Int edge) const {
        const auto &[start, finish] = edges[edge].vertices;
        if (start >= 0)
            return vertices[start];
        if (finish >= 0)
            return vertices[finish];
        const auto [a, b] = edges[edge].sites;
        return (sites[a] + sites[b]) / Rational(2);
    }

    // 辺の方向。無限な辺の向きも vertices の始点側から終点側へ向かう。
    P edge_direction(Int edge) const {
        const auto [a, b] = edges[edge].sites;
        auto delta = sites[b] - sites[a];
        return {-delta.y, delta.x};
    }
};

namespace detail {
template <class T> VoronoiDiagram::Rational voronoiCoordinate(T x) {
    using R = VoronoiDiagram::Rational;
    if constexpr (FloatingPoint<T>) {
        auto value = delaunayBinaryCoordinate(x);
        auto numerator = delaunaySignificand(value);
        if (value.exponent >= 0)
            return R(numerator << value.exponent);
        return R(std::move(numerator), BigInt(1) << -value.exponent, false);
    } else {
        static_assert(std::is_integral_v<T> && sizeof(T) <= sizeof(Int));
        return R(BigInt(x));
    }
}

template <class T, class Id> VoronoiDiagram::P voronoiPoint(const Point<T, Id> &p) {
    return {voronoiCoordinate(p.x), voronoiCoordinate(p.y)};
}

inline VoronoiDiagram::P voronoiCircumcenter(const VoronoiDiagram::P &a,
                                            const VoronoiDiagram::P &b, const VoronoiDiagram::P &c) {
    using R = VoronoiDiagram::Rational;
    auto integral = [](const VoronoiDiagram::P &p) { return p.x.den == 1 && p.y.den == 1; };
    if (integral(a) && integral(b) && integral(c)) {
        Point<BigInt> u{b.x.num - a.x.num, b.y.num - a.y.num}, v{c.x.num - a.x.num, c.y.num - a.y.num};
        BigInt d = BigInt(2) * cross(u, v), uu = norm(u), vv = norm(v);
        return {R(a.x.num * d + uu * v.y - vv * u.y, d),
                R(a.y.num * d + u.x * vv - v.x * uu, d)};
    }
    auto u = b - a, v = c - a;
    R d = R(2) * cross(u, v), uu = norm(u), vv = norm(v);
    return a + VoronoiDiagram::P{(uu * v.y - vv * u.y) / d, (u.x * vv - v.x * uu) / d};
}
}

// ユークリッド・ボロノイ図。O(N log N) 回の算術・比較と O(N) 領域。
// 有限頂点は有理数で厳密に表す。共円点に対する長さ 0 の双対辺は縮約する。
template <class T, class Id = void>
VoronoiDiagram voronoi_diagram(std::span<const Point<T, Id>> points) {
    auto mesh = delaunay_triangulation(points);
    VoronoiDiagram diagram;
    diagram.sites.reserve(points.size());
    for (const auto &p : points)
        diagram.sites.push_back(detail::voronoiPoint(p));
    diagram.representative = std::move(mesh.representative);
    diagram.cell_edges.resize(points.size());

    std::vector<std::pair<VoronoiDiagram::P, Int>> centers;
    centers.reserve(mesh.triangles.size());
    for (auto [a, b, c] : mesh.triangles)
        centers.emplace_back(detail::voronoiCircumcenter(diagram.sites[a], diagram.sites[b], diagram.sites[c]),
                             centers.size());
    std::sort(centers.begin(), centers.end(), [](const auto &a, const auto &b) { return a.first < b.first; });
    std::vector<Int> face_vertex(centers.size());
    for (auto &[center, face] : centers) {
        if (diagram.vertices.empty() || !exact_equal(diagram.vertices.back(), center))
            diagram.vertices.push_back(std::move(center));
        face_vertex[face] = static_cast<Int>(diagram.vertices.size()) - 1;
    }
    auto vertex_of = [&](Int face) { return face < 0 ? Int(-1) : face_vertex[face]; };
    for (std::size_t i = 0; i < mesh.edges.size(); ++i) {
        auto [left, right] = mesh.edge_faces[i];
        Int start = vertex_of(right), finish = vertex_of(left);
        if (start >= 0 && start == finish)
            continue;
        auto [a, b] = mesh.edges[i];
        Int edge = diagram.edges.size();
        diagram.edges.push_back({{a, b}, {start, finish}});
        diagram.cell_edges[a].push_back(edge);
        diagram.cell_edges[b].push_back(edge);
    }
    return diagram;
}

template <class T, class Id = void>
VoronoiDiagram voronoi_diagram(const std::vector<Point<T, Id>> &points) {
    return voronoi_diagram(std::span<const Point<T, Id>>(points));
}

}
