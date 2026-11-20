#pragma once
#include <cplib/geometry/half_plane_count.hpp>
#include <tuple>

namespace cplib {

// 三角形の頂点候補 N 点と、数える M 点を前処理する。頂点は入力添字で指定する。
// 時間 O((N+M) log(N+M+1) + N(N+M) log(M+2))、領域 O(N^2+M)。
// 64ビット以下の符号付き整数座標。重複点は別々に数える。
template <detail::half_plane_count::Coordinate T = Int> class TrianglePointCounter {
    using P = Point<T>;
    using Wide = __int128;
    using Count = Int;

    struct Counts {
        Count left = 0, on = 0;
    };
    struct Ray {
        Wide dx, dy;
        bool operator<(const Ray &other) const {
            return detail::half_plane_count::compare_products(dx, other.dy, other.dx, dy) < 0;
        }
    };
    std::vector<P> vertices_;
    std::vector<Counts> rows_, edges_;

    static bool yx_less(P a, P b) { return std::tie(a.y, a.x) < std::tie(b.y, b.x); }
    static Ray ray(P a, P b) { return {Wide(b.x) - a.x, Wide(b.y) - a.y}; }
    const Counts &edge(std::size_t a, std::size_t b) const { return edges_[a * vertices_.size() + b]; }

    Count open_segment(std::size_t a, std::size_t b) const {
        if (vertices_[a].y != vertices_[b].y)
            return edge(a, b).on;
        return rows_[b].left - rows_[a].left - rows_[a].on;
    }

public:
    TrianglePointCounter(std::span<const P> vertices, std::span<const P> points)
        : vertices_(vertices.begin(), vertices.end()), rows_(vertices.size()),
          edges_(vertices.size() * vertices.size()) {
        std::vector<P> sorted(points.begin(), points.end());
        std::sort(sorted.begin(), sorted.end(), yx_less);
        std::vector<std::size_t> order(vertices.size());
        std::iota(order.begin(), order.end(), std::size_t(0));
        std::sort(order.begin(), order.end(),
                  [&](auto a, auto b) { return yx_less(vertices_[a], vertices_[b]); });

        for (auto a : order) {
            const auto origin = vertices_[a];
            auto row_begin = std::lower_bound(sorted.begin(), sorted.end(), P{std::numeric_limits<T>::min(), origin.y},
                                              yx_less);
            auto before = std::lower_bound(row_begin, sorted.end(), origin, yx_less);
            auto after = std::upper_bound(before, sorted.end(), origin, yx_less);
            rows_[a] = {before - row_begin, after - before};

            auto first = std::upper_bound(after, sorted.end(), P{std::numeric_limits<T>::max(), origin.y}, yx_less);
            std::vector<Ray> rays;
            rays.reserve(sorted.end() - first);
            for (auto it = first; it != sorted.end(); ++it)
                rays.push_back(ray(origin, *it));
            std::sort(rays.begin(), rays.end());
            rays.erase(std::unique(rays.begin(), rays.end(), [](Ray l, Ray r) { return !(l < r) && !(r < l); }), rays.end());

            // 上端の y が増える順に点を加え、傾きの区間和を取る。
            std::vector<Count> fenwick(rays.size() + 1);
            auto prefix = [&](std::size_t n) {
                Count result = 0;
                for (; n != 0; n -= n & -n)
                    result += fenwick[n];
                return result;
            };
            auto next = first;
            for (auto b : order) {
                if (vertices_[b].y <= origin.y)
                    continue;
                while (next != sorted.end() && next->y < vertices_[b].y) {
                    auto rank = std::lower_bound(rays.begin(), rays.end(), ray(origin, *next)) - rays.begin();
                    for (std::size_t i = rank + 1; i < fenwick.size(); i += i & -i)
                        ++fenwick[i];
                    ++next;
                }
                auto [l, r] = std::equal_range(rays.begin(), rays.end(), ray(origin, vertices_[b]));
                auto left = prefix(l - rays.begin());
                edges_[a * vertices_.size() + b] = {left, prefix(r - rays.begin()) - left};
            }
        }
    }

    // O(1)。strict=true は境界を除く。false なら退化した三角形の点・線分も含める。
    Count count(std::size_t a, std::size_t b, std::size_t c, bool strict = true) const {
        assert(a < vertices_.size() && b < vertices_.size() && c < vertices_.size());
        auto less = [&](auto i, auto j) { return yx_less(vertices_[i], vertices_[j]); };
        if (less(b, a))
            std::swap(a, b);
        if (less(c, b))
            std::swap(b, c);
        if (less(b, a))
            std::swap(a, b);
        auto ac = ray(vertices_[a], vertices_[c]), ab = ray(vertices_[a], vertices_[b]);
        auto side = detail::half_plane_count::compare_products(ac.dx, ab.dy, ab.dx, ac.dy);
        if (side == 0) {
            if (strict)
                return 0;
            if (exact_equal(vertices_[a], vertices_[c]))
                return rows_[a].on;
            return rows_[a].on + rows_[c].on + open_segment(a, c);
        }

        auto long_side = edge(a, c);
        Counts short_side{edge(a, b).left + edge(b, c).left, edge(a, b).on + edge(b, c).on};
        const Count short_edges_on = short_side.on;
        if (vertices_[a].y != vertices_[b].y && vertices_[b].y != vertices_[c].y) {
            short_side.left += rows_[b].left;
            short_side.on += rows_[b].on;
        }
        const auto &left = side > 0 ? short_side : long_side;
        const auto &right = side > 0 ? long_side : short_side;
        Count result = right.left - left.left - left.on;
        if (!strict)
            result += edge(a, c).on + short_edges_on + rows_[a].on + rows_[b].on + rows_[c].on
                    + (vertices_[a].y == vertices_[b].y ? open_segment(a, b) : 0)
                    + (vertices_[b].y == vertices_[c].y ? open_segment(b, c) : 0);
        return result;
    }
};

template <class T>
TrianglePointCounter(const std::vector<Point<T>> &, const std::vector<Point<T>> &) -> TrianglePointCounter<T>;

}
