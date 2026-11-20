#pragma once
#include <cplib/geometry/private/persistent_edge_order.hpp>
#include <cplib/geometry/voronoi_diagram.hpp>
#include <optional>

namespace cplib {

// ボロノイ図の各 x 区間の辺の並びを、走査線と永続 AVL 木で保存する。
// 構築 O(N log N)、領域 O(N log N)、最悪問い合わせ O(log N) 回の算術・比較。
// BigInt の各演算の費用は座標の桁数に依存する。無限領域を含め、問い合わせ点の座標範囲に制約はない。
class VoronoiNearestNeighbor {
    using P = VoronoiDiagram::P;
    using Rational = VoronoiDiagram::Rational;

    struct Bisector {
        BigInt a, b, c; // ax+by=c、b>0。
        Int above;
    };
    struct Event {
        Int vertex;
        int edge;
        bool start;
    };
    std::vector<Bisector> bisectors_;
    detail::PersistentEdgeOrder order_;
    int initial_root_ = -1;
    std::vector<Rational> xs_, bottom_breaks_;
    std::vector<int> roots_;
    std::vector<Int> bottom_sites_;

    static int sign(const BigInt &x) { return (x > 0) - (x < 0); }

    int compare_point(int edge, const P &p) const {
        const auto &[a, b, c, above] = bisectors_[edge];
        return sign((c * p.x.den - a * p.x.num) * p.y.den - b * p.y.num * p.x.den);
    }

    int compare_slopes(int edge, int other) const {
        const auto &a = bisectors_[edge], &b = bisectors_[other];
        return sign(b.a * a.b - a.a * b.b);
    }

    int compare_initial(int edge, int other) const {
        if (auto slope = compare_slopes(edge, other); slope != 0)
            return -slope;
        const auto &a = bisectors_[edge], &b = bisectors_[other];
        return sign(a.c * b.b - b.c * a.b);
    }

public:
    explicit VoronoiNearestNeighbor(const VoronoiDiagram &diagram) {
        std::vector<Event> events;
        for (const auto &edge : diagram.edges) {
            auto [a, b] = edge.sites;
            const bool right_to_left = diagram.sites[a].y < diagram.sites[b].y;
            if (diagram.sites[a].y == diagram.sites[b].y)
                continue; // 垂直な境界は x の直後の区間を使えばよい。
            const auto [lower, upper] = right_to_left ? std::pair{a, b} : std::pair{b, a};
            const auto &p = diagram.sites[lower], &q = diagram.sites[upper];
            Rational dx = (q.x - p.x) * Rational(2);
            Rational dy = (q.y - p.y) * Rational(2);
            Rational c = dot(q - p, q + p);
            int index = static_cast<int>(bisectors_.size());
            bisectors_.push_back({dx.num * dy.den * c.den, dy.num * dx.den * c.den,
                                  c.num * dx.den * dy.den, upper});
            auto [left, right] = edge.vertices;
            if (right_to_left)
                std::swap(left, right);
            if (left >= 0)
                events.push_back({left, index, true});
            else
                initial_root_ = order_.insert(initial_root_, index,
                                              [&](int first, int second) { return compare_initial(first, second); });
            if (right >= 0)
                events.push_back({right, index, false});
        }

        auto event_x = [&](const Event &event) -> const Rational & { return diagram.vertices[event.vertex].x; };
        std::sort(events.begin(), events.end(), [&](const Event &a, const Event &b) {
            if (event_x(a) != event_x(b))
                return event_x(a) < event_x(b);
            return a.start < b.start;
        });
        int root = initial_root_;
        for (std::size_t i = 0; i < events.size(); ++i) {
            const auto &event = events[i];
            auto compare = [&](int a, int b) {
                if (a == b)
                    return 0;
                int height = -compare_point(b, diagram.vertices[event.vertex]);
                return height != 0 ? height : (event.start ? 1 : -1) * compare_slopes(a, b);
            };
            root = event.start ? order_.insert(root, event.edge, compare) : order_.erase(root, event.edge, compare);
            if (i + 1 == events.size() || event_x(events[i]) != event_x(events[i + 1])) {
                xs_.push_back(event_x(events[i]));
                roots_.push_back(root);
            }
        }

        // y が十分小さい領域では y が最小の母点が勝ち、同じ y の中では最も近い x が勝つ。
        for (Int i = 0; i < Int(diagram.sites.size()); ++i) {
            if (diagram.representative[i] != i)
                continue;
            if (!bottom_sites_.empty() && diagram.sites[i].y < diagram.sites[bottom_sites_.front()].y)
                bottom_sites_.clear();
            if (bottom_sites_.empty() || diagram.sites[i].y == diagram.sites[bottom_sites_.front()].y)
                bottom_sites_.push_back(i);
        }
        std::sort(bottom_sites_.begin(), bottom_sites_.end(),
                  [&](Int a, Int b) { return diagram.sites[a].x < diagram.sites[b].x; });
        for (std::size_t i = 1; i < bottom_sites_.size(); ++i)
            bottom_breaks_.push_back((diagram.sites[bottom_sites_[i - 1]].x + diagram.sites[bottom_sites_[i]].x)
                                     / Rational(2));
    }

    template <class T, class Id>
    explicit VoronoiNearestNeighbor(std::span<const Point<T, Id>> points)
        : VoronoiNearestNeighbor(voronoi_diagram(points)) {}

    template <class T, class Id>
    explicit VoronoiNearestNeighbor(const std::vector<Point<T, Id>> &points)
        : VoronoiNearestNeighbor(std::span<const Point<T, Id>>(points)) {}

    // 最も近い母点の入力添字を返す。同距離ならいずれか 1 つ。母点が空なら nullopt。
    std::optional<Int> nearest(const P &query) const {
        if (bottom_sites_.empty())
            return std::nullopt;
        auto x = std::upper_bound(xs_.begin(), xs_.end(), query.x);
        int root = x == xs_.begin() ? initial_root_ : roots_[x - xs_.begin() - 1];
        int below = -1;
        while (root >= 0) {
            const auto &node = order_[root];
            if (compare_point(node.edge, query) <= 0) {
                below = node.edge;
                root = node.right;
            } else {
                root = node.left;
            }
        }
        if (below >= 0)
            return bisectors_[below].above;
        auto split = std::upper_bound(bottom_breaks_.begin(), bottom_breaks_.end(), query.x);
        return bottom_sites_[split - bottom_breaks_.begin()];
    }

    template <class T, class Id> std::optional<Int> nearest(const Point<T, Id> &query) const {
        return nearest(detail::voronoiPoint(query));
    }
};

}
