#pragma once
#include <cplib/geometry/polygon.hpp>
#include <cplib/geometry/private/persistent_edge_order.hpp>
#include <algorithm>
#include <cassert>
#include <cstdint>
#include <type_traits>
#include <vector>

namespace cplib {

// 単純多角形のオンライン点位置探索。整数座標専用。
// 構築 O(N log N)、領域 O(N log N)、包含判定 O(log N)。
// 整数幅: 差・外積・面積累積を__int128へ昇格。B=max(8M^2,2NM^2)、W(B)<=128が条件。
// Mには問い合わせ点も含む。M=10^9なら入力はintで可、10^12ならlong longで可。
// Nと内部ノード数はintに収まること。このNの範囲でM<=10^12なら上記Bは128ビットに収まる。
template <class T> class PersistentPolygonPointLocator {
    static_assert(std::is_integral_v<T> && sizeof(T) <= 8);
    using Wide = __int128_t;
    using P = Point<T>;
    struct Edge {
        P a, b; // 入力多角形の周回方向
        P left, right;
        Wide dx, dy;
    };
    struct Vertical {
        T x, lo, hi;
        bool operator<(const Vertical &other) const {
            return x == other.x ? lo < other.lo : x < other.x;
        }
    };
    std::vector<Edge> edges_;
    detail::PersistentEdgeOrder edge_order_;
    std::vector<T> xs_;
    std::vector<int> roots_;
    std::vector<P> vertices_;
    std::vector<Vertical> vertical_;
    bool ccw_ = true;

    static int sign(Wide x) { return (x > 0) - (x < 0); }
    // y_edge(x) - p.y の符号。dx > 0 で割る前の整数演算。
    int compare_point(int edge, P p) const {
        const Edge &e = edges_[edge];
        return sign(e.dy * (Wide(p.x) - e.left.x) - e.dx * (Wide(p.y) - e.left.y));
    }
    // 更新対象の辺には x に端点がある。x の直前/直後での高さ順を比較する。
    int compare_update(int edge, int other, int side) const {
        if (edge == other)
            return 0;
        const Edge &a = edges_[edge], &b = edges_[other];
        P endpoint = side < 0 ? a.right : a.left;
        int at_x = -compare_point(other, endpoint);
        if (at_x)
            return at_x;
        int slope = sign(a.dy * b.dx - b.dy * a.dx);
        assert(slope != 0); // 同じ側に重なる辺は単純多角形にはない
        return side * slope;
    }

public:
    // 整数幅: クラス冒頭のB=max(8M^2,2NM^2)。内部は128ビット、M=10^9なら入力はint可。
    template <class Id>
    explicit PersistentPolygonPointLocator(const Polygon<T, Id> &polygon) {
        const int n = int(polygon.size());
        assert(n >= 3);
        Wide twice_area = 0;
        vertices_.assign(polygon.begin(), polygon.end());
        std::sort(vertices_.begin(), vertices_.end());
        for (int i = 0; i < n; ++i) {
            P a = polygon[i], b = polygon[(i + 1) % n];
            assert(!exact_equal(a, b));
            twice_area += Wide(a.x) * b.y - Wide(a.y) * b.x;
            if (a.x == b.x) {
                vertical_.push_back({a.x, std::min(a.y, b.y), std::max(a.y, b.y)});
            } else {
                P left = a.x < b.x ? a : b, right = a.x < b.x ? b : a;
                edges_.push_back({a, b, left, right, Wide(right.x) - left.x,
                                  Wide(right.y) - left.y});
            }
        }
        assert(twice_area != 0);
        ccw_ = twice_area > 0;
        std::sort(vertical_.begin(), vertical_.end());
        xs_.reserve(vertices_.size());
        for (P p : vertices_)
            if (xs_.empty() || xs_.back() != p.x)
                xs_.push_back(p.x);
        std::vector<std::vector<int>> starts(xs_.size()), ends(xs_.size());
        for (int i = 0; i < int(edges_.size()); ++i) {
            int l = int(std::lower_bound(xs_.begin(), xs_.end(), edges_[i].left.x) - xs_.begin());
            int r = int(std::lower_bound(xs_.begin(), xs_.end(), edges_[i].right.x) - xs_.begin());
            starts[l].push_back(i);
            ends[r].push_back(i);
        }
        roots_.reserve(xs_.size());
        int root = -1;
        for (int i = 0; i < int(xs_.size()); ++i) {
            for (int edge : ends[i])
                root = edge_order_.erase(root, edge, [&](int a, int b) { return compare_update(a, b, -1); });
            for (int edge : starts[i])
                root = edge_order_.insert(root, edge, [&](int a, int b) { return compare_update(a, b, +1); });
            roots_.push_back(root);
        }
    }

    // 0: 外部、1: 境界、2: 内部。
    // 整数幅: 頂点とpの座標上限をMとしてB=8M^2。内部は128ビット、M=10^9なら入力はint可。
    Int contains(P p) const {
        if (std::binary_search(vertices_.begin(), vertices_.end(), p))
            return 1;
        auto it = std::upper_bound(vertical_.begin(), vertical_.end(), Vertical{p.x, p.y, T(0)});
        if (it != vertical_.begin()) {
            --it;
            if (it->x == p.x && it->lo <= p.y && p.y <= it->hi)
                return 1;
        }
        auto xi = std::upper_bound(xs_.begin(), xs_.end(), p.x);
        if (xi == xs_.begin())
            return 0;
        int root = roots_[xi - xs_.begin() - 1];
        int below = -1;
        while (root >= 0) {
            auto node = edge_order_[root];
            int cmp = compare_point(node.edge, p);
            if (cmp == 0)
                return 1;
            if (cmp < 0) {
                below = node.edge;
                root = node.right;
            } else {
                root = node.left;
            }
        }
        if (below < 0)
            return 0;
        bool inside_above = (edges_[below].b.x > edges_[below].a.x) == ccw_;
        return inside_above ? 2 : 0;
    }
};

// 整数幅: 構築と同じB=max(8M^2,2NM^2)。内部は128ビット、M=10^9なら入力はint可。
template <class T, class Id = void> auto initPersistentPolygonPointLocator(const Polygon<T, Id> &polygon) {
    return PersistentPolygonPointLocator<T>(polygon);
}

// 整数幅: メンバーcontainsと同じB=8M^2。内部は128ビット、M=10^9なら入力はint可。
template <class T> Int contains(const PersistentPolygonPointLocator<T> &locator, Point<T> p) {
    return locator.contains(p);
}

} // namespace cplib
