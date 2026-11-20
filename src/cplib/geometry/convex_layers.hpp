#pragma once
#include <cplib/geometry/base.hpp>
#include <stdexcept>

namespace cplib {
namespace detail {
// (y,x)順の相異なる点の左側凸包。各節点には左右の凸包を結ぶ接線だけを保持する。
// 反転した座標・逆順の添字を使えば右側凸包も同じ処理で管理できる。
class ConvexLayersHull {
    using WidePoint = Point<__int128>;
    struct Node {
        Int first, last;
        Int bridgeLeft, bridgeRight;
        Int left = -1, right = -1;
    };
    const std::vector<WidePoint> &points_;
    bool reversed_;
    std::vector<Node> nodes_;
    Int root_;

    WidePoint point(Int i) const {
        return reversed_ ? -points_[points_.size() - 1 - i] : points_[i];
    }
    static __int128 turn(WidePoint a, WidePoint b, WidePoint c) {
        return cross(b - a, c - a);
    }
    bool leaf(Int node) const {
        return nodes_[node].left == -1;
    }

    // 両子の接線をたどって共通接線を探す。毎回どちらかの子へ降りるのでO(log N)。
    // 厳密な不等号で判定し、接線上の共線点も凸包に残す。
    void pull(Int node) {
        Int left = nodes_[node].left, right = nodes_[node].right;
        const __int128 split = point(nodes_[right].first).y;
        while (!leaf(left) || !leaf(right)) {
            const auto a = point(nodes_[left].bridgeLeft);
            const auto b = point(nodes_[left].bridgeRight);
            const auto c = point(nodes_[right].bridgeLeft);
            const auto d = point(nodes_[right].bridgeRight);
            if (!leaf(left) && turn(a, b, c) > 0) {
                left = nodes_[left].left;
            } else if (!leaf(right) && turn(b, c, d) > 0) {
                right = nodes_[right].right;
            } else if (leaf(left)) {
                right = nodes_[right].left;
            } else if (leaf(right)) {
                left = nodes_[left].right;
            } else {
                const __int128 s = turn(a, b, c), t = turn(b, a, d);
                // 接線の交点と分割線の位置を、除算せず比較する。
                // 絶対座標を掛けず差を使うため、必要な整数幅は座標幅だけに依存する。
                if (s + t == 0 || s * (d.y - split) + t * (c.y - split) < 0)
                    left = nodes_[left].right;
                else
                    right = nodes_[right].left;
            }
        }
        nodes_[node].bridgeLeft = nodes_[left].first;
        nodes_[node].bridgeRight = nodes_[right].first;
    }
    Int build(Int first, Int last) {
        Int node = nodes_.size();
        nodes_.push_back({first, last, first, first});
        if (last - first > 1) {
            Int mid = first + (last - first) / 2;
            Int left = build(first, mid), right = build(mid, last);
            nodes_[node].left = left;
            nodes_[node].right = right;
            pull(node);
        }
        return node;
    }
    Int erase(Int node, Int index) {
        if (index < nodes_[node].first || index >= nodes_[node].last)
            return node;
        if (leaf(node))
            return -1;
        Int left = erase(nodes_[node].left, index);
        Int right = erase(nodes_[node].right, index);
        // 空になった子を取り除き、残る子を昇格させる。
        if (left == -1)
            return right;
        if (right == -1)
            return left;
        nodes_[node].left = left;
        nodes_[node].right = right;
        pull(node);
        return node;
    }
    void collect(Int node, Int first, Int last, std::vector<Int> &out) const {
        const auto &current = nodes_[node];
        if (leaf(node)) {
            Int index = current.first;
            out.push_back(reversed_ ? Int(points_.size()) - 1 - index : index);
        } else if (last <= current.bridgeLeft) {
            collect(current.left, first, last, out);
        } else if (first >= current.bridgeRight) {
            collect(current.right, first, last, out);
        } else {
            collect(current.left, first, current.bridgeLeft, out);
            collect(current.right, current.bridgeRight, last, out);
        }
    }

public:
    explicit ConvexLayersHull(const std::vector<WidePoint> &points, bool reversed)
        : points_(points), reversed_(reversed), root_(-1) {
        if (!points.empty()) {
            nodes_.reserve(2 * points.size() - 1);
            root_ = build(0, points.size());
        }
    }
    // 元のpoints内の添字を削除する。存在する点を一度だけ指定すること。O(log^2 N)。
    void erase(Int index) {
        root_ = erase(root_, reversed_ ? Int(points_.size()) - 1 - index : index);
    }
    // 凸包上のK点を元のpoints内の添字で追記する。O(K log(N+1))。
    void append_hull(std::vector<Int> &out) const {
        if (root_ != -1)
            collect(root_, 0, Int(points_.size()) - 1, out);
    }
};
}

// 凸包とその辺上の点を同時に削除し、各点が削除される層番号を入力順で返す。最外層は1。
// 入力条件: 順不同のPoint<T>列。Tは64ビット以下の符号付き整数、各軸の座標幅は2*10^9以下。
// 事前の凸包構築は不要。同一点も許し同じ層にする。空入力は空、全点共線なら全て1。
// 座標幅が制約を超える場合はinvalid_argument。入力は変更しない。
// 時間 O(N log^2(N+1))、領域 O(N)。削除可能な凸包で各点を一度ずつ削除する。
// 整数幅: 内部で__int128に拡張。D=max(1,各軸の座標幅)として外積と接線比較は最大4D^3。
// D<=2*10^9なら128ビットで十分。Point<int>もそのまま渡せる。
template <class T, class Id = void> requires (std::is_integral_v<T> && std::is_signed_v<T> && sizeof(T) <= sizeof(Int))
std::vector<Int> convex_layers(std::span<const Point<T, Id>> input) {
    const Int n = input.size();
    if (n == 0)
        return {};
    std::vector<Int> order(n), group(n);
    std::iota(order.begin(), order.end(), Int(0));
    std::sort(order.begin(), order.end(), [&](Int a, Int b) {
        return input[a].y == input[b].y ? input[a].x < input[b].x : input[a].y < input[b].y;
    });
    std::vector<Point<__int128>> points;
    points.reserve(n);
    __int128 minX = input[0].x, maxX = input[0].x;
    for (Int id : order) {
        Point<__int128> p{input[id].x, input[id].y};
        if (points.empty() || !exact_equal(points.back(), p))
            points.push_back(p);
        group[id] = Int(points.size()) - 1;
        minX = std::min(minX, p.x);
        maxX = std::max(maxX, p.x);
    }
    if (maxX - minX > 2000000000 || points.back().y - points.front().y > 2000000000)
        throw std::invalid_argument("各軸の座標幅は2*10^9以下である必要があります");

    detail::ConvexLayersHull left(points, false), right(points, true);
    std::vector<Int> layers(points.size()), boundary;
    boundary.reserve(2 * points.size());
    Int remaining = points.size();
    for (Int layer = 1; remaining > 0; ++layer) {
        boundary.clear();
        // 必ず両側を取得してから削除する。途中で現れた点は次の層に属する。
        left.append_hull(boundary);
        right.append_hull(boundary);
        for (Int index : boundary) {
            if (layers[index] != 0)
                continue;
            layers[index] = layer;
            --remaining;
            left.erase(index);
            right.erase(index);
        }
    }
    std::vector<Int> result(n);
    for (Int i = 0; i < n; ++i)
        result[i] = layers[group[i]];
    return result;
}

// vector版。入力条件・返り値はspan版と同じ。時間 O(N log^2(N+1))、領域 O(N)。
// 順不同・重複可、64ビット以下の符号付き整数座標、各軸の座標幅2*10^9以下。最外層は1。
template <class T, class Id = void> requires (std::is_integral_v<T> && std::is_signed_v<T> && sizeof(T) <= sizeof(Int))
std::vector<Int> convex_layers(const std::vector<Point<T, Id>> &points) {
    return convex_layers<T, Id>(std::span<const Point<T, Id>>(points));
}
}
