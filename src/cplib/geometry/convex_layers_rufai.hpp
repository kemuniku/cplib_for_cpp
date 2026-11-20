#pragma once
#include <cplib/geometry/base.hpp>
#include <stdexcept>

namespace cplib {
namespace detail {
// Rufai–Richards, A Simple Convex Layers Algorithm (arXiv:1702.06829) の
// 4方向のhull treeを使う実験版。各点は1つの双方向リストにだけ属する。
// 原論文との差: 静的再帰構築、below判定を接線走査と統合、削除済み点を根で除去してから層を確定。
// 走査方法を変更しているため、原論文のO(N log N)の償却保証はこの版では未保証。
class RufaiQuarterHull {
    using P = Point<__int128>;
    struct Node { int head = -1, tail = -1, cursorLeft = -1, cursorRight = -1; };
    struct Range {
        int first = -1, last = -1;
        explicit operator bool() const { return first >= 0; }
    };
    std::vector<P> points_;
    std::vector<int> original_, prev_, next_, owner_, pending_;
    std::vector<Node> nodes_;
    const std::vector<Int> &layers_;

    __int128 turn(int a, int b, int c) const {
        return cross(points_[b] - points_[a], points_[c] - points_[a]);
    }
    void build(int k, int first, int last, const std::vector<int> &ids) {
        if (ids.empty()) return;
        std::vector<int> hull;
        for (int i : ids) {
            if (!hull.empty() && points_[i].y < points_[hull.back()].y) continue;
            while (hull.size() > 1 && turn(hull[hull.size() - 2], hull.back(), i) > 0)
                hull.pop_back();
            hull.push_back(i);
        }
        nodes_[k] = {hull.front(), hull.back(), hull.front(), hull.back()};
        for (std::size_t j = 0; j < hull.size(); ++j) {
            int i = hull[j];
            owner_[i] = k;
            prev_[i] = j ? hull[j - 1] : -1;
            next_[i] = j + 1 < hull.size() ? hull[j + 1] : -1;
        }
        if (last - first == 1) return;
        int mid = first + (last - first) / 2;
        std::vector<int> left, right;
        for (int i : ids)
            if (owner_[i] != k) (i < mid ? left : right).push_back(i);
        build(k * 2, first, mid, left);
        build(k * 2 + 1, mid, last, right);
    }
    void cursor_left(int k, int i) {
        nodes_[k].cursorLeft = i;
        nodes_[k].cursorRight = std::max(nodes_[k].cursorRight, i);
    }
    void cursor_right(int k, int i) {
        nodes_[k].cursorRight = i;
        nodes_[k].cursorLeft = std::min(nodes_[k].cursorLeft, i);
    }
    int tangent_left(int k, int p) {
        if (p < 0) return nodes_[k].head;
        int tail = nodes_[k].tail;
        if (tail < p || points_[tail].y < points_[p].y) return -1;
        int i = nodes_[k].cursorLeft;
        while (i < p || points_[i].y < points_[p].y) i = next_[i];
        while (prev_[i] > p && points_[prev_[i]].y >= points_[p].y && turn(p, i, prev_[i]) >= 0)
            i = prev_[i];
        while (next_[i] >= 0 && turn(p, i, next_[i]) > 0) i = next_[i];
        cursor_left(k, i);
        return i;
    }
    int tangent_right(int k, int q) {
        if (q < 0) return nodes_[k].tail;
        int head = nodes_[k].head;
        if (head > q || points_[head].y > points_[q].y) return -1;
        int i = nodes_[k].cursorRight;
        while (i > q || points_[i].y > points_[q].y) i = prev_[i];
        while (next_[i] >= 0 && next_[i] < q && points_[next_[i]].y <= points_[q].y && turn(i, q, next_[i]) >= 0)
            i = next_[i];
        while (prev_[i] >= 0 && turn(i, q, prev_[i]) > 0) i = prev_[i];
        cursor_right(k, i);
        return i;
    }
    Range exposed(int k, int p, int q) {
        if (std::size_t(k) >= nodes_.size() || nodes_[k].head < 0) return {};
        int a = tangent_left(k, p), b = tangent_right(k, q);
        if (a < 0 || b < 0 || a > b) return {};
        if (p >= 0 && q >= 0 && turn(p, q, a) < 0) return {};
        return {a, b};
    }
    void bridge(Range &a, Range &b, int p, int q) const {
        if (!a || !b) return;
        if (p >= 0 && turn(p, b.first, a.first) < 0) { a = {}; return; }
        if ((q >= 0 && turn(a.last, q, b.last) < 0) || points_[a.last].y > points_[b.last].y) {
            b = {}; return;
        }
        int l = a.first, r = b.last;
        for (;;) {
            if (next_[l] >= 0 && next_[l] <= a.last && turn(l, r, next_[l]) >= 0) l = next_[l];
            else if (prev_[l] >= a.first && turn(l, r, prev_[l]) > 0) l = prev_[l];
            else if (next_[r] >= 0 && next_[r] <= b.last && turn(l, r, next_[r]) > 0) r = next_[r];
            else if (prev_[r] >= b.first && turn(l, r, prev_[r]) >= 0) r = prev_[r];
            else break;
        }
        a.last = l;
        b.first = r;
    }
    // 根の部分点列を切り出し、露出した子の点列を持ち上げて穴を埋める。
    void remove(int k, int a, int b) {
        int p = prev_[a], q = next_[b];
        for (int i = a;; i = next_[i]) {
            owner_[i] = -1;
            if (i == b) break;
        }
        if (p >= 0) next_[p] = q; else nodes_[k].head = q;
        if (q >= 0) prev_[q] = p; else nodes_[k].tail = p;
        prev_[a] = next_[b] = -1;
        Range left = exposed(k * 2, p, q), right = exposed(k * 2 + 1, p, q);
        bridge(left, right, p, q);
        if (left) remove(k * 2, left.first, left.last);
        if (right) remove(k * 2 + 1, right.first, right.last);
        if (left && right) {
            next_[left.last] = right.first;
            prev_[right.first] = left.last;
        }
        if (left || right) {
            int first = left ? left.first : right.first;
            int last = right ? right.last : left.last;
            if (p >= 0) next_[p] = first; else nodes_[k].head = first;
            if (q >= 0) prev_[q] = last; else nodes_[k].tail = last;
            prev_[first] = p; next_[last] = q;
            for (int i = first;; i = next_[i]) {
                owner_[i] = k;
                if (k == 1 && layers_[original_[i]] != 0) pending_.push_back(i);
                if (i == last) break;
            }
        }
        nodes_[k].cursorLeft = p >= 0 ? p : nodes_[k].head;
        nodes_[k].cursorRight = q >= 0 ? q : nodes_[k].tail;
    }

public:
    RufaiQuarterHull(const std::vector<P> &points, int direction, const std::vector<Int> &layers)
        : points_(points), original_(points.size()), prev_(points.size(), -1),
          next_(points.size(), -1), owner_(points.size(), -1), nodes_(4 * points.size() + 4), layers_(layers) {
        for (auto &p : points_)
            for (int d = 0; d < direction; ++d) p = {-p.y, p.x};
        std::iota(original_.begin(), original_.end(), 0);
        std::sort(original_.begin(), original_.end(), [&](int a, int b) { return points_[a] < points_[b]; });
        std::vector<P> sorted;
        sorted.reserve(points.size());
        for (int i : original_) sorted.push_back(points_[i]);
        points_.swap(sorted);
        std::vector<int> ids(points.size());
        std::iota(ids.begin(), ids.end(), 0);
        build(1, 0, points.size(), ids);
    }
    void append_hull(std::vector<int> &out) const {
        for (int i = nodes_[1].head; i >= 0; i = next_[i]) out.push_back(original_[i]);
    }
    void peel() {
        if (nodes_[1].head >= 0) remove(1, nodes_[1].head, nodes_[1].tail);
        // 他方向ですでに削除した点が新たな根に現れたら、次の層を確定する前に取り除く。
        // 根全体の再走査はせず、新たに根へ移った点だけをキューへ入れる。
        while (!pending_.empty()) {
            int i = pending_.back(); pending_.pop_back();
            if (owner_[i] == 1) remove(1, i, i);
        }
    }
};
}

// Rufai–Richardsの4方向hull treeを基にした比較用実験版。原論文そのままの移植ではない。
// 入力条件: 順不同・重複可。64ビット以下の符号付き整数、各軸の座標幅2*10^9以下。
// 入力順の層番号（最外層1）を返す。空入力は空。各層では辺上の点も同時に除去する。
// 領域 O(N)。原論文は時間 O(N log N)だが、変更した走査・結合処理の償却計算量は未保証。
// 安全側の時間上限は O(N^2 log(N+1))。実測値はaudit/convex_layersで比較する。
// 座標幅超過はinvalid_argument。差・外積は__int128で計算する（D=max(1,各軸幅)として2D^2以下）。
template <class T, class Id = void> requires (std::is_integral_v<T> && std::is_signed_v<T> && sizeof(T) <= sizeof(Int))
std::vector<Int> convex_layers_rufai(std::span<const Point<T, Id>> input) {
    if (input.empty()) return {};
    if (input.size() > std::size_t(std::numeric_limits<int>::max() / 4 - 1))
        throw std::length_error("点数が多すぎます");
    const int n = input.size();
    std::vector<int> order(n), group(n);
    std::iota(order.begin(), order.end(), 0);
    std::sort(order.begin(), order.end(), [&](int a, int b) { return input[a] < input[b]; });
    std::vector<Point<__int128>> points;
    __int128 minY = input[0].y, maxY = input[0].y;
    for (int i : order) {
        Point<__int128> p{input[i].x, input[i].y};
        if (points.empty() || !exact_equal(points.back(), p)) points.push_back(p);
        group[i] = int(points.size()) - 1;
        minY = std::min(minY, p.y); maxY = std::max(maxY, p.y);
    }
    if (points.back().x - points.front().x > 2000000000 || maxY - minY > 2000000000)
        throw std::invalid_argument("各軸の座標幅は2*10^9以下である必要があります");
    std::vector<Int> layers(points.size());
    std::vector<detail::RufaiQuarterHull> trees;
    trees.reserve(4);
    for (int d = 0; d < 4; ++d) trees.emplace_back(points, d, layers);
    std::vector<int> boundary;
    int remaining = points.size();
    for (Int layer = 1; remaining; ++layer) {
        boundary.clear();
        for (const auto &tree : trees) tree.append_hull(boundary);
        for (int i : boundary)
            if (layers[i] == 0) { layers[i] = layer; --remaining; }
        for (auto &tree : trees) tree.peel();
    }
    std::vector<Int> result(n);
    for (int i = 0; i < n; ++i) result[i] = layers[group[i]];
    return result;
}
// vector版。入力条件・返り値・実験版としての制約はspan版と同じ。
template <class T, class Id = void> requires (std::is_integral_v<T> && std::is_signed_v<T> && sizeof(T) <= sizeof(Int))
std::vector<Int> convex_layers_rufai(const std::vector<Point<T, Id>> &points) {
    return convex_layers_rufai<T, Id>(std::span<const Point<T, Id>>(points));
}
}
