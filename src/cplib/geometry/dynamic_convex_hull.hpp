#pragma once
#include <cplib/geometry/convex_polygon.hpp>
#include <cplib/geometry/private/convex_hull_bridge.hpp>
#include <memory>
#include <optional>
#include <stdexcept>
#include <tuple>

namespace cplib {

// 座標の多重集合を保持する完全動的凸包。座標の事前登録は不要。
// 符号付き整数座標、各軸の活動点の座標幅は 2*10^9 以下。内部演算は __int128。
// 挿入・削除 O(log^2(N+1))、支持点 O(log(N+1))、領域 O(N)。
template <class T = Int> class DynamicConvexHull {
    static_assert(std::is_integral_v<T> && std::is_signed_v<T> && sizeof(T) <= sizeof(Int));
    using P = Point<__int128>;
    using Bridge = std::pair<P, P>;
    struct Node {
        std::unique_ptr<Node> left, right;
        P first, last;
        __int128 xmin, xmax;
        std::array<Bridge, 2> bridges;
        Int height = 1;
        std::size_t copies = 1;

        explicit Node(P p) : first(p), last(p), xmin(p.x), xmax(p.x), bridges{{{p, p}, {-p, -p}}} {}
        bool leaf() const { return !left; }
        const Node *child(bool next, bool reversed) const {
            return next != reversed ? right.get() : left.get();
        }
        P minimum(bool reversed) const { return reversed ? -last : first; }
        P maximum(bool reversed) const { return reversed ? -first : last; }
    };
    using Link = std::unique_ptr<Node>;
    Link root_;
    std::size_t size_ = 0;

    static bool less(P a, P b) { return std::tie(a.y, a.x) < std::tie(b.y, b.x); }
    static P wide(Point<T> p) { return {p.x, p.y}; }
    static Point<T> narrow(P p) { return {T(p.x), T(p.y)}; }

    static void pull(Node &node) {
        node.height = 1 + std::max(node.left->height, node.right->height);
        node.first = node.left->first;
        node.last = node.right->last;
        node.xmin = std::min(node.left->xmin, node.right->xmin);
        node.xmax = std::max(node.left->xmax, node.right->xmax);
        for (bool reversed : {false, true}) {
            auto l = node.child(false, reversed), r = node.child(true, reversed);
            auto [a, b] = detail::convexHullBridge(
                l, r, r->minimum(reversed).y, [](const Node *p) { return p->leaf(); },
                [&](const Node *p, bool next) { return p->child(next, reversed); },
                [&](const Node *p) { return p->bridges[reversed]; });
            node.bridges[reversed] = {a->minimum(reversed), b->minimum(reversed)};
        }
    }

    static void rotate_left(Link &root) {
        auto top = std::move(root->right);
        root->right = std::move(top->left);
        pull(*root);
        top->left = std::move(root);
        pull(*top);
        root = std::move(top);
    }
    static void rotate_right(Link &root) {
        auto top = std::move(root->left);
        root->left = std::move(top->right);
        pull(*root);
        top->right = std::move(root);
        pull(*top);
        root = std::move(top);
    }
    static void rebalance(Link &root) {
        auto &left = root->left, &right = root->right;
        if (left->height - right->height > 1) {
            if (left->right->height > left->left->height)
                rotate_left(left);
            rotate_right(root);
        } else if (right->height - left->height > 1) {
            if (right->left->height > right->right->height)
                rotate_right(right);
            rotate_left(root);
        } else {
            pull(*root);
        }
    }
    static void insert(Link &root, P p) {
        if (!root) {
            root = std::make_unique<Node>(p);
        } else if (root->leaf()) {
            if (exact_equal(root->first, p)) {
                ++root->copies;
                return;
            }
            auto parent = std::make_unique<Node>(p);
            auto leaf = std::make_unique<Node>(p);
            if (less(p, root->first)) {
                parent->left = std::move(leaf);
                parent->right = std::move(root);
            } else {
                parent->left = std::move(root);
                parent->right = std::move(leaf);
            }
            pull(*parent);
            root = std::move(parent);
        } else {
            insert(less(p, root->right->first) ? root->left : root->right, p);
            rebalance(root);
        }
    }
    static bool erase(Link &root, P p) {
        if (!root)
            return false;
        if (root->leaf()) {
            if (!exact_equal(root->first, p))
                return false;
            if (--root->copies == 0)
                root.reset();
            return true;
        }
        auto &child = less(p, root->right->first) ? root->left : root->right;
        if (!erase(child, p))
            return false;
        if (!root->left) {
            auto next = std::move(root->right);
            root = std::move(next);
        } else if (!root->right) {
            auto next = std::move(root->left);
            root = std::move(next);
        } else {
            rebalance(root);
        }
        return true;
    }
    static void append(const Node *node, bool reversed, P first, P last, Polygon<__int128> &out) {
        if (node->leaf()) {
            if (out.empty() || !exact_equal(out.back(), node->first))
                out.push_back(node->first);
            return;
        }
        const auto &[l, r] = node->bridges[reversed];
        if (!less(l, last))
            append(node->child(false, reversed), reversed, first, last, out);
        else if (!less(first, r))
            append(node->child(true, reversed), reversed, first, last, out);
        else {
            append(node->child(false, reversed), reversed, first, l, out);
            append(node->child(true, reversed), reversed, r, last, out);
        }
    }

public:
    DynamicConvexHull() = default;
    DynamicConvexHull(DynamicConvexHull &&other) noexcept
        : root_(std::move(other.root_)), size_(std::exchange(other.size_, 0)) {}
    DynamicConvexHull &operator=(DynamicConvexHull &&other) noexcept {
        if (this != &other) {
            root_ = std::move(other.root_);
            size_ = std::exchange(other.size_, 0);
        }
        return *this;
    }
    explicit DynamicConvexHull(std::span<const Point<T>> points) {
        for (auto point : points)
            insert(point);
    }

    std::size_t size() const { return size_; }
    bool empty() const { return !root_; }

    // 同じ座標を複数回挿入できる。座標幅の上限を超える挿入は invalid_argument。
    void insert(Point<T> point) {
        auto p = wide(point);
        constexpr __int128 limit = 2000000000;
        if (root_ && (std::max(p.x, root_->xmax) - std::min(p.x, root_->xmin) > limit ||
                      std::max(p.y, root_->last.y) - std::min(p.y, root_->first.y) > limit))
            throw std::invalid_argument("DynamicConvexHull coordinate span exceeds 2e9");
        insert(root_, p);
        ++size_;
    }

    // 存在する場合は 1 個削除し true。存在しなければ false。
    bool erase(Point<T> point) {
        if (!erase(root_, wide(point)))
            return false;
        --size_;
        return true;
    }

    // 内積を最大にする点。同値の場合は任意。空なら nullopt。
    std::optional<Point<T>> support_point(Point<T> direction) const {
        const Node *node = root_.get();
        if (!node)
            return std::nullopt;
        if (direction.x == 0)
            return narrow(direction.y >= 0 ? node->last : node->first);
        auto d = wide(direction);
        bool reversed = d.x > 0;
        if (reversed)
            d = -d;
        while (!node->leaf()) {
            const auto &[l, r] = node->bridges[reversed];
            node = node->child(dot(d, r - l) > 0, reversed);
        }
        return narrow(node->first);
    }

    // 正規化済みの凸包を返す。辺上の点を含め境界に K 個の相異なる点があるとき O(K log(N+1))。
    ConvexPolygon<T> hull() const {
        if (!root_)
            return {};
        Polygon<__int128> boundary;
        for (bool reversed : {false, true})
            append(root_.get(), reversed, root_->minimum(reversed), root_->maximum(reversed), boundary);
        auto normalized = initConvexPolygon(boundary);
        ConvexPolygon<T> result;
        result.reserve(normalized.size());
        for (auto p : normalized)
            result.push_back(narrow(p));
        return result;
    }
};

}
