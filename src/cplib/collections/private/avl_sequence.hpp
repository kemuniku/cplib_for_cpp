#pragma once
#include <cplib/utils/backwards_index.hpp>
#include <functional>
#include <memory>
#include <iterator>
#include <sstream>

namespace cplib::detail {
struct NoAvlAction {};

template <class S, class F, bool Aggregate, bool Lazy> class AvlSequence {
    struct Config {
        std::function<S(S, S)> op;
        S e;
        std::function<S(F, S)> mapping;
        std::function<F(F, F)> composition;
        F id;
    };
    struct Node {
        std::unique_ptr<Node> left, right;
        Int height = 1, size = 1;
        S value, product;
        F lazy;
        bool pending = false;
        Node(S value, F id) : value(value), product(value), lazy(id) {
        }
    };
    using Ptr = std::unique_ptr<Node>;
    Ptr root_;
    std::shared_ptr<const Config> config_;

    static Int nodeSize(const Ptr &node) {
        return node ? node->size : 0;
    }
    static Int height(const Ptr &node) {
        return node ? node->height : 0;
    }
    S product(const Ptr &node) const {
        return node ? node->product : config_->e;
    }
    void pull(Node &node) {
        node.size = 1 + nodeSize(node.left) + nodeSize(node.right);
        node.height = 1 + std::max(height(node.left), height(node.right));
        if constexpr (Aggregate)
            node.product = config_->op(config_->op(product(node.left), node.value), product(node.right));
    }
    void allApply(Ptr &node, F action) const {
        if (!node)
            return;
        node->value = config_->mapping(action, node->value);
        if constexpr (Aggregate)
            node->product = config_->mapping(action, node->product);
        node->lazy = node->pending ? config_->composition(action, node->lazy) : action;
        node->pending = true;
    }
    void push(Node &node) const {
        if constexpr (Lazy) {
            if (node.pending) {
                allApply(node.left, node.lazy);
                allApply(node.right, node.lazy);
                node.lazy = config_->id;
                node.pending = false;
            }
        }
    }
    Ptr rotateLeft(Ptr node) {
        push(*node);
        auto child = std::move(node->right);
        push(*child);
        node->right = std::move(child->left);
        pull(*node);
        child->left = std::move(node);
        pull(*child);
        return child;
    }
    Ptr rotateRight(Ptr node) {
        push(*node);
        auto child = std::move(node->left);
        push(*child);
        node->left = std::move(child->right);
        pull(*node);
        child->right = std::move(node);
        pull(*child);
        return child;
    }
    Ptr rebalance(Ptr node) {
        pull(*node);
        if (height(node->left) > height(node->right) + 1) {
            push(*node->left);
            if (height(node->left->left) < height(node->left->right))
                node->left = rotateLeft(std::move(node->left));
            return rotateRight(std::move(node));
        }
        if (height(node->right) > height(node->left) + 1) {
            push(*node->right);
            if (height(node->right->right) < height(node->right->left))
                node->right = rotateRight(std::move(node->right));
            return rotateLeft(std::move(node));
        }
        return node;
    }
    Ptr join(Ptr left, Ptr pivot, Ptr right) {
        if (height(left) > height(right) + 1) {
            push(*left);
            left->right = join(std::move(left->right), std::move(pivot), std::move(right));
            return rebalance(std::move(left));
        }
        if (height(right) > height(left) + 1) {
            push(*right);
            right->left = join(std::move(left), std::move(pivot), std::move(right->left));
            return rebalance(std::move(right));
        }
        pivot->left = std::move(left);
        pivot->right = std::move(right);
        pull(*pivot);
        return pivot;
    }
    std::pair<Ptr, Ptr> splitRoot(Ptr node, Int k) {
        if (!node)
            return {};
        push(*node);
        Int leftSize = nodeSize(node->left);
        auto left = std::move(node->left), right = std::move(node->right);
        if (k <= leftSize) {
            auto [first, middle] = splitRoot(std::move(left), k);
            return {std::move(first), join(std::move(middle), std::move(node), std::move(right))};
        }
        auto [middle, last] = splitRoot(std::move(right), k - leftSize - 1);
        return {join(std::move(left), std::move(node), std::move(middle)), std::move(last)};
    }
    Ptr mergeRoot(Ptr left, Ptr right) {
        if (!left)
            return right;
        if (!right)
            return left;
        Int length = nodeSize(left);
        auto [rest, pivot] = splitRoot(std::move(left), length - 1);
        return join(std::move(rest), std::move(pivot), std::move(right));
    }
    Ptr build(const std::vector<S> &values, Int l, Int r) {
        if (l == r)
            return {};
        Int mid = l + (r - l) / 2;
        auto node = std::make_unique<Node>(values[mid], config_->id);
        node->left = build(values, l, mid);
        node->right = build(values, mid + 1, r);
        pull(*node);
        return node;
    }
    static Ptr clone(const Ptr &node) {
        if (!node)
            return {};
        auto copy = std::make_unique<Node>(node->value, node->lazy);
        copy->height = node->height;
        copy->size = node->size;
        copy->product = node->product;
        copy->pending = node->pending;
        copy->left = clone(node->left);
        copy->right = clone(node->right);
        return copy;
    }
    void checkRange(Int l, Int r) const {
        assert(0 <= l && l <= r && r <= size());
    }
    S foldNode(const Ptr &node, Int l, Int r) const {
        if (!node || l == r)
            return config_->e;
        if (l == 0 && r == node->size)
            return node->product;
        push(*node);
        Int leftSize = nodeSize(node->left);
        S result = config_->e;
        if (l < leftSize)
            result = foldNode(node->left, l, std::min(r, leftSize));
        if (l <= leftSize && leftSize < r)
            result = config_->op(result, node->value);
        if (r > leftSize + 1)
            result = config_->op(result, foldNode(node->right, std::max<Int>(0, l - leftSize - 1), r - leftSize - 1));
        return result;
    }
    template <class Predicate>
    Int maxRightNode(const Ptr &node, Int offset, Int l, S &acc, Predicate &pred) const {
        if (!node || offset + node->size <= l)
            return -1;
        if (l <= offset) {
            S next = config_->op(acc, node->product);
            if (pred(next)) {
                acc = std::move(next);
                return -1;
            }
        }
        push(*node);
        Int mid = offset + nodeSize(node->left);
        Int result = maxRightNode(node->left, offset, l, acc, pred);
        if (result >= 0)
            return result;
        if (mid >= l) {
            S next = config_->op(acc, node->value);
            if (!pred(next))
                return mid;
            acc = std::move(next);
        }
        return maxRightNode(node->right, mid + 1, l, acc, pred);
    }
    template <class Predicate>
    Int minLeftNode(const Ptr &node, Int offset, Int r, S &acc, Predicate &pred) const {
        if (!node || r <= offset)
            return -1;
        if (offset + node->size <= r) {
            S next = config_->op(node->product, acc);
            if (pred(next)) {
                acc = std::move(next);
                return -1;
            }
        }
        push(*node);
        Int mid = offset + nodeSize(node->left);
        Int result = minLeftNode(node->right, mid + 1, r, acc, pred);
        if (result >= 0)
            return result;
        if (mid < r) {
            S next = config_->op(node->value, acc);
            if (!pred(next))
                return mid + 1;
            acc = std::move(next);
        }
        return minLeftNode(node->left, offset, r, acc, pred);
    }
    explicit AvlSequence(std::shared_ptr<const Config> config) : config_(std::move(config)) {
    }
public:
    // 配列からAVL列を構築する。演算・作用がO(1)なら時間・領域 O(N)。
    AvlSequence(const std::vector<S> &values, std::function<S(S, S)> op, S e,
                std::function<S(F, S)> mapping, std::function<F(F, F)> composition, F id)
        : config_(std::make_shared<Config>(Config{op, e, mapping, composition, id})) {
        root_ = build(values, 0, values.size());
    }
    // 独立した木を複製する。時間・領域 O(N)。
    AvlSequence(const AvlSequence &other) : root_(clone(other.root_)), config_(other.config_) {
    }
    // 木を複製して代入する。新旧の要素数の和をNとして時間 O(N)、領域 O(新要素数)。
    AvlSequence &operator=(const AvlSequence &other) {
        if (this != &other) {
            auto copy = clone(other.root_);
            root_ = std::move(copy);
            config_ = other.config_;
        }
        return *this;
    }
    // 木の所有権を移す。計算量 O(1)。
    AvlSequence(AvlSequence &&) noexcept = default;
    // 木の所有権を移す。旧要素の破棄に O(N)。
    AvlSequence &operator=(AvlSequence &&) noexcept = default;
    // 要素数を返す。計算量 O(1)。
    Int size() const {
        return nodeSize(root_);
    }
    // 先頭k要素で分割し、元を空にする。時間 O(log(N+1))。
    auto split(Int k) {
        checkRange(0, k);
        auto [left, right] = splitRoot(std::move(root_), k);
        AvlSequence first(config_), second(config_);
        first.root_ = std::move(left);
        second.root_ = std::move(right);
        return std::pair{std::move(first), std::move(second)};
    }
    // 同じ演算・単位元を持つ列を連結し、元を両方空にする。時間 O(log(N+1))。
    static AvlSequence merge(AvlSequence &left, AvlSequence &right) {
        assert(&left != &right);
        AvlSequence result(left.config_);
        result.root_ = result.mergeRoot(std::move(left.root_), std::move(right.root_));
        return result;
    }
    // index直前に挿入する。過去の作用は新要素に適用しない。時間 O(log(N+1))。
    void insert(Int index, S value) {
        checkRange(0, index);
        auto [left, right] = splitRoot(std::move(root_), index);
        root_ = join(std::move(left), std::make_unique<Node>(value, config_->id), std::move(right));
    }
    // ソート済み列の同値要素の末尾に挿入し位置を返す。時間 O(log(N+1))。
    template <class Less = std::less<S>> Int insort(S value, Less less = {}) {
        Int index = 0;
        Node *node = root_.get();
        while (node) {
            push(*node);
            if (less(value, node->value)) {
                node = node->left.get();
            } else {
                index += nodeSize(node->left) + 1;
                node = node->right.get();
            }
        }
        insert(index, value);
        return index;
    }
    // 半開区間を削除する。木の操作 O(log(N+1))、K要素の破棄 O(K)。
    void erase(Int l, Int r) {
        checkRange(l, r);
        auto [left, rest] = splitRoot(std::move(root_), l);
        auto [removed, right] = splitRoot(std::move(rest), r - l);
        root_ = mergeRoot(std::move(left), std::move(right));
    }
    // 一要素を削除する。時間 O(log N)。
    void erase(Int index) {
        erase(index, index + 1);
    }
    // 指定位置の要素を返す。時間 O(log N)。
    S get(Int index) const {
        assert(0 <= index && index < size());
        Node *node = root_.get();
        while (true) {
            push(*node);
            Int leftSize = nodeSize(node->left);
            if (index < leftSize)
                node = node->left.get();
            else if (index == leftSize)
                return node->value;
            else {
                index -= leftSize + 1;
                node = node->right.get();
            }
        }
    }
    struct Reference {
        AvlSequence *owner;
        Int index;
        // 対象要素を読み出す。時間 O(log N)。
        operator S() const {
            return owner->get(index);
        }
        // 対象要素を上書きする。時間 O(log N)。
        Reference &operator=(S value) {
            owner->update(index, std::move(value));
            return *this;
        }
        // 他の位置を読み出して代入する。時間 O(log N)。
        Reference &operator=(const Reference &other) {
            return *this = S(other);
        }
    };
    // 指定位置を読み書きする代理参照を返す。取得 O(1)、読み書き O(log N)。
    Reference operator[](Int index) {
        assert(0 <= index && index < size());
        return {this, index};
    }
    // 末尾からの位置を読み書きする代理参照を返す。取得 O(1)、読み書き O(log N)。
    Reference operator[](BackwardsIndex index) {
        return (*this)[size() - index.value];
    }
    // 指定位置の要素を返す。時間 O(log N)。
    S operator[](Int index) const {
        return get(index);
    }
    // 末尾から数えた位置の要素を返す。時間 O(log N)。
    S operator[](BackwardsIndex index) const {
        return get(size() - index.value);
    }
    // 指定位置を上書きする。時間 O(log N)。
    void update(Int index, S value) {
        checkRange(index, index + 1);
        auto [left, rest] = splitRoot(std::move(root_), index);
        auto [middle, right] = splitRoot(std::move(rest), 1);
        middle->value = value;
        pull(*middle);
        root_ = join(std::move(left), std::move(middle), std::move(right));
    }
    // 半開区間に作用を適用する。時間 O(log(N+1))。
    void apply(Int l, Int r, F action) requires Lazy {
        checkRange(l, r);
        auto [left, rest] = splitRoot(std::move(root_), l);
        auto [middle, right] = splitRoot(std::move(rest), r - l);
        allApply(middle, action);
        root_ = mergeRoot(std::move(left), mergeRoot(std::move(middle), std::move(right)));
    }
    // 一要素に作用を適用する。時間 O(log N)。
    void apply(Int index, F action) requires Lazy {
        apply(index, index + 1, action);
    }
    // 半開区間の積を順序を保って返す。時間 O(log(N+1))。
    S fold(Int l, Int r) const requires Aggregate {
        checkRange(l, r);
        return foldNode(root_, l, r);
    }
    // 全要素の積を返す。計算量 O(1)。
    S fold() const requires Aggregate {
        return product(root_);
    }
    // 半開区間の積を返す。時間 O(log(N+1))。
    S get(Int l, Int r) const requires Aggregate {
        return fold(l, r);
    }
    // 全要素の積を返す。計算量 O(1)。
    S get_all() const requires Aggregate {
        return fold();
    }
    // predが真の最大右端を返す。pred(e)=trueで右端について単調であること。時間 O(log(N+1))。
    template <class Predicate> Int max_right(Int l, Predicate pred) const requires Aggregate {
        checkRange(l, size());
        S acc = config_->e;
        assert(pred(acc));
        Int result = maxRightNode(root_, 0, l, acc, pred);
        return result < 0 ? size() : result;
    }
    // predが真の最小左端を返す。pred(e)=trueで左端について単調であること。時間 O(log(N+1))。
    template <class Predicate> Int min_left(Int r, Predicate pred) const requires Aggregate {
        checkRange(0, r);
        S acc = config_->e;
        assert(pred(acc));
        Int result = minLeftNode(root_, 0, r, acc, pred);
        return result < 0 ? 0 : result;
    }
    // 末尾基準で指定した位置を上書きする。時間 O(log N)。
    void update(BackwardsIndex index, S value) {
        update(size() - index.value, std::move(value));
    }
    // 両端を含む区間を削除する。木の操作 O(log(N+1))、K要素の破棄 O(K)。
    template <class L, class R> void erase(ClosedSlice<L, R> range) {
        erase(resolve_index(size(), range.a), resolve_index(size(), range.b) + 1);
    }
    // 両端を含む区間に作用を適用する。時間 O(log(N+1))。
    template <class L, class R> void apply(ClosedSlice<L, R> range, F action) requires Lazy {
        apply(resolve_index(size(), range.a), resolve_index(size(), range.b) + 1, action);
    }
    // 両端を含む区間の積を返す。時間 O(log(N+1))。
    template <class L, class R> S get(ClosedSlice<L, R> range) const requires Aggregate {
        return fold(resolve_index(size(), range.a), resolve_index(size(), range.b) + 1);
    }
    // 両端を含む区間の積を返す。時間 O(log(N+1))。
    template <class L, class R> S fold(ClosedSlice<L, R> range) const requires Aggregate {
        return get(range);
    }
    // 両端を含む区間の積を返す。時間 O(log(N+1))。
    template <class L, class R> S operator[](ClosedSlice<L, R> range) const requires Aggregate {
        return get(range);
    }
    class Iterator {
        const AvlSequence *owner_ = nullptr;
        std::vector<Node *> stack_;
        void descend(Node *node) {
            while (node) {
                owner_->push(*node);
                stack_.push_back(node);
                node = node->left.get();
            }
        }
    public:
        using iterator_category = std::input_iterator_tag;
        using value_type = S;
        using difference_type = std::ptrdiff_t;
        using reference = const S &;
        // 終端イテレータを作る。計算量 O(1)。
        Iterator() = default;
        // 指定部分木の先頭へ進める。時間・領域 O(log(N+1))。
        Iterator(const AvlSequence *owner, Node *node) : owner_(owner) {
            descend(node);
        }
        // 現在の要素を参照する。計算量 O(1)。
        const S &operator*() const {
            return stack_.back()->value;
        }
        // 次の要素へ進む。一回の最悪 O(log(N+1))、全列挙 O(N)。
        Iterator &operator++() {
            Node *node = stack_.back();
            stack_.pop_back();
            descend(node->right.get());
            return *this;
        }
        // 更新前のイテレータを返して次へ進む。スタックのコピーを含め O(log(N+1))。
        Iterator operator++(int) {
            auto previous = *this;
            ++*this;
            return previous;
        }
        // 同じ走査位置か返す。計算量 O(1)。
        bool operator==(const Iterator &other) const {
            if (stack_.empty() || other.stack_.empty())
                return stack_.empty() == other.stack_.empty();
            return owner_ == other.owner_ && stack_.back() == other.stack_.back();
        }
    };
    // 先頭イテレータを返す。時間・領域 O(log(N+1))。走査中の構造変更は禁止する。
    Iterator begin() const {
        return Iterator(this, root_.get());
    }
    // 終端イテレータを返す。計算量 O(1)。
    Iterator end() const {
        return Iterator();
    }
    // 要素を列挙する範囲を返す。取得 O(1)、全列挙 O(N)、走査用領域 O(log(N+1))。
    const AvlSequence &items() const {
        return *this;
    }
    // 全要素を文字列化する。時間 O(N+出力長)と要素の変換費用、領域 O(出力長+log(N+1))。
    std::string to_string() const {
        std::ostringstream output;
        output << "@[";
        bool first = true;
        for (const S &value : *this) {
            if (!first)
                output << ", ";
            first = false;
            output << value;
        }
        output << ']';
        return output.str();
    }
    // 全要素を順番にコピーする。時間・出力領域 O(N)。
    std::vector<S> toSeq() const {
        std::vector<S> result;
        result.reserve(size());
        auto visit = [&](auto &&self, const Ptr &node) -> void {
            if (!node)
                return;
            push(*node);
            self(self, node->left);
            result.push_back(node->value);
            self(self, node->right);
        };
        visit(visit, root_);
        return result;
    }
};
// 同じ設定のAVL列を破壊的に結合する。時間 O(log(N+1))。
template <class S, class F, bool A, bool L>
auto merge(AvlSequence<S, F, A, L> &left, AvlSequence<S, F, A, L> &right) {
    return AvlSequence<S, F, A, L>::merge(left, right);
}
}
