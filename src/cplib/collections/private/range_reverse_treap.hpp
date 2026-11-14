#pragma once
#include <cplib/utils/backwards_index.hpp>
#include <memory>
#include <random>
#include <sstream>

namespace cplib::detail {
struct NoTreapAction {};

template <class S, bool Enabled> struct TreapAggregate {
    explicit TreapAggregate(const S &) {};
};

template <class S> struct TreapAggregate<S, true> {
    S prod, rprod;

    explicit TreapAggregate(const S &value) : prod(value), rprod(value) {
    }
};

template <class F, bool Enabled> struct TreapLazy {
    explicit TreapLazy(const F &) {};
};

template <class F> struct TreapLazy<F, true> {
    F lazy;

    explicit TreapLazy(const F &id) : lazy(id) {
    }
};

inline std::mt19937_64 range_reverse_random(std::random_device{}());

// 通常・モノイド・双対・遅延の4種類は別の節点型を生成する。
// モノイド版の直接探索と挿入、遅延版のsplit/mergeによる操作をそれぞれ保持。
template <class S, class F, bool Aggregate, bool Lazy> class RangeReverseTreap {
    struct Node;
    using Ptr = std::shared_ptr<Node>;

    struct Node : TreapAggregate<S, Aggregate>, TreapLazy<F, Lazy> {
        Ptr left, right;
        UInt priority;
        Int size = 1;
        bool rev = false;
        S value;

        Node(const S &value, UInt priority, const F &id)
            : TreapAggregate<S, Aggregate>(value), TreapLazy<F, Lazy>(id), priority(priority),
              value(value) {
        }
    };

    Ptr root_;
    Int length_ = 0;
    std::function<S(S, S)> op_;
    S identity_;
    std::function<S(F, S)> mapping_;
    std::function<F(F, F)> composition_;
    F id_;

    static Int size(const Ptr &node) {
        return node ? node->size : 0;
    }

    S product(const Ptr &node) const
        requires(Aggregate)
    {
        return node ? node->prod : identity_;
    }

    S reverse_product(const Ptr &node) const
        requires(Aggregate)
    {
        return node ? node->rprod : identity_;
    }

    void pull(const Ptr &node) {
        if (!node)
            return;
        node->size = 1 + size(node->left) + size(node->right);
        if constexpr (Aggregate) {
            node->prod = op_(op_(product(node->left), node->value), product(node->right));
            node->rprod =
                op_(op_(reverse_product(node->right), node->value), reverse_product(node->left));
        }
    }

    static void toggle(const Ptr &node) {
        if (!node)
            return;
        node->rev = !node->rev;
        if constexpr (Aggregate)
            std::swap(node->prod, node->rprod);
    }

    void all_apply(const Ptr &node, const F &action)
        requires(Lazy)
    {
        if (!node)
            return;
        node->value = mapping_(action, node->value);
        if constexpr (Aggregate) {
            node->prod = mapping_(action, node->prod);
            node->rprod = mapping_(action, node->rprod);
        }
        node->lazy = composition_(action, node->lazy);
    }

    void push(const Ptr &node) {
        if (!node)
            return;
        if (node->rev) {
            std::swap(node->left, node->right);
            toggle(node->left);
            toggle(node->right);
            node->rev = false;
        }
        if constexpr (Lazy) {
            all_apply(node->left, node->lazy);
            all_apply(node->right, node->lazy);
            node->lazy = id_;
        }
    }

    Ptr make_node(const S &value) {
        return std::make_shared<Node>(value, range_reverse_random(), id_);
    }

    void update_all(const Ptr &node) {
        if (!node)
            return;
        update_all(node->left);
        update_all(node->right);
        pull(node);
    }

    template <class V> void build(const V &values) {
        std::vector<Ptr> stack;
        for (const S &value : values) {
            auto node = make_node(value);
            Ptr last;
            while (!stack.empty() && stack.back()->priority < node->priority) {
                last = stack.back();
                stack.pop_back();
            }
            node->left = last;
            if (!stack.empty())
                stack.back()->right = node;
            stack.push_back(node);
        }
        if (!stack.empty()) {
            root_ = stack.front();
            update_all(root_);
        }
    }

    Ptr merge(Ptr left, Ptr right) {
        if (!left)
            return right;
        if (!right)
            return left;
        if (left->priority > right->priority) {
            push(left);
            left->right = merge(left->right, right);
            pull(left);
            return left;
        }
        push(right);
        right->left = merge(left, right->left);
        pull(right);
        return right;
    }

    std::pair<Ptr, Ptr> split(Ptr node, Int k) {
        if (!node)
            return {};
        push(node);
        Int left = size(node->left);
        if (k <= left) {
            auto [a, b] = split(node->left, k);
            node->left = b;
            pull(node);
            return {a, node};
        }
        auto [a, b] = split(node->right, k - left - 1);
        node->right = a;
        pull(node);
        return {node, b};
    }

    Ptr insert_node(Ptr root, Ptr node, Int k)
        requires(!Lazy)
    {
        if (!root)
            return node;
        if (node->priority > root->priority) {
            std::tie(node->left, node->right) = split(root, k);
            pull(node);
            return node;
        }
        push(root);
        Int left = size(root->left);
        if (k <= left)
            root->left = insert_node(root->left, node, k);
        else
            root->right = insert_node(root->right, node, k - left - 1);
        pull(root);
        return root;
    }

    Ptr erase_node(Ptr node, Int k)
        requires(!Lazy)
    {
        push(node);
        Int left = size(node->left);
        if (k == left)
            return merge(node->left, node->right);
        if (k < left)
            node->left = erase_node(node->left, k);
        else
            node->right = erase_node(node->right, k - left - 1);
        pull(node);
        return node;
    }

    S get_node(const Ptr &node, Int l, Int r)
        requires(Aggregate && !Lazy)
    {
        if (l == 0 && r == node->size)
            return node->prod;
        push(node);
        Int mid = size(node->left);
        if (r <= mid)
            return get_node(node->left, l, r);
        if (l > mid)
            return get_node(node->right, l - mid - 1, r - mid - 1);
        S result = node->value;
        if (l < mid)
            result = op_(get_node(node->left, l, mid), result);
        if (r > mid + 1)
            result = op_(result, get_node(node->right, 0, r - mid - 1));
        return result;
    }

    void update_node(const Ptr &node, Int k, const S &value)
        requires(Aggregate && !Lazy)
    {
        push(node);
        Int left = size(node->left);
        if (k == left)
            node->value = value;
        else if (k < left)
            update_node(node->left, k, value);
        else
            update_node(node->right, k - left - 1, value);
        pull(node);
    }

    template <class Pred>
    Int search_right(const Ptr &node, Int start, Int l, S &acc, Pred &predicate)
        requires(Aggregate && !Lazy)
    {
        Int finish = start + size(node);
        if (!node || finish <= l)
            return finish;
        if (l <= start) {
            S next = op_(acc, node->prod);
            if (predicate(next)) {
                acc = next;
                return finish;
            }
        }
        push(node);
        Int mid = start + size(node->left),
            result = search_right(node->left, start, l, acc, predicate);
        if (result < mid)
            return result;
        if (l <= mid) {
            S next = op_(acc, node->value);
            if (!predicate(next))
                return mid;
            acc = next;
        }
        return search_right(node->right, mid + 1, l, acc, predicate);
    }

    template <class Pred>
    Int search_left(const Ptr &node, Int start, Int r, S &acc, Pred &predicate)
        requires(Aggregate && !Lazy)
    {
        if (!node || r <= start)
            return start;
        if (start + node->size <= r) {
            S next = op_(node->prod, acc);
            if (predicate(next)) {
                acc = next;
                return start;
            }
        }
        push(node);
        Int mid = start + size(node->left),
            result = search_left(node->right, mid + 1, r, acc, predicate);
        if (result > mid + 1)
            return result;
        if (mid < r) {
            S next = op_(node->value, acc);
            if (!predicate(next))
                return mid + 1;
            acc = next;
        }
        return search_left(node->left, start, r, acc, predicate);
    }

    static Ptr clone(const Ptr &node) {
        if (!node)
            return {};
        auto out = std::make_shared<Node>(*node);
        out->left = clone(node->left);
        out->right = clone(node->right);
        return out;
    }

    void check_range(Int l, Int r) const {
        assert(0 <= l && l <= r && r <= length_);
    }

public:
    // Cartesian treeの単調スタックでO(N)構築。各操作は期待O(log N)。
    template <class V>
    RangeReverseTreap(const V &values, std::function<S(S, S)> op, S identity,
                      std::function<S(F, S)> mapping, std::function<F(F, F)> composition, F id)
        : length_(values.size()), op_(op), identity_(identity), mapping_(mapping),
          composition_(composition), id_(id) {
        build(values);
    }

    // 全ノードを複製する。要素数をnとしてO(n)。
    RangeReverseTreap(const RangeReverseTreap &other)
        : root_(clone(other.root_)), length_(other.length_), op_(other.op_),
          identity_(other.identity_), mapping_(other.mapping_), composition_(other.composition_),
          id_(other.id_) {
    }

    // 列をムーブ構築する。O(1)。
    RangeReverseTreap(RangeReverseTreap &&) = default;
    // 列をムーブ代入する。代入先の要素数をnとしてO(n)。
    RangeReverseTreap &operator=(RangeReverseTreap &&) = default;

    // 列を複製して代入する。両者の要素数の合計をnとしてO(n)。
    RangeReverseTreap &operator=(const RangeReverseTreap &other) {
        if (this != &other) {
            RangeReverseTreap copy(other);
            *this = std::move(copy);
        }
        return *this;
    }

    // 要素数を返す。O(1)。
    Int size() const {
        return length_;
    }

    // 要素数をnとして期待O(log(n+1))。
    // 半開区間[l, r)を反転します。
    void reverse(Int l, Int r) {
        check_range(l, r);
        auto [left, rest] = split(root_, l);
        auto [middle, right] = split(rest, r - l);
        toggle(middle);
        root_ = merge(left, merge(middle, right));
    }

    // 要素数をnとして期待O(log(n+1))。
    // 閉区間segmentを反転します。
    template <class L, class R> void reverse(ClosedSlice<L, R> range) {
        reverse(resolve_index(length_, range.a), resolve_index(length_, range.b) + 1);
    }

    // index番目の値を取得する。要素数をnとして期待O(log(n+1))。
    S get(Int index) {
        assert(0 <= index && index < length_);
        auto node = root_;
        for (;;) {
            push(node);
            Int left = size(node->left);
            if (index < left)
                node = node->left;
            else if (index == left)
                return node->value;
            else {
                index -= left + 1;
                node = node->right;
            }
        }
    }

    // index番目の値を設定する。要素数をnとして期待O(log(n+1))。
    void update(Int index, const S &value) {
        assert(0 <= index && index < length_);
        if constexpr (Lazy) {
            auto [left, rest] = split(root_, index);
            auto [middle, right] = split(rest, 1);
            middle->value = value;
            middle->lazy = id_;
            middle->rev = false;
            if constexpr (Aggregate)
                middle->prod = middle->rprod = value;
            root_ = merge(left, merge(middle, right));
        } else if constexpr (Aggregate)
            update_node(root_, index, value);
        else {
            auto node = root_;
            for (;;) {
                push(node);
                Int left = size(node->left);
                if (index < left)
                    node = node->left;
                else if (index == left) {
                    node->value = value;
                    return;
                } else {
                    index -= left + 1;
                    node = node->right;
                }
            }
        }
    }

    struct Reference {
        RangeReverseTreap *owner;
        Int index;

        // 対象の値を取得する。要素数をnとして期待O(log(n+1))。
        operator S() const {
            return owner->get(index);
        }

        // 対象の値を設定する。要素数をnとして期待O(log(n+1))。
        Reference &operator=(const S &value) {
            owner->update(index, value);
            return *this;
        }

        // 他の位置の値を読み取り代入する。要素数をnとして期待O(log(n+1))。
        Reference &operator=(const Reference &value) {
            return *this = S(value);
        }
    };

    // i番目を読み書きする参照を返す。O(1)。
    Reference operator[](Int i) {
        return {this, i};
    }

    // 末尾からi番目を読み書きする参照を返す。O(1)。
    Reference operator[](BackwardsIndex i) {
        return (*this)[length_ - i.value];
    }

    // 先頭k要素で分割し、元の列を空にする。期待時間 O(log(N+1))。
    auto split(Int k) {
        check_range(0, k);
        RangeReverseTreap left(std::vector<S>{}, op_, identity_, mapping_, composition_, id_);
        RangeReverseTreap right(std::vector<S>{}, op_, identity_, mapping_, composition_, id_);
        auto roots = split(root_, k);
        left.root_ = std::move(roots.first);
        right.root_ = std::move(roots.second);
        left.length_ = k;
        right.length_ = length_ - k;
        root_.reset();
        length_ = 0;
        return std::pair{std::move(left), std::move(right)};
    }
    // 同じ演算・作用の列を結合し、両方を空にする。期待時間 O(log(N+1))。
    static RangeReverseTreap merge(RangeReverseTreap &left, RangeReverseTreap &right) {
        assert(&left != &right);
        RangeReverseTreap result(std::vector<S>{}, left.op_, left.identity_, left.mapping_, left.composition_, left.id_);
        result.root_ = result.merge(left.root_, right.root_);
        result.length_ = left.length_ + right.length_;
        left.root_.reset();
        right.root_.reset();
        left.length_ = right.length_ = 0;
        return result;
    }
    // 昇順の列で同値の末尾に挿入して位置を返す。期待時間 O(log(N+1))。
    template <class Less = std::less<S>> Int insort(const S &value, Less less = {}) {
        Int index = 0;
        auto node = root_;
        while (node) {
            push(node);
            if (less(value, node->value)) {
                node = node->left;
            } else {
                index += size(node->left) + 1;
                node = node->right;
            }
        }
        insert(index, value);
        return index;
    }
    // indexの直前にvalueを挿入します。index=size()なら末尾です。期待O(log N)。
    // 遅延版では挿入前の区間更新は新しい要素には作用しません。
    void insert(Int index, const S &value)
    {
        assert(0 <= index && index <= length_);
        auto node = make_node(value);
        if constexpr (!Lazy)
            root_ = insert_node(root_, node, index);
        else {
            auto [left, right] = split(root_, index);
            root_ = merge(left, merge(node, right));
        }
        ++length_;
    }

    // index番目の要素を削除する。要素数をnとして期待O(log(n+1))。
    void erase(Int index)
    {
        assert(0 <= index && index < length_);
        if constexpr (!Lazy) {
            root_ = erase_node(root_, index);
            --length_;
        } else
            erase(index, index + 1);
    }

    // 半開区間 [l, r) を削除する。木の操作は期待 O(log N)。
    // 削除した K 個のノードの解放には別途 O(K) かかりうる。
    void erase(Int l, Int r)
    {
        check_range(l, r);
        if (l == r)
            return;
        auto [left, rest] = split(root_, l);
        auto [middle, right] = split(rest, r - l);
        root_ = merge(left, right);
        length_ -= r - l;
    }

    // 両端を含む区間を削除する。要素数n、削除数kとして期待O(log(n+1)+k)。
    template <class L, class R>
    void erase(ClosedSlice<L, R> range)
    {
        erase(resolve_index(length_, range.a), resolve_index(length_, range.b) + 1);
    }

    // 要素数をnとして期待O(log(n+1))。
    // 半開区間 [l, r) の総積を返す。
    S get(Int l, Int r)
        requires(Aggregate)
    {
        check_range(l, r);
        if constexpr (!Lazy) {
            if (l == r)
                return identity_;
            return get_node(root_, l, r);
        } else {
            auto [left, rest] = split(root_, l);
            auto [middle, right] = split(rest, r - l);
            S result = product(middle);
            root_ = merge(left, merge(middle, right));
            return result;
        }
    }

    // 要素数をnとして期待O(log(n+1))。
    // 閉区間segmentの総積を返します。
    template <class L, class R>
    S get(ClosedSlice<L, R> range)
        requires(Aggregate)
    {
        return get(resolve_index(length_, range.a), resolve_index(length_, range.b) + 1);
    }

    // 両端を含む区間の総積を返す。要素数をnとして期待O(log(n+1))。
    template <class L, class R>
    S operator[](ClosedSlice<L, R> range)
        requires(Aggregate)
    {
        return get(range);
    }

    // 全要素の総積を返す。O(1)。
    S get_all() const
        requires(Aggregate)
    {
        return product(root_);
    }

    // 全要素の総積を返す。O(1)。
    S fold() const
        requires(Aggregate)
    {
        return get_all();
    }

    // 半開区間[l,r)の総積を返す。要素数をnとして期待O(log(n+1))。
    S fold(Int l, Int r)
        requires(Aggregate)
    {
        return get(l, r);
    }

    // 両端を含む区間の総積を返す。要素数をnとして期待O(log(n+1))。
    template <class L, class R>
    S fold(ClosedSlice<L, R> range)
        requires(Aggregate)
    {
        return get(range);
    }

    // 半開区間[l,r)にactionを作用させる。要素数をnとして期待O(log(n+1))。
    void apply(Int l, Int r, const F &action)
        requires(Lazy)
    {
        check_range(l, r);
        auto [left, rest] = split(root_, l);
        auto [middle, right] = split(rest, r - l);
        all_apply(middle, action);
        root_ = merge(left, merge(middle, right));
    }

    // index番目にactionを作用させる。要素数をnとして期待O(log(n+1))。
    void apply(Int index, const F &action)
        requires(Lazy)
    {
        assert(0 <= index && index < length_);
        apply(index, index + 1, action);
    }

    // 両端を含む区間にactionを作用させる。要素数をnとして期待O(log(n+1))。
    template <class L, class R>
    void apply(ClosedSlice<L, R> range, const F &action)
        requires(Lazy)
    {
        apply(resolve_index(length_, range.a), resolve_index(length_, range.b) + 1, action);
    }

    // 要素数をnとして期待O(log(n+1))回の演算と述語呼出し。
    // predicate(get(l,r))が真となる最大のrを返します。
    // predicate(identity_)=trueで、区間を伸ばしたとき真から偽への変化が単調であること。
    template <class Pred>
    Int max_right(Int l, Pred predicate)
        requires(Aggregate)
    {
        assert(0 <= l && l <= length_ && predicate(identity_));
        S acc = identity_;
        if constexpr (!Lazy)
            return search_right(root_, 0, l, acc, predicate);
        else {
            auto [left, right] = split(root_, l);
            auto node = right;
            Int result = l;
            while (node) {
                push(node);
                S next = op_(acc, product(node->left));
                if (!predicate(next))
                    node = node->left;
                else {
                    acc = next;
                    result += size(node->left);
                    S value = op_(acc, node->value);
                    if (!predicate(value))
                        break;
                    acc = value;
                    ++result;
                    node = node->right;
                }
            }
            root_ = merge(left, right);
            return result;
        }
    }

    // 要素数をnとして期待O(log(n+1))回の演算と述語呼出し。
    // predicate(get(l,r))が真となる最小のlを返します。
    // predicate(identity_)=trueで、区間を伸ばしたとき真から偽への変化が単調であること。
    template <class Pred>
    Int min_left(Int r, Pred predicate)
        requires(Aggregate)
    {
        assert(0 <= r && r <= length_ && predicate(identity_));
        S acc = identity_;
        if constexpr (!Lazy)
            return search_left(root_, 0, r, acc, predicate);
        else {
            auto [left, right] = split(root_, r);
            auto node = left;
            Int result = r;
            while (node) {
                push(node);
                S next = op_(product(node->right), acc);
                if (!predicate(next))
                    node = node->right;
                else {
                    acc = next;
                    result -= size(node->right);
                    S value = op_(node->value, acc);
                    if (!predicate(value))
                        break;
                    acc = value;
                    --result;
                    node = node->left;
                }
            }
            root_ = merge(left, right);
            return result;
        }
    }

    struct Iterator {
        RangeReverseTreap *owner;
        std::vector<Ptr> stack;

        // nodeから左端まで遅延を伝播し、走査用スタックに積む。訪問ノード数をkとしてO(k)。
        void descend(Ptr node) {
            while (node) {
                owner->push(node);
                stack.push_back(node);
                node = node->left;
            }
        }

        // 現在の要素を返す。O(1)。
        S operator*() const {
            return stack.back()->value;
        }

        // 次の要素に進む。期待O(log(n+1))、全要素の走査では合計O(n)。
        Iterator &operator++() {
            auto node = stack.back();
            stack.pop_back();
            descend(node->right);
            return *this;
        }

        // 同じ位置を指すか比較する。O(1)。
        bool operator==(const Iterator &other) const {
            return owner == other.owner &&
                   (stack.empty() ? other.stack.empty()
                                  : !other.stack.empty() && stack.back() == other.stack.back());
        }
    };

    // 先頭のイテレータを返す。要素数をnとして期待O(log(n+1))。
    Iterator begin() {
        Iterator it{this, {}};
        it.descend(root_);
        return it;
    }

    // 末尾の次のイテレータを返す。O(1)。
    Iterator end() {
        return {this, {}};
    }

    // 列挙用に自身への参照を返す。O(1)。
    RangeReverseTreap &items() {
        return *this;
    }

    // 全要素を配列にして返す。要素数をnとしてO(n)。
    std::vector<S> toSeq() {
        std::vector<S> result;
        result.reserve(length_);
        for (S x : *this)
            result.push_back(x);
        return result;
    }

    // 全要素を空白区切りで文字列化する。O(n+出力長)と各要素の文字列化の計算量。
    std::string str() {
        std::ostringstream out;
        bool first = true;
        for (S x : *this) {
            if (!first)
                out << ' ';
            first = false;
            out << x;
        }
        return out.str();
    }
};
}

namespace cplib {
#define CPLIB_REVERSE_TEMPLATE template <class S, class F, bool A, bool L>
#define CPLIB_REVERSE_TYPE detail::RangeReverseTreap<S, F, A, L>

// 先頭k要素で分割して元を空にする。期待時間 O(log(N+1))。
CPLIB_REVERSE_TEMPLATE auto split(CPLIB_REVERSE_TYPE &tree, Int k) {
    return tree.split(k);
}
// 同じ演算・作用の列を結合して両方を空にする。期待時間 O(log(N+1))。
CPLIB_REVERSE_TEMPLATE auto merge(CPLIB_REVERSE_TYPE &left, CPLIB_REVERSE_TYPE &right) {
    return CPLIB_REVERSE_TYPE::merge(left, right);
}
// 昇順の列へ同値の末尾に挿入し位置を返す。期待時間 O(log(N+1))。
CPLIB_REVERSE_TEMPLATE Int insort(CPLIB_REVERSE_TYPE &tree, const S &value) {
    return tree.insort(value);
}
// 半開区間[l,r)を反転する。要素数nとして期待O(log(n+1))。
CPLIB_REVERSE_TEMPLATE void reverse(CPLIB_REVERSE_TYPE &tree, Int l, Int r) {
    tree.reverse(l, r);
}

// 両端を含む区間を反転する。要素数nとして期待O(log(n+1))。
template <class S, class F, bool A, bool L, class I, class J>
void reverse(CPLIB_REVERSE_TYPE &tree, ClosedSlice<I, J> range) {
    tree.reverse(range);
}

// i番目の値を返す。要素数nとして期待O(log(n+1))。
CPLIB_REVERSE_TEMPLATE S get(CPLIB_REVERSE_TYPE &tree, Int i) {
    return tree.get(i);
}

// 半開区間[l,r)の総積を返す。要素数nとして期待O(log(n+1))。
CPLIB_REVERSE_TEMPLATE S get(CPLIB_REVERSE_TYPE &tree, Int l, Int r)
    requires(A)
{
    return tree.get(l, r);
}

// 両端を含む区間の総積を返す。要素数nとして期待O(log(n+1))。
template <class S, class F, bool A, bool L, class I, class J>
S get(CPLIB_REVERSE_TYPE &tree, ClosedSlice<I, J> range)
    requires(A)
{
    return tree.get(range);
}

// i番目の値を更新する。要素数nとして期待O(log(n+1))。
CPLIB_REVERSE_TEMPLATE void update(CPLIB_REVERSE_TYPE &tree, Int i, const S &value) {
    tree.update(i, value);
}

// i番目の直前に値を挿入する。要素数nとして期待O(log(n+1))。
CPLIB_REVERSE_TEMPLATE void insert(CPLIB_REVERSE_TYPE &tree, Int i, const S &value)
{
    tree.insert(i, value);
}

// i番目を削除する。要素数nとして期待O(log(n+1))。
CPLIB_REVERSE_TEMPLATE void erase(CPLIB_REVERSE_TYPE &tree, Int i)
{
    tree.erase(i);
}

// 半開区間[l,r)を削除する。要素数n、削除数kとして期待O(log(n+1)+k)。
CPLIB_REVERSE_TEMPLATE void erase(CPLIB_REVERSE_TYPE &tree, Int l, Int r)
{
    tree.erase(l, r);
}

// 両端を含む区間を削除する。要素数nとして期待O(log(n+1)+k)、kは削除数。
template <class S, class F, bool A, bool L, class I, class J>
void erase(CPLIB_REVERSE_TYPE &tree, ClosedSlice<I, J> range)
{
    tree.erase(range);
}

// 列挙用に自身への参照を返す。O(1)。
CPLIB_REVERSE_TEMPLATE CPLIB_REVERSE_TYPE &items(CPLIB_REVERSE_TYPE &tree) {
    return tree;
}

// 全要素を配列にして返す。要素数nとしてO(n)。
CPLIB_REVERSE_TEMPLATE auto toSeq(CPLIB_REVERSE_TYPE &tree) {
    return tree.toSeq();
}

// 全要素を文字列化する。O(n+出力長)と各要素の文字列化の計算量。
CPLIB_REVERSE_TEMPLATE std::string to_string(CPLIB_REVERSE_TYPE &tree) {
    return tree.str();
}

// 全要素の総積を返す。O(1)。
CPLIB_REVERSE_TEMPLATE S get_all(const CPLIB_REVERSE_TYPE &tree)
    requires(A)
{
    return tree.get_all();
}

// 全要素の総積を返す。O(1)。
CPLIB_REVERSE_TEMPLATE S fold(const CPLIB_REVERSE_TYPE &tree)
    requires(A)
{
    return tree.fold();
}

// 半開区間[l,r)の総積を返す。要素数nとして期待O(log(n+1))。
CPLIB_REVERSE_TEMPLATE S fold(CPLIB_REVERSE_TYPE &tree, Int l, Int r)
    requires(A)
{
    return tree.fold(l, r);
}

// 両端を含む区間の総積を返す。要素数nとして期待O(log(n+1))。
template <class S, class F, bool A, bool L, class I, class J>
S fold(CPLIB_REVERSE_TYPE &tree, ClosedSlice<I, J> range)
    requires(A)
{
    return tree.fold(range);
}

// i番目に作用を適用する。要素数nとして期待O(log(n+1))。
CPLIB_REVERSE_TEMPLATE void apply(CPLIB_REVERSE_TYPE &tree, Int i, const F &f)
    requires(L)
{
    tree.apply(i, f);
}

// 半開区間[l,r)に作用を適用する。要素数nとして期待O(log(n+1))。
CPLIB_REVERSE_TEMPLATE void apply(CPLIB_REVERSE_TYPE &tree, Int l, Int r, const F &f)
    requires(L)
{
    tree.apply(l, r, f);
}

// 両端を含む区間に作用を適用する。要素数nとして期待O(log(n+1))。
template <class S, class F, bool A, bool L, class I, class J>
void apply(CPLIB_REVERSE_TYPE &tree, ClosedSlice<I, J> range, const F &f)
    requires(L)
{
    tree.apply(range, f);
}

// 単調な述語を満たす最大の右端を返す。要素数nとして期待O(log(n+1))回の演算と述語呼出し。
template <class S, class F, bool A, bool L, class Pred>
Int max_right(CPLIB_REVERSE_TYPE &tree, Int l, Pred predicate)
    requires(A)
{
    return tree.max_right(l, predicate);
}

// 単調な述語を満たす最小の左端を返す。要素数nとして期待O(log(n+1))回の演算と述語呼出し。
template <class S, class F, bool A, bool L, class Pred>
Int min_left(CPLIB_REVERSE_TYPE &tree, Int r, Pred predicate)
    requires(A)
{
    return tree.min_left(r, predicate);
}

#undef CPLIB_REVERSE_TEMPLATE
#undef CPLIB_REVERSE_TYPE
}
