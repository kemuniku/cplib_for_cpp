#pragma once
#include <cplib/common.hpp>
#include <memory>
#include <string_view>
#include <stdexcept>
#include <ranges>
#include <type_traits>

namespace cplib {
struct PalindromicTreeNode;
class PalindromicTree;
class PalindromicTreeLinks {
    Int alphabet_ = 0, first_ = 0, count_ = 0, single_key_ = 0;
    PalindromicTreeNode *single_child_ = nullptr;
    std::vector<PalindromicTreeNode *> branches_;
    friend class PalindromicTree;
    void set(Int value, PalindromicTreeNode *child) {
        assert(__int128(value) >= first_ && __int128(value) - first_ < alphabet_ && !hasKey(value));
        if (count_ == 0) {
            single_key_ = value;
            single_child_ = child;
        } else {
            if (count_ == 1) {
                branches_.resize(alphabet_);
                branches_[single_key_ - first_] = single_child_;
                single_child_ = nullptr;
            }
            branches_[value - first_] = child;
        }
        ++count_;
    }
public:
    // 指定範囲の空の遷移ビューを作る。計算量 O(1)。
    explicit PalindromicTreeLinks(Int alphabet = 0, Int first = 0) : alphabet_(alphabet), first_(first) {
    }
    // 登録済みの子の個数を返す。計算量 O(1)。
    Int size() const {
        return count_;
    }
    // 対応する子を返す。範囲外・不在はnullptr。計算量 O(1)。
    PalindromicTreeNode *getOrDefault(Int value) const {
        if (__int128(value) < first_ || __int128(value) - first_ >= alphabet_ || count_ == 0)
            return nullptr;
        if (count_ == 1)
            return value == single_key_ ? single_child_ : nullptr;
        return branches_[value - first_];
    }
    // 対応する子が存在するか返す。計算量 O(1)。
    bool hasKey(Int value) const {
        return getOrDefault(value) != nullptr;
    }
    // 対応する子を返す。不在はout_of_range。計算量 O(1)。
    PalindromicTreeNode *operator[](Int value) const {
        auto child = getOrDefault(value);
        if (!child)
            throw std::out_of_range("回文木の遷移が存在しません");
        return child;
    }
    // 要素と子の組を返す。時間 O(σ+1)、出力領域 O(子の数)。
    std::vector<std::pair<Int, PalindromicTreeNode *>> pairs() const {
        if (count_ == 0)
            return {};
        if (count_ == 1)
            return {{single_key_, single_child_}};
        std::vector<std::pair<Int, PalindromicTreeNode *>> result;
        result.reserve(count_);
        for (Int value = 0; value < alphabet_; ++value)
            if (branches_[value])
                result.emplace_back(first_ + value, branches_[value]);
        return result;
    }
};
struct PalindromicTreeNode {
    PalindromicTreeLinks link;
    PalindromicTreeNode *parent = nullptr, *suffix_link = nullptr, *series_link = nullptr;
    Int length = 0, number = 0, character = 0;
    Int terminal_count = 0, total_count = 0;
    std::shared_ptr<bool> updated;

    // 回文の長さを返す。計算量 O(1)。
    Int size() const {
        return length;
    }
    // 既存コード用に回文の長さを返す。計算量 O(1)。
    Int len() const {
        return size();
    }
    // 再集計後は総出現数、それ以外は最長回文接尾辞としての出現数を返す。計算量 O(1)。
    Int count() const {
        return *updated ? total_count : terminal_count;
    }
    // 節点番号を返す。長さ-1と0の根はそれぞれ0と1。計算量 O(1)。
    Int id() const {
        return number;
    }
    // suffix link先との長さの差を返す。根は0。計算量 O(1)。
    Int diff() const {
        return number < 2 ? 0 : length - suffix_link->length;
    }
    // 親リンクをたどって回文を復元する。時間・出力領域 O(回文長+1)。
    std::vector<Int> get_palindrome() const {
        std::vector<Int> result(std::max<Int>(0, length));
        auto node = this;
        Int left = 0, right = Int(result.size()) - 1;
        while (node->length > 0) {
            result[left++] = node->character;
            result[right--] = node->character;
            node = node->parent;
        }
        return result;
    }
};
class PalindromicTree {
    Int alphabet_, first_;
    std::vector<Int> data_;
    std::shared_ptr<bool> updated_ = std::make_shared<bool>(false);
    PalindromicTreeNode *newNode(Int length) {
        auto node = std::make_unique<PalindromicTreeNode>();
        node->length = length;
        node->number = nodes.size();
        node->updated = updated_;
        node->link = PalindromicTreeLinks(alphabet_, first_);
        nodes.push_back(std::move(node));
        return nodes.back().get();
    }
    PalindromicTreeNode *findLongest(PalindromicTreeNode *node) const {
        Int pos = Int(data_.size()) - 1;
        while (true) {
            Int left = pos - node->length - 1;
            if (left >= 0 && data_[left] == data_[pos])
                return node;
            node = node->suffix_link;
        }
    }
public:
    std::vector<std::unique_ptr<PalindromicTreeNode>> nodes;
    PalindromicTreeNode *last_node = nullptr;
    std::vector<PalindromicTreeNode *> prefix_nodes;

    // 空の回文木を作る。文字はfirst..<first+alphabet、事前確保はcapacity個分。時間・領域 O(capacity+σ)。
    explicit PalindromicTree(Int alphabet = 26, Int capacity = 0, Int first = 0)
        : alphabet_(alphabet), first_(first) {
        assert(alphabet >= 0 && capacity >= 0);
        assert(__int128(first) + std::max<Int>(0, alphabet - 1) <= std::numeric_limits<Int>::max());
        data_.reserve(capacity);
        nodes.reserve(capacity + 2);
        prefix_nodes.reserve(capacity);
        newNode(-1);
        newNode(0);
        nodes[1]->suffix_link = nodes[1]->series_link = nodes[0].get();
        last_node = nodes[1].get();
    }
    // 整数列から回文木を作る。amax=-1なら文字範囲を推定する。分岐節点数Bとして時間・領域 O(N+σB+1)。
    explicit PalindromicTree(const std::vector<Int> &a, Int amax = -1)
        : PalindromicTree(amax < 0 ? (a.empty() ? 0 : *std::max_element(a.begin(), a.end()) + 1) : amax, a.size()) {
        for (Int value : a)
            add(value);
    }
    // 独立した木を複製し、すべてのリンクを張り直す。時間・領域 O((N+1)(σ+1))。
    PalindromicTree(const PalindromicTree &other)
        : alphabet_(other.alphabet_), first_(other.first_), data_(other.data_), updated_(std::make_shared<bool>(*other.updated_)) {
        for (const auto &node : other.nodes) {
            nodes.push_back(std::make_unique<PalindromicTreeNode>(*node));
            nodes.back()->updated = updated_;
        }
        auto remap = [&](PalindromicTreeNode *node) {
            return node ? nodes[node->number].get() : nullptr;
        };
        for (auto &node : nodes) {
            node->link.single_child_ = remap(node->link.single_child_);
            for (auto &child : node->link.branches_)
                child = remap(child);
            node->parent = remap(node->parent);
            node->suffix_link = remap(node->suffix_link);
            node->series_link = remap(node->series_link);
        }
        last_node = remap(other.last_node);
        for (auto node : other.prefix_nodes)
            prefix_nodes.push_back(remap(node));
    }
    // 独立した木を複製して代入する。新旧要素数の和をNとして時間 O((N+1)(σ+1))。
    PalindromicTree &operator=(const PalindromicTree &other) {
        if (this != &other) {
            PalindromicTree copy(other);
            *this = std::move(copy);
        }
        return *this;
    }
    // 木の所有権を移す。計算量 O(1)。
    PalindromicTree(PalindromicTree &&) = default;
    // 木の所有権を移す。旧木の破棄に O((N+1)(σ+1))。
    PalindromicTree &operator=(PalindromicTree &&) = default;
    // 現在の列の長さを返す。計算量 O(1)。
    Int size() const {
        return data_.size();
    }
    // 一要素追加し最長回文接尾辞を返す。償却 O(σ+1)、一回の最悪 O(N+σ)。
    PalindromicTreeNode *add(Int value) {
        assert(__int128(value) >= first_ && __int128(value) - first_ < alphabet_);
        data_.push_back(value);
        auto parent = findLongest(last_node);
        auto node = parent->link.getOrDefault(value);
        if (!node) {
            node = newNode(parent->length + 2);
            node->parent = parent;
            node->character = value;
            node->suffix_link = parent == nodes[0].get()
                                    ? nodes[1].get()
                                    : findLongest(parent->suffix_link)->link[value];
            auto suffix = node->suffix_link;
            node->series_link = node->diff() == suffix->diff() ? suffix->series_link : suffix;
            parent->link.set(value, node);
        }
        ++node->terminal_count;
        prefix_nodes.push_back(node);
        last_node = node;
        *updated_ = false;
        return node;
    }
    // 指定節点の回文を復元する。時間・出力領域 O(回文長+1)。
    auto get_palindrome(const PalindromicTreeNode *node) const {
        return node->get_palindrome();
    }
    // 所有節点の回文を復元する。時間・出力領域 O(回文長+1)。
    auto get_palindrome(const std::unique_ptr<PalindromicTreeNode> &node) const {
        return node->get_palindrome();
    }
    // 総出現回数を再集計する。繰り返し呼び出せる。時間 O(節点数)、追加領域 O(1)。
    void update_count() {
        for (auto &node : nodes)
            node->total_count = node->terminal_count;
        for (std::size_t i = nodes.size(); i-- > 1;)
            nodes[i]->suffix_link->total_count += nodes[i]->total_count;
        *updated_ = true;
    }
};
// 整数列から回文木を作る。時間・領域 O((N+1)(σ+1))。
inline auto initPalindromicTree(const std::vector<Int> &a, Int amax = -1) {
    return PalindromicTree(a, amax);
}
// 指定文字種数の空の回文木を作る。時間・領域 O(capacity+σ)。
inline auto initPalindromicTree(Int alphabet = 26, Int capacity = 0) {
    return PalindromicTree(alphabet, capacity);
}
// 英小文字列を0..25へ変換して回文木を作る。時間・領域 O(N+1)。
inline auto initPalindromicTree(std::string_view s, char first = 'a') {
    PalindromicTree tree(26, s.size());
    for (unsigned char c : s)
        tree.add(Int(c) - static_cast<unsigned char>(first));
    return tree;
}
// 節点の回文を復元する。時間・出力領域 O(回文長+1)。
inline auto get_palindrome(const PalindromicTreeNode *node) {
    return node->get_palindrome();
}
// 指定節点の回文を復元する。時間・出力領域 O(回文長+1)。
inline auto get_palindrome(const PalindromicTree &tree, const PalindromicTreeNode *node) {
    return tree.get_palindrome(node);
}
// 総出現回数を再集計する。時間 O(節点数)。
inline void update_count(PalindromicTree &tree) {
    tree.update_count();
}
// 既存コード用に回文長を返す。計算量 O(1)。
inline Int len(const PalindromicTreeNode &node) {
    return node.size();
}
// 節点の出現数を返す。計算量 O(1)。
inline Int count(const PalindromicTreeNode &node) {
    return node.count();
}
// 節点番号を返す。計算量 O(1)。
inline Int id(const PalindromicTreeNode &node) {
    return node.id();
}
}

namespace cplib {
namespace detail {
template <class T> constexpr Int palindromicOrdinal(T value) {
    if constexpr (std::is_same_v<T, char>)
        return static_cast<unsigned char>(value);
    else
        return Int(value);
}
}
// コンパイル時の要素範囲を持つ回文木。文字・整数・列挙型を扱える。
template <auto Low, auto High> class ArrayPalindromicTree : public PalindromicTree {
    using T = decltype(Low);
    static_assert(std::is_same_v<T, decltype(High)>);
    static_assert(std::is_integral_v<T> || std::is_enum_v<T>);
    static constexpr Int first = detail::palindromicOrdinal(Low);
    static constexpr Int last = detail::palindromicOrdinal(High);
    static_assert(first <= last && __int128(last) - first < std::numeric_limits<Int>::max());
public:
    // 静的な要素範囲の空の木を作る。時間・領域 O(capacity+1)。
    explicit ArrayPalindromicTree(Int capacity = 0) : PalindromicTree(last - first + 1, capacity, first) {
    }
    // 静的な要素範囲で列から構築する。分岐節点数Bとして時間・領域 O(N+σB+1)。
    template <std::ranges::input_range Range> explicit ArrayPalindromicTree(const Range &values) : ArrayPalindromicTree() {
        for (T value : values)
            add(value);
    }
    // NUL終端文字列から構築する。時間・領域 O(N+σB+1)。
    explicit ArrayPalindromicTree(const char *values) requires std::is_same_v<T, char>
        : ArrayPalindromicTree(std::string_view(values)) {
    }
    // 末尾に要素を追加し最長回文接尾辞を返す。償却 O(σ+1)、一回の最悪 O(N+σ)。
    PalindromicTreeNode *add(T value) {
        return PalindromicTree::add(detail::palindromicOrdinal(value));
    }
    // 元の要素型で回文を復元する。時間・出力領域 O(回文長+1)。
    std::vector<T> get_palindrome(const PalindromicTreeNode *node) const {
        auto values = node->get_palindrome();
        std::vector<T> result;
        result.reserve(values.size());
        for (Int value : values)
            result.push_back(static_cast<T>(value));
        return result;
    }
    // 所有節点から元の要素型で回文を復元する。時間・出力領域 O(回文長+1)。
    std::vector<T> get_palindrome(const std::unique_ptr<PalindromicTreeNode> &node) const {
        return get_palindrome(node.get());
    }
};
// 静的な要素範囲の空の回文木を作る。時間・領域 O(capacity+1)。
template <auto Low, auto High> auto initPalindromicTree(Int capacity = 0) {
    return ArrayPalindromicTree<Low, High>(capacity);
}
// 静的な要素範囲の列から回文木を作る。分岐節点数Bとして時間・領域 O(N+σB+1)。
template <auto Low, auto High, std::ranges::input_range Range> auto initPalindromicTree(const Range &values) {
    return ArrayPalindromicTree<Low, High>(values);
}
// 文字範囲と事前確保長から空の回文木を作る。時間・領域 O(capacity+1)。
inline auto initPalindromicTree(char first, char last, Int capacity = 0) {
    Int low = static_cast<unsigned char>(first), high = static_cast<unsigned char>(last);
    assert(low <= high);
    return PalindromicTree(high - low + 1, capacity, low);
}
}
