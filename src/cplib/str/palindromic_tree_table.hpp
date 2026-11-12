#pragma once
#include <cplib/common.hpp>
#include <memory>
#include <string_view>
#include <unordered_map>

namespace cplib {
template <class T> struct TablePalindromicTreeNode {
    std::unordered_map<T, TablePalindromicTreeNode<T> *> link;
    TablePalindromicTreeNode<T> *parent = nullptr, *suffix_link = nullptr, *series_link = nullptr;
    Int length = 0, number = 0;
    T character{};
    Int terminal_count = 0, total_count = 0;
    std::shared_ptr<bool> updated;

    // 回文の長さを返す。計算量 O(1)。
    Int size() const {
        return length;
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
    std::vector<T> get_palindrome() const {
        std::vector<T> result(std::max<Int>(0, length));
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
template <class T = char> class TablePalindromicTree {
    std::vector<T> data_;
    std::shared_ptr<bool> updated_ = std::make_shared<bool>(false);
    TablePalindromicTreeNode<T> *newNode(Int length) {
        auto node = std::make_unique<TablePalindromicTreeNode<T>>();
        node->length = length;
        node->number = nodes.size();
        node->updated = updated_;
        nodes.push_back(std::move(node));
        return nodes.back().get();
    }
    TablePalindromicTreeNode<T> *findLongest(TablePalindromicTreeNode<T> *node) const {
        Int pos = Int(data_.size()) - 1;
        while (true) {
            Int left = pos - node->length - 1;
            if (left >= 0 && data_[left] == data_[pos])
                return node;
            node = node->suffix_link;
        }
    }
public:
    std::vector<std::unique_ptr<TablePalindromicTreeNode<T>>> nodes;
    TablePalindromicTreeNode<T> *last_node = nullptr;
    std::vector<TablePalindromicTreeNode<T> *> prefix_nodes;

    // 根だけを持つ空の回文木を作る。時間・領域 O(1)。
    TablePalindromicTree() {
        newNode(-1);
        newNode(0);
        nodes[1]->suffix_link = nodes[1]->series_link = nodes[0].get();
        last_node = nodes[1].get();
    }
    // 列から回文木を作る。期待時間・領域 O(N+1)。
    template <class Range> explicit TablePalindromicTree(const Range &a) : TablePalindromicTree() {
        for (const T &value : a)
            add(value);
    }
    // 独立した木を複製し、すべてのリンクを張り直す。時間・領域 O(N+1)。
    TablePalindromicTree(const TablePalindromicTree &other)
        : data_(other.data_), updated_(std::make_shared<bool>(*other.updated_)) {
        for (const auto &node : other.nodes) {
            nodes.push_back(std::make_unique<TablePalindromicTreeNode<T>>(*node));
            nodes.back()->updated = updated_;
        }
        auto remap = [&](TablePalindromicTreeNode<T> *node) {
            return node ? nodes[node->number].get() : nullptr;
        };
        for (auto &node : nodes) {
            for (auto &[value, child] : node->link)
                child = remap(child);
            node->parent = remap(node->parent);
            node->suffix_link = remap(node->suffix_link);
            node->series_link = remap(node->series_link);
        }
        last_node = remap(other.last_node);
        for (auto node : other.prefix_nodes)
            prefix_nodes.push_back(remap(node));
    }
    // 独立した木を複製して代入する。新旧要素数の和をNとして時間 O(N+1)。
    TablePalindromicTree &operator=(const TablePalindromicTree &other) {
        if (this != &other) {
            TablePalindromicTree copy(other);
            *this = std::move(copy);
        }
        return *this;
    }
    // 木の所有権を移す。計算量 O(1)。
    TablePalindromicTree(TablePalindromicTree &&) = default;
    // 木の所有権を移す。旧木の破棄に O(N+1)。
    TablePalindromicTree &operator=(TablePalindromicTree &&) = default;
    // 現在の列の長さを返す。計算量 O(1)。
    Int size() const {
        return data_.size();
    }
    // NUL終端文字列から構築する。期待時間・領域 O(N+1)。
    explicit TablePalindromicTree(const char *s) requires std::is_same_v<T, char>
        : TablePalindromicTree(std::string_view(s)) {
    }
    // 一要素追加し最長回文接尾辞を返す。期待償却 O(1)、一回の最悪 O(N)。
    TablePalindromicTreeNode<T> *add(T value) {
        data_.push_back(value);
        auto parent = findLongest(last_node);
        auto node = parent->link[value];
        if (!node) {
            node = newNode(parent->length + 2);
            node->parent = parent;
            node->character = value;
            node->suffix_link = parent == nodes[0].get()
                                    ? nodes[1].get()
                                    : findLongest(parent->suffix_link)->link[value];
            auto suffix = node->suffix_link;
            node->series_link = node->diff() == suffix->diff() ? suffix->series_link : suffix;
            parent->link[value] = node;
        }
        ++node->terminal_count;
        prefix_nodes.push_back(node);
        last_node = node;
        *updated_ = false;
        return node;
    }
    // 指定節点の回文を復元する。時間・出力領域 O(回文長+1)。
    auto get_palindrome(const TablePalindromicTreeNode<T> *node) const {
        return node->get_palindrome();
    }
    // 所有節点の回文を復元する。時間・出力領域 O(回文長+1)。
    auto get_palindrome(const std::unique_ptr<TablePalindromicTreeNode<T>> &node) const {
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
// 空のTable版回文木を作る。時間・領域 O(1)。
template <class T = char> auto initTablePalindromicTree() {
    return TablePalindromicTree<T>();
}
// 列からTable版回文木を作る。期待時間・領域 O(N+1)。
template <class Range> auto initTablePalindromicTree(const Range &a) {
    using T = std::remove_cvref_t<decltype(*std::begin(a))>;
    return TablePalindromicTree<T>(a);
}
// 文字列からTable版回文木を作る。期待時間・領域 O(N+1)。
inline auto initTablePalindromicTree(const char *s) {
    return TablePalindromicTree<char>(std::string_view(s));
}
// 節点の回文を復元する。時間・出力領域 O(回文長+1)。
template <class T> auto get_palindrome(const TablePalindromicTreeNode<T> *node) {
    return node->get_palindrome();
}
// 指定節点の回文を復元する。時間・出力領域 O(回文長+1)。
template <class T> auto get_palindrome(const TablePalindromicTree<T> &tree, const TablePalindromicTreeNode<T> *node) {
    return tree.get_palindrome(node);
}
// 総出現回数を再集計する。時間 O(節点数)。
template <class T> void update_count(TablePalindromicTree<T> &tree) {
    tree.update_count();
}
}
