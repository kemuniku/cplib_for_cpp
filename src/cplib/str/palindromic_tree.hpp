#pragma once
#include <cplib/common.hpp>
#include <memory>
#include <string_view>

namespace cplib {
struct PalindromicTreeNode {
    std::vector<PalindromicTreeNode *> link;
    PalindromicTreeNode *suffix_link = nullptr;
    Int length = 0, occurrences = 0, number = 0;

    Int len() const {
        return length;
    }

    Int count() const {
        return occurrences;
    }

    Int id() const {
        return number;
    }
};

class PalindromicTree {
    Int amax_ = 0;

    PalindromicTreeNode *new_node(Int length) {
        auto n = std::make_unique<PalindromicTreeNode>();
        n->link.resize(amax_);
        n->length = length;
        n->number = nodes.size();
        nodes.push_back(std::move(n));
        return nodes.back().get();
    }

public:
    std::vector<std::unique_ptr<PalindromicTreeNode>> nodes;
    PalindromicTreeNode *last_node = nullptr;

    explicit PalindromicTree(const std::vector<Int> &a, Int amax = -1)
        : amax_(amax < 0 ? (a.empty() ? 0 : *std::max_element(a.begin(), a.end()) + 1) : amax) {
        new_node(-1);
        new_node(0);
        nodes[1]->suffix_link = nodes[0].get();
        auto find_longest = [&](Int pos, PalindromicTreeNode *n) {
            while (pos - n->length - 1 < 0 || a[pos - n->length - 1] != a[pos])
                n = n->suffix_link;
            return n;
        };
        auto current = nodes[0].get();
        for (Int i = 0; i < Int(a.size()); ++i) {
            assert(0 <= a[i] && a[i] < amax_);
            current = find_longest(i, current);
            if (!current->link[a[i]])
                current->link[a[i]] = new_node(current->length + 2);
            if (current == nodes[0].get())
                current->link[a[i]]->suffix_link = nodes[1].get();
            else
                current->link[a[i]]->suffix_link =
                    find_longest(i, current->suffix_link)->link[a[i]];
            current = current->link[a[i]];
            ++current->occurrences;
        }
        last_node = current;
    }

    PalindromicTree(const PalindromicTree &other) : amax_(other.amax_) {
        for (const auto &n : other.nodes)
            nodes.push_back(std::make_unique<PalindromicTreeNode>(*n));
        for (auto &n : nodes) {
            for (auto &p : n->link)
                if (p)
                    p = nodes[p->number].get();
            if (n->suffix_link)
                n->suffix_link = nodes[n->suffix_link->number].get();
        }
        last_node = other.last_node ? nodes[other.last_node->number].get() : nullptr;
    }

    PalindromicTree &operator=(const PalindromicTree &other) {
        if (this != &other) {
            PalindromicTree copy(other);
            *this = std::move(copy);
        }
        return *this;
    }

    PalindromicTree(PalindromicTree &&) = default;
    PalindromicTree &operator=(PalindromicTree &&) = default;

    // 元と同じリンク木探索で回文を復元。O(ノード数×文字種類数)。
    std::vector<Int> get_palindrome(const PalindromicTreeNode *target) const {
        if (target == nodes[0].get() || target == nodes[1].get())
            return {};
        std::vector<Int> ans;
        auto dfs = [&](auto &&self, const PalindromicTreeNode *x) -> bool {
            if (x == target)
                return true;
            for (Int i = 0; i < amax_; ++i)
                if (x->link[i] && self(self, x->link[i])) {
                    ans.push_back(i);
                    return true;
                }
            return false;
        };
        dfs(dfs, nodes[0].get());
        dfs(dfs, nodes[1].get());
        for (Int i = ans.size(); i < target->length; ++i)
            ans.push_back(ans[target->length - 1 - i]);
        return ans;
    }

    auto get_palindrome(const std::unique_ptr<PalindromicTreeNode> &n) const {
        return get_palindrome(n.get());
    }

    // 最長接尾回文の出現数を suffix link に加算。元どおり呼ぶたびに加算する。
    void update_count() {
        for (std::size_t i = nodes.size(); i-- > 1;)
            nodes[i]->suffix_link->occurrences += nodes[i]->occurrences;
    }
};

inline auto initPalindromicTree(const std::vector<Int> &a, Int amax = -1) {
    return PalindromicTree(a, amax);
}

inline auto initPalindromicTree(std::string_view s, char c = 'a') {
    std::vector<Int> a;
    for (unsigned char x : s)
        a.push_back(Int(x) - static_cast<unsigned char>(c));
    return PalindromicTree(a, 26);
}

inline auto get_palindrome(const PalindromicTree &t, const PalindromicTreeNode *n) {
    return t.get_palindrome(n);
}

inline void update_count(PalindromicTree &t) {
    t.update_count();
}

inline Int len(const PalindromicTreeNode &n) {
    return n.len();
}

inline Int count(const PalindromicTreeNode &n) {
    return n.count();
}

inline Int id(const PalindromicTreeNode &n) {
    return n.id();
}
}
