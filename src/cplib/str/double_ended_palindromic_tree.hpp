#pragma once
#include <cplib/common.hpp>
#include <deque>
#include <stdexcept>
#include <string_view>

namespace cplib {
class DoubleEndedPalindromicTree {
    struct Node {
        Int length = 0, parent = 0, suffix = 0, depth = 0, longestCount = 0, suffixChildren = 0;
    };

    struct Entry {
        Int value, prefixSurface = 1, suffixSurface = 1;
    };

    Int sigma, freeHead = 0, freeCount = 0, total = 0;
    char offset;
    std::vector<Node> nodes;
    std::vector<std::int32_t> children, direct;
    std::deque<Entry> data;

    Int addNode(Int parent, Int suffix, Int value, Int preceding) {
        Int out;
        if (freeCount) {
            out = freeHead;
            freeHead = nodes[out].parent;
            --freeCount;
        } else {
            out = nodes.size();
            assert(out <= std::numeric_limits<std::int32_t>::max());
            nodes.emplace_back();
            children.resize(nodes.size() * sigma);
            direct.resize(nodes.size() * sigma);
        }
        nodes[out] = {nodes[parent].length + 2, parent, suffix, nodes[suffix].depth + 1, 0, 0};
        std::fill_n(children.begin() + out * sigma, sigma, 0);
        std::copy_n(direct.begin() + suffix * sigma, sigma, direct.begin() + out * sigma);
        direct[out * sigma + preceding] = suffix;
        children[parent * sigma + value] = out;
        ++nodes[suffix].suffixChildren;
        return out;
    }

    template <bool Front> Entry &entry(Int i) {
        return data[Front ? i : len() - 1 - i];
    }

    template <bool Front> Int &near(Int i) {
        if constexpr (Front)
            return entry<Front>(i).prefixSurface;
        else
            return entry<Front>(i).suffixSurface;
    }

    template <bool Front> Int &far(Int i) {
        if constexpr (Front)
            return entry<Front>(i).suffixSurface;
        else
            return entry<Front>(i).prefixSurface;
    }

    template <bool Front> void push(Int value) {
        assert(0 <= value && value < sigma);
        Int parent = data.empty() ? 1 : near<Front>(0);
        if constexpr (Front)
            data.push_front({value});
        else
            data.push_back({value});
        Int opposite = nodes[parent].length + 1;
        if (opposite >= len() || entry<Front>(opposite).value != value)
            parent = direct[parent * sigma + value];
        Int node = children[parent * sigma + value];
        if (!node) {
            Int suffix = parent == 0 ? 1 : children[direct[parent * sigma + value] * sigma + value];
            node = addNode(parent, suffix, value, entry<Front>(nodes[suffix].length).value);
        }
        Int length = nodes[node].length, suffix = nodes[node].suffix, sl = nodes[suffix].length;
        near<Front>(0) = node;
        far<Front>(length - 1) = node;
        if (sl > 0 && near<Front>(length - sl) == suffix)
            near<Front>(length - sl) = 1;
        ++nodes[node].longestCount;
        total += nodes[node].depth;
    }

    template <bool Front> void pop() {
        if (data.empty())
            throw std::out_of_range("empty palindromic tree");
        Int node = near<Front>(0), suffix = nodes[node].suffix, length = nodes[node].length,
            sl = nodes[suffix].length;
        far<Front>(length - 1) = 1;
        if (sl > 0 && nodes[near<Front>(length - sl)].length < sl) {
            near<Front>(length - sl) = suffix;
            far<Front>(length - 1) = suffix;
        }
        --nodes[node].longestCount;
        total -= nodes[node].depth;
        if (!nodes[node].longestCount && !nodes[node].suffixChildren) {
            children[nodes[node].parent * sigma + entry<Front>(0).value] = 0;
            --nodes[suffix].suffixChildren;
            nodes[node].parent = freeHead;
            freeHead = node;
            ++freeCount;
        }
        if constexpr (Front)
            data.pop_front();
        else
            data.pop_back();
    }

public:
    // 空の回文木を作り、capacity文字分を事前確保する。O(amax * (capacity + 1))。
    // 整数は 0..<amax、文字は ord(ch)-ord(c) として扱う。capacityを超えても追加可能。
    explicit DoubleEndedPalindromicTree(Int amax = 26, char c = 'a', Int capacity = 0)
        : sigma(amax), offset(c) {
        assert(amax > 0 && capacity >= 0 &&
               capacity <= std::numeric_limits<std::int32_t>::max() - 2);
        nodes.reserve(capacity + 2);
        nodes.resize(2);
        nodes[0].length = -1;
        children.reserve((capacity + 2) * sigma);
        direct.reserve((capacity + 2) * sigma);
        children.resize(2 * sigma);
        direct.resize(2 * sigma);
    }

    // 現在の文字列の長さを O(1) で返す。
    Int len() const {
        return data.size();
    }

    Int size() const {
        return len();
    }

    // 現在の文字列が空かを O(1) で返す。
    bool isEmpty() const {
        return data.empty();
    }

    // 現在の文字列に含まれる異なる非空回文の個数を O(1) で返す。
    Int count_distinct_palindromes() const {
        return nodes.size() - freeCount - 2;
    }

    // 現在の文字列の非空回文の出現総数を O(1) で返す。位置が異なる出現も数える。
    Int count_palindromes() const {
        return total;
    }

    // 最長回文接頭辞の長さを O(1) で返す。空文字列なら0。
    Int longest_prefix_palindrome() const {
        return data.empty() ? 0 : nodes[data.front().prefixSurface].length;
    }

    // 最長回文接尾辞の長さを O(1) で返す。空文字列なら0。
    Int longest_suffix_palindrome() const {
        return data.empty() ? 0 : nodes[data.back().suffixSurface].length;
    }

    // 先頭に0 <= v < amaxの整数を追加する。償却O(σ)。σは文字種数。
    void push_front(Int v) {
        push<true>(v);
    }

    // 末尾に0 <= v < amaxの整数を追加する。償却O(σ)。σは文字種数。
    void push_back(Int v) {
        push<false>(v);
    }

    // 先頭に文字cを追加する。初期化時の文字との差が[0,amax)に入る必要がある。償却O(σ)。
    void push_front(char c) {
        push_front(Int(static_cast<unsigned char>(c)) - static_cast<unsigned char>(offset));
    }

    // 末尾に文字cを追加する。初期化時の文字との差が[0,amax)に入る必要がある。償却O(σ)。
    void push_back(char c) {
        push_back(Int(static_cast<unsigned char>(c)) - static_cast<unsigned char>(offset));
    }

    // 先頭を一文字削除する。O(1)。空ならstd::out_of_range。
    void pop_front() {
        pop<true>();
    }

    // 末尾を一文字削除する。O(1)。空ならstd::out_of_range。
    void pop_back() {
        pop<false>();
    }
};

inline auto initDoubleEndedPalindromicTree(Int amax = 26, char c = 'a', Int capacity = 0) {
    return DoubleEndedPalindromicTree(amax, c, capacity);
}

// 非負整数列から O(σ(|a| + 1)) で構築する。amax省略時は最大値+1（空なら1）。
inline auto initDoubleEndedPalindromicTree(std::span<const Int> a, Int amax = -1) {
    if (amax < 0) {
        amax = 1;
        for (Int x : a) {
            assert(0 <= x && x < std::numeric_limits<Int>::max());
            amax = std::max(amax, x + 1);
        }
    }
    DoubleEndedPalindromicTree t(amax, 'a', a.size());
    for (Int x : a)
        t.push_back(x);
    return t;
}

// 文字列から O(amax(|s| + 1)) で構築する。既定の文字範囲は 'a'..'z'。
inline auto initDoubleEndedPalindromicTree(std::string_view s, char c = 'a', Int amax = 26) {
    DoubleEndedPalindromicTree t(amax, c, s.size());
    for (char x : s)
        t.push_back(x);
    return t;
}
}
