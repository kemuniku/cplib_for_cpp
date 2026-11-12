#pragma once
#include <cplib/common.hpp>
#include <array>
#include <stdexcept>
#include <string_view>

namespace cplib {
template <unsigned char Low, unsigned char High> struct SuffixAutomatonArrayNode {
    static_assert(Low <= High);
    Int length = 0, link = -1;
    std::array<Int, High - Low + 1> next;
    // 全遷移を不在で初期化する。計算量 O(σ)。
    SuffixAutomatonArrayNode(Int length = 0, Int link = -1) : length(length), link(link) {
        next.fill(-1);
    }
    // この状態が表す文字列の最大長を返す。計算量 O(1)。
    Int size() const {
        return length;
    }
};
template <unsigned char Low = 'a', unsigned char High = 'z'> class SuffixAutomatonArray {
public:
    std::vector<SuffixAutomatonArrayNode<Low, High>> nodes;
    Int last = 0;
    // 根だけを持つ空のSAMを作る。時間・領域 O(σ)。
    SuffixAutomatonArray() : nodes(1) {
    }
    // 列からSAMを作る。時間・領域 O((N+1)σ)。
    template <class Range> explicit SuffixAutomatonArray(const Range &s) : SuffixAutomatonArray() {
        for (unsigned char c : s)
            extend(c);
    }
    // NUL終端文字列から構築する。時間・領域 O((N+1)σ)。
    explicit SuffixAutomatonArray(const char *s) : SuffixAutomatonArray(std::string_view(s)) {
    }
    // 一要素追加し全体の状態IDを返す。N回で O((N+1)σ)、一回の最悪 O((N+1)σ)。
    Int extend(unsigned char c) {
        if (c < Low || c > High)
            throw std::invalid_argument("文字がアルファベット範囲外です");
        Int key = c - Low;
        Int current = nodes.size();
        nodes.emplace_back(nodes[last].length + 1, 0);
        Int p = last;
        while (p >= 0 && nodes[p].next[key] < 0) {
            nodes[p].next[key] = current;
            p = nodes[p].link;
        }
        if (p >= 0) {
            Int q = nodes[p].next[key];
            if (nodes[p].length + 1 == nodes[q].length) {
                nodes[current].link = q;
            } else {
                Int clone = nodes.size();
                auto copy = nodes[q];
                copy.length = nodes[p].length + 1;
                nodes.push_back(std::move(copy));
                while (p >= 0 && next(p, c) == q) {
                    nodes[p].next[key] = clone;
                    p = nodes[p].link;
                }
                nodes[q].link = nodes[current].link = clone;
            }
        }
        last = current;
        return current;
    }
    // 根の状態IDを返す。計算量 O(1)。
    Int root() const {
        return 0;
    }
    // 根とcloneを含む状態数を返す。計算量 O(1)。
    Int nodeCount() const {
        return nodes.size();
    }
    // 入力列の長さを返す。計算量 O(1)。
    Int size() const {
        return nodes[last].length;
    }
    // 状態からの遷移先を返す。不在は-1、表は変更しない。O(1)。
    Int next(Int node, unsigned char c) const {
        return c < Low || c > High ? -1 : nodes[node].next[c - Low];
    }
    // 連続部分文字列の状態を返す。空は根、不在は-1。O(M+1)。
    template <class Range> Int findNode(const Range &s) const {
        Int node = 0;
        for (unsigned char c : s) {
            node = next(node, c);
            if (node < 0)
                break;
        }
        return node;
    }
    // 連続部分文字列が存在するか返す。空も含む。O(M+1)。
    template <class Range> bool contains(const Range &s) const {
        return findNode(s) >= 0;
    }
    // NUL終端文字列の状態を返す。時間 O(M+1)（Table版は期待時間）。
    Int findNode(const char *s) const {
        return findNode(std::string_view(s));
    }
    // NUL終端文字列が含まれるか返す。時間 O(M+1)（Table版は期待時間）。
    bool contains(const char *s) const {
        return findNode(std::string_view(s)) >= 0;
    }
    // 空を除く異なる連続部分文字列の個数を返す。時間 O(状態数)、追加領域 O(1)。
    Int countDistinctSubstrings() const {
        Int result = 0;
        for (Int v = 1; v < nodeCount(); ++v)
            result += nodes[v].length - nodes[nodes[v].link].length;
        return result;
    }
};
// 指定文字範囲の空のSAMを作る。時間・領域 O(σ)。
template <unsigned char Low = 'a', unsigned char High = 'z'> auto initSuffixAutomatonArray() {
    return SuffixAutomatonArray<Low, High>();
}
// 指定文字範囲の列からSAMを作る。時間・領域 O((N+1)σ)。
template <unsigned char Low = 'a', unsigned char High = 'z', class Range>
auto initSuffixAutomatonArray(const Range &s) {
    return SuffixAutomatonArray<Low, High>(s);
}
// 文字列からSAMを作る。時間・領域 O((N+1)σ)。
template <unsigned char Low = 'a', unsigned char High = 'z'>
auto initSuffixAutomatonArray(const char *s) {
    return SuffixAutomatonArray<Low, High>(std::string_view(s));
}
}
