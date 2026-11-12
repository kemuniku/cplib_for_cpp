#pragma once
#include <cplib/common.hpp>
#include <unordered_map>
#include <string_view>

namespace cplib {
template <class T> struct SuffixAutomatonTableNode {
    Int length = 0, link = -1;
    std::unordered_map<T, Int> next;
    // この状態が表す文字列の最大長を返す。計算量 O(1)。
    Int size() const {
        return length;
    }
};
template <class T = char> class SuffixAutomatonTable {
public:
    std::vector<SuffixAutomatonTableNode<T>> nodes;
    Int last = 0;
    // 根だけを持つ空のSAMを作る。時間・領域 O(1)。
    SuffixAutomatonTable() : nodes(1) {
    }
    // 列からSAMを作る。期待時間・領域 O(N+1)。
    template <class Range> explicit SuffixAutomatonTable(const Range &s) : SuffixAutomatonTable() {
        for (T c : s)
            extend(c);
    }
    // NUL終端文字列から構築する。期待時間・領域 O(N+1)。
    explicit SuffixAutomatonTable(const char *s) requires std::is_same_v<T, char> : SuffixAutomatonTable(std::string_view(s)) {
    }
    // 一要素追加し全体の状態IDを返す。N回で期待 O(N)、一回の最悪 O(N)。
    Int extend(T c) {
        Int current = nodes.size();
        nodes.push_back({nodes[last].length + 1, 0, {}});
        Int p = last;
        while (p >= 0 && !nodes[p].next.contains(c)) {
            nodes[p].next[c] = current;
            p = nodes[p].link;
        }
        if (p >= 0) {
            Int q = nodes[p].next.at(c);
            if (nodes[p].length + 1 == nodes[q].length) {
                nodes[current].link = q;
            } else {
                Int clone = nodes.size();
                auto copy = nodes[q];
                copy.length = nodes[p].length + 1;
                nodes.push_back(std::move(copy));
                while (p >= 0 && next(p, c) == q) {
                    nodes[p].next[c] = clone;
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
    // 状態からの遷移先を返す。不在は-1、表は変更しない。期待 O(1)。
    Int next(Int node, T c) const {
        auto it = nodes[node].next.find(c);
        return it == nodes[node].next.end() ? -1 : it->second;
    }
    // 連続部分文字列の状態を返す。空は根、不在は-1。期待 O(M+1)。
    template <class Range> Int findNode(const Range &s) const {
        Int node = 0;
        for (T c : s) {
            node = next(node, c);
            if (node < 0)
                break;
        }
        return node;
    }
    // 連続部分文字列が存在するか返す。空も含む。期待 O(M+1)。
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
// 空のSAMを作る。時間・領域 O(1)。
template <class T = char> auto initSuffixAutomatonTable() {
    return SuffixAutomatonTable<T>();
}
// 列からSAMを作る。期待時間・領域 O(N+1)。
template <class Range> auto initSuffixAutomatonTable(const Range &s) {
    using T = std::remove_cvref_t<decltype(*std::begin(s))>;
    return SuffixAutomatonTable<T>(s);
}
// 文字列からSAMを作る。期待時間・領域 O(N+1)。
inline auto initSuffixAutomatonTable(const char *s) {
    return SuffixAutomatonTable<char>(std::string_view(s));
}
}
