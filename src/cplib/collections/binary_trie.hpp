#pragma once
#include <cplib/collections/private/binary_trie_base.hpp>
#include <cplib/utils/backwards_index.hpp>

namespace cplib {
class BinaryTrie : public detail::BinaryTrieBase {
public:
    // ビット幅hの空のTrieを構築する。O(1)。
    explicit BinaryTrie(Int h) : BinaryTrieBase(h) {
    }

    // 全ノードを複製する。確保済みノード数をMとしてO(M)。
    BinaryTrie(const BinaryTrie &other) : BinaryTrieBase(other.h) {
        root = clone(other.root);
    }

    // 全ノードを複製して代入する。両者の確保済みノード数の合計をMとしてO(M)。
    BinaryTrie &operator=(const BinaryTrie &other) {
        if (this != &other) {
            root = clone(other.root);
            h = other.h;
        }
        return *this;
    }

    // Trieをムーブ構築する。O(1)。
    BinaryTrie(BinaryTrie &&) = default;
    // Trieをムーブ代入する。代入先の確保済みノード数をMとしてO(M)。
    BinaryTrie &operator=(BinaryTrie &&) = default;

    // 要素の多重度を更新する。O(H)。空になったノードも元実装と同様に残す。
    void incl(Int x, Int v = 1) {
        assert(x >= 0);
        auto now = root;
        now->value += v;
        for (Int i = h - 1; i >= 0; --i) {
            auto &child = ((UInt(x) >> i) & 1) ? now->one : now->zero;
            if (!child)
                child = std::make_shared<detail::BinaryTrieNode>();
            now = child;
            now->value += v;
        }
    }

    // xの多重度をv減らす。ビット幅Hに対しO(H+1)。
    void excl(Int x, Int v = 1) {
        assert(x >= 0);
        auto now = root;
        now->value -= v;
        assert(now->value >= 0);
        for (Int i = h - 1; i >= 0; --i) {
            auto &child = ((UInt(x) >> i) & 1) ? now->one : now->zero;
            if (!child)
                child = std::make_shared<detail::BinaryTrieNode>();
            now = child;
            now->value -= v;
            assert(now->value >= 0);
        }
    }

    // XOR後の値がx以下である要素数を、重複込みで返します。O(H)。
    Int upperBound(Int x, Int xor_value = 0) const {
        return bound(x, xor_value, true);
    }

    // XOR後の値がx未満である要素数を、重複込みで返します。O(H)。
    Int lowerBound(Int x, Int xor_value = 0) const {
        return bound(x, xor_value, false);
    }

    using BinaryTrieBase::operator[];

    // 末尾からk番目の要素を返す。ビット幅Hに対しO(H+1)。
    Int operator[](BackwardsIndex k) const {
        return get_kth(size() - k.value);
    }
};

// ビット幅hの空のTrieを生成する。O(1)。
inline auto initBineryTrie(Int h) {
    return BinaryTrie(h);
}

// ビット幅hの空のTrieを生成する。O(1)。
inline auto initBinaryTrie(Int h) {
    return BinaryTrie(h);
}

// xの多重度をv増やす。ビット幅Hに対しO(H+1)。
inline void incl(BinaryTrie &t, Int x, Int v = 1) {
    t.incl(x, v);
}

// xの多重度をv減らす。ビット幅Hに対しO(H+1)。
inline void excl(BinaryTrie &t, Int x, Int v = 1) {
    t.excl(x, v);
}

// xの多重度を返す。ビット幅Hに対しO(H+1)。
inline Int count(const BinaryTrie &t, Int x) {
    return t.count(x);
}

// xが存在するか返す。ビット幅Hに対しO(H+1)。
inline bool contains(const BinaryTrie &t, Int x) {
    return t.contains(x);
}

// 値にxをXORした順でk番目の元の値を返す。O(H+1)。
inline Int get_kth(const BinaryTrie &t, Int k, Int x = 0) {
    return t.get_kth(k, x);
}

// maskとのXORがx以下の要素数を返す。O(H+1)。
inline Int upperBound(const BinaryTrie &t, Int x, Int mask = 0) {
    return t.upperBound(x, mask);
}

// maskとのXORがx未満の要素数を返す。O(H+1)。
inline Int lowerBound(const BinaryTrie &t, Int x, Int mask = 0) {
    return t.lowerBound(x, mask);
}

// 値と多重度の組を文字列化する。確保済みノード数をM、出力長をLとしてO(M+L)。
inline auto to_string(const BinaryTrie &t) {
    return t.to_string();
}
}
