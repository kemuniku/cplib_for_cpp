#pragma once
#include <cplib/collections/private/binary_trie_base.hpp>
#include <cplib/utils/backwards_index.hpp>

namespace cplib {
class BinaryTrie : public detail::BinaryTrieBase {
public:
    explicit BinaryTrie(Int h) : BinaryTrieBase(h) {
    }

    BinaryTrie(const BinaryTrie &other) : BinaryTrieBase(other.h) {
        root = clone(other.root);
    }

    BinaryTrie &operator=(const BinaryTrie &other) {
        if (this != &other) {
            root = clone(other.root);
            h = other.h;
        }
        return *this;
    }

    BinaryTrie(BinaryTrie &&) = default;
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

    Int operator[](BackwardsIndex k) const {
        return get_kth(len() - k.value);
    }
};

inline auto initBineryTrie(Int h) {
    return BinaryTrie(h);
}

inline auto initBinaryTrie(Int h) {
    return BinaryTrie(h);
}

inline void incl(BinaryTrie &t, Int x, Int v = 1) {
    t.incl(x, v);
}

inline void excl(BinaryTrie &t, Int x, Int v = 1) {
    t.excl(x, v);
}

inline Int count(const BinaryTrie &t, Int x) {
    return t.count(x);
}

inline bool contains(const BinaryTrie &t, Int x) {
    return t.contains(x);
}

inline Int get_kth(const BinaryTrie &t, Int k, Int x = 0) {
    return t.get_kth(k, x);
}

inline Int upperBound(const BinaryTrie &t, Int x, Int mask = 0) {
    return t.upperBound(x, mask);
}

inline Int lowerBound(const BinaryTrie &t, Int x, Int mask = 0) {
    return t.lowerBound(x, mask);
}

inline Int len(const BinaryTrie &t) {
    return t.len();
}

inline auto to_string(const BinaryTrie &t) {
    return t.to_string();
}
}
