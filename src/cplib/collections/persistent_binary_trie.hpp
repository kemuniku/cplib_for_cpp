#pragma once
#include <cplib/collections/private/binary_trie_base.hpp>

namespace cplib {
class PersistentBinaryTrie : public detail::BinaryTrieBase {
    PersistentBinaryTrie updated(Int x, Int v, bool remove) const {
        assert(x >= 0);
        PersistentBinaryTrie out = *this;
        out.root = std::make_shared<detail::BinaryTrieNode>(*root);
        auto now = out.root;
        now->value += v;
        if (remove)
            assert(now->value >= 0);
        for (Int i = h - 1; i >= 0; --i) {
            auto &child = ((UInt(x) >> i) & 1) ? now->one : now->zero;
            if (remove)
                assert(child);
            child = child ? std::make_shared<detail::BinaryTrieNode>(*child)
                          : std::make_shared<detail::BinaryTrieNode>();
            now = child;
            now->value += v;
            if (remove)
                assert(now->value >= 0);
        }
        return out;
    }

public:
    explicit PersistentBinaryTrie(Int h) : BinaryTrieBase(h) {
    }

    // 更新経路のみ複製する永続版。時間・新規領域O(H)。
    auto incl(Int x, Int v = 1) const {
        return updated(x, v, false);
    }

    auto excl(Int x, Int v = 1) const {
        return updated(x, -v, true);
    }

    auto set_value(Int x, Int v) const {
        return incl(x, v - count(x));
    }

    // x以下の要素数
    Int upperBound(Int x) const {
        return bound(x, 0, true);
    }

    // x未満の要素数
    Int lowerBound(Int x) const {
        return bound(x, 0, false);
    }
};

inline auto initPersistentBineryTrie(Int h) {
    return PersistentBinaryTrie(h);
}

inline auto initPersistentBinaryTrie(Int h) {
    return PersistentBinaryTrie(h);
}

inline auto incl(const PersistentBinaryTrie &t, Int x, Int v = 1) {
    return t.incl(x, v);
}

inline auto excl(const PersistentBinaryTrie &t, Int x, Int v = 1) {
    return t.excl(x, v);
}

inline auto set_value(const PersistentBinaryTrie &t, Int x, Int v) {
    return t.set_value(x, v);
}

inline Int count(const PersistentBinaryTrie &t, Int x) {
    return t.count(x);
}

inline bool contains(const PersistentBinaryTrie &t, Int x) {
    return t.contains(x);
}

inline Int get_kth(const PersistentBinaryTrie &t, Int k, Int x = 0) {
    return t.get_kth(k, x);
}

inline Int upperBound(const PersistentBinaryTrie &t, Int x) {
    return t.upperBound(x);
}

inline Int lowerBound(const PersistentBinaryTrie &t, Int x) {
    return t.lowerBound(x);
}

inline auto RLE(const PersistentBinaryTrie &t) {
    return t.RLE();
}

inline auto to_string(const PersistentBinaryTrie &t) {
    return t.to_string();
}
}
