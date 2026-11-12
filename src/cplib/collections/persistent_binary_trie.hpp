#pragma once
#include <cplib/collections/private/binary_trie_base.hpp>

namespace cplib {
// 値は非負整数とし、キーの比較には下位Hビットを使う。Hは構築時のビット幅。
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
    // ビット幅hの空の永続Trieを構築する。O(1)。
    explicit PersistentBinaryTrie(Int h) : BinaryTrieBase(h) {
    }

    // xの多重度をv増やした新版を返す。ビット幅Hとして時間・追加領域O(H+1)。
    auto incl(Int x, Int v = 1) const {
        return updated(x, v, false);
    }

    // xの多重度をv減らした新版を返す。O(H+1)時間・追加領域。
    auto excl(Int x, Int v = 1) const {
        return updated(x, -v, true);
    }

    // xの多重度をvにした新版を返す。O(H+1)時間・追加領域。
    auto set_value(Int x, Int v) const {
        return incl(x, v - count(x));
    }

    // 計算量はビット幅Hに対しO(H+1)。
    // x以下の要素数
    Int upperBound(Int x) const {
        return bound(x, 0, true);
    }

    // 計算量はビット幅Hに対しO(H+1)。
    // x未満の要素数
    Int lowerBound(Int x) const {
        return bound(x, 0, false);
    }
};

// ビット幅hの空の永続Trieを生成する。O(1)。
inline auto initPersistentBineryTrie(Int h) {
    return PersistentBinaryTrie(h);
}

// ビット幅hの空の永続Trieを生成する。O(1)。
inline auto initPersistentBinaryTrie(Int h) {
    return PersistentBinaryTrie(h);
}

// xの多重度をv増やした新版を返す。O(H+1)時間・追加領域。
inline auto incl(const PersistentBinaryTrie &t, Int x, Int v = 1) {
    return t.incl(x, v);
}

// xの多重度をv減らした新版を返す。O(H+1)時間・追加領域。
inline auto excl(const PersistentBinaryTrie &t, Int x, Int v = 1) {
    return t.excl(x, v);
}

// xの多重度をvにした新版を返す。O(H+1)時間・追加領域。
inline auto set_value(const PersistentBinaryTrie &t, Int x, Int v) {
    return t.set_value(x, v);
}

// xの多重度を返す。O(H+1)。
inline Int count(const PersistentBinaryTrie &t, Int x) {
    return t.count(x);
}

// xが存在するか返す。O(H+1)。
inline bool contains(const PersistentBinaryTrie &t, Int x) {
    return t.contains(x);
}

// 値にxをXORした順でk番目の元の値を返す。O(H+1)。
inline Int get_kth(const PersistentBinaryTrie &t, Int k, Int x = 0) {
    return t.get_kth(k, x);
}

// x以下の要素数を返す。O(H+1)。
inline Int upperBound(const PersistentBinaryTrie &t, Int x) {
    return t.upperBound(x);
}

// x未満の要素数を返す。O(H+1)。
inline Int lowerBound(const PersistentBinaryTrie &t, Int x) {
    return t.lowerBound(x);
}

// 値と多重度の組を昇順に列挙する。確保済みノード数をMとしてO(M)。
inline auto RLE(const PersistentBinaryTrie &t) {
    return t.RLE();
}

// 値と多重度の組を文字列化する。確保済みノード数をM、出力長をLとしてO(M+L)。
inline auto to_string(const PersistentBinaryTrie &t) {
    return t.to_string();
}
}
