#pragma once
#include <cplib/common.hpp>
#include <cplib/collections/private/wordsizetree_avx2_kernel.hpp>
#include <string_view>

namespace cplib {
inline constexpr Int WordsizeTreeAvx2Capacity = Int(1) << 24;

class WordsizeTreeAvx2 {
    Int count_ = 0;
    std::array<std::uint64_t, 1 << 18> leaf_{};
    std::array<std::uint64_t, 1 << 10> middle_{};
    std::array<std::uint64_t, 4> top_{};

    void initialize(const void *p, std::size_t n, unsigned char one) {
        assert(n <= std::size_t(WordsizeTreeAvx2Capacity));
        if (n) {
            detail::wordsize_avx2::wst_init(p, n, one, leaf_.data(), middle_.data(), top_.data());
            for (std::size_t i = 0; i < (n + 63) / 64; ++i)
                count_ += std::popcount(leaf_[i]);
        }
    }

public:
    // 固定容量2^24の空集合を作る。固定サイズの領域を初期化するためO(1)。
    // 256 分岐・3段。元の AVX2 比較・movemask カーネルを使用。
    WordsizeTreeAvx2() = default;

    // 文字'1'の位置を集合に登録する。入力長NとしてO(N+1)。
    explicit WordsizeTreeAvx2(std::string_view v) {
        initialize(v.data(), v.size(), '1');
    }

    // 真の位置を集合に登録する。入力長NとしてO(N+1)。
    explicit WordsizeTreeAvx2(std::span<const bool> v) {
        static_assert(sizeof(bool) == 1);
        initialize(v.data(), v.size(), 1);
    }

    // 真の位置を集合に登録する。入力長NとしてO(N+1)。
    explicit WordsizeTreeAvx2(const std::vector<bool> &v) {
        std::vector<unsigned char> bytes(v.size());
        for (std::size_t i = 0; i < v.size(); ++i)
            bytes[i] = v[i];
        initialize(bytes.data(), bytes.size(), 1);
    }

    // 登録されている要素数。固定容量とは異なる。O(1)。
    Int size() const {
        return count_;
    }

    // xを集合に追加する。O(1)。
    void incl(Int x) {
        assert(0 <= x && x < WordsizeTreeAvx2Capacity);
        if (!(*this)[x])
            ++count_;
        detail::wordsize_avx2::wst_incl(leaf_.data(), middle_.data(), top_.data(), x);
    }

    // xを集合から削除する。O(1)。
    void excl(Int x) {
        assert(0 <= x && x < WordsizeTreeAvx2Capacity);
        if ((*this)[x])
            --count_;
        detail::wordsize_avx2::wst_excl(leaf_.data(), middle_.data(), top_.data(), x);
    }

    // xが集合に含まれるか返す。O(1)。
    bool operator[](Int x) const {
        assert(0 <= x && x < WordsizeTreeAvx2Capacity);
        return (leaf_[x >> 6] >> (x & 63)) & 1;
    }

    // 計算量O(1)。
    // x以上の最小の要素を返し、存在しなければ-1を返します。
    Int ge(Int x) const {
        if (x >= WordsizeTreeAvx2Capacity)
            return -1;
        return detail::wordsize_avx2::wst_ge(leaf_.data(), middle_.data(), top_.data(),
                                             std::max(x, Int(0)));
    }

    // 計算量O(1)。
    // x以下の最大の要素を返し、存在しなければ-1を返します。
    Int le(Int x) const {
        if (x < 0)
            return -1;
        return detail::wordsize_avx2::wst_le(leaf_.data(), middle_.data(), top_.data(),
                                             std::min(x, WordsizeTreeAvx2Capacity - 1));
    }
};

// 固定容量2^24の空集合を作る。固定サイズの領域を初期化するためO(1)。
inline WordsizeTreeAvx2 initWordsizeTreeAvx2() {
    return {};
}

// 文字'1'の位置を集合に登録する。入力長NとしてO(N+1)。
inline WordsizeTreeAvx2 initWordsizeTreeAvx2(std::string_view v) {
    return WordsizeTreeAvx2(v);
}

// 真の位置を集合に登録する。入力長NとしてO(N+1)。
inline WordsizeTreeAvx2 initWordsizeTreeAvx2(std::span<const bool> v) {
    return WordsizeTreeAvx2(v);
}

// 真の位置を集合に登録する。入力長NとしてO(N+1)。
inline WordsizeTreeAvx2 initWordsizeTreeAvx2(const std::vector<bool> &v) {
    return WordsizeTreeAvx2(v);
}
}
