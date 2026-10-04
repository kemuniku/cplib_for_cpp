#pragma once
#include <cplib/common.hpp>
#include <cplib/collections/private/wordsizetree_avx2_kernel.hpp>
#include <string_view>

namespace cplib {
inline constexpr Int WordsizeTreeAvx2Capacity = Int(1) << 24;

class WordsizeTreeAvx2 {
    std::array<std::uint64_t, 1 << 18> leaf_{};
    std::array<std::uint64_t, 1 << 10> middle_{};
    std::array<std::uint64_t, 4> top_{};

    void initialize(const void *p, std::size_t n, unsigned char one) {
        assert(n <= std::size_t(WordsizeTreeAvx2Capacity));
        if (n)
            detail::wordsize_avx2::wst_init(p, n, one, leaf_.data(), middle_.data(), top_.data());
    }

public:
    // 256 分岐・3段。元の AVX2 比較・movemask カーネルを使用。
    WordsizeTreeAvx2() = default;

    explicit WordsizeTreeAvx2(std::string_view v) {
        initialize(v.data(), v.size(), '1');
    }

    explicit WordsizeTreeAvx2(std::span<const bool> v) {
        static_assert(sizeof(bool) == 1);
        initialize(v.data(), v.size(), 1);
    }

    explicit WordsizeTreeAvx2(const std::vector<bool> &v) {
        std::vector<unsigned char> bytes(v.size());
        for (std::size_t i = 0; i < v.size(); ++i)
            bytes[i] = v[i];
        initialize(bytes.data(), bytes.size(), 1);
    }

    void incl(Int x) {
        assert(0 <= x && x < WordsizeTreeAvx2Capacity);
        detail::wordsize_avx2::wst_incl(leaf_.data(), middle_.data(), top_.data(), x);
    }

    void excl(Int x) {
        assert(0 <= x && x < WordsizeTreeAvx2Capacity);
        detail::wordsize_avx2::wst_excl(leaf_.data(), middle_.data(), top_.data(), x);
    }

    bool operator[](Int x) const {
        assert(0 <= x && x < WordsizeTreeAvx2Capacity);
        return (leaf_[x >> 6] >> (x & 63)) & 1;
    }

    // x以上の最小の要素を返し、存在しなければ-1を返します。
    Int ge(Int x) const {
        if (x >= WordsizeTreeAvx2Capacity)
            return -1;
        return detail::wordsize_avx2::wst_ge(leaf_.data(), middle_.data(), top_.data(),
                                             std::max(x, Int(0)));
    }

    // x以下の最大の要素を返し、存在しなければ-1を返します。
    Int le(Int x) const {
        if (x < 0)
            return -1;
        return detail::wordsize_avx2::wst_le(leaf_.data(), middle_.data(), top_.data(),
                                             std::min(x, WordsizeTreeAvx2Capacity - 1));
    }
};

inline WordsizeTreeAvx2 initWordsizeTree() {
    return {};
}

inline WordsizeTreeAvx2 initWordsizeTree(std::string_view v) {
    return WordsizeTreeAvx2(v);
}

inline WordsizeTreeAvx2 initWordsizeTree(std::span<const bool> v) {
    return WordsizeTreeAvx2(v);
}

inline WordsizeTreeAvx2 initWordsizeTree(const std::vector<bool> &v) {
    return WordsizeTreeAvx2(v);
}
}
