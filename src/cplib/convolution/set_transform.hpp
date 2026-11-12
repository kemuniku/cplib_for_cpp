#pragma once
#include <cplib/common.hpp>

namespace cplib {
namespace detail {
template <bool Superset, bool Inverse, class T>
void setTransform(std::vector<T> &a) {
    assert(a.empty() || std::has_single_bit(a.size()));
    for (std::size_t bit = 1; bit < a.size(); bit *= 2) {
        for (std::size_t base = 0; base < a.size(); base += bit * 2) {
            for (std::size_t j = 0; j < bit; ++j) {
                auto &target = a[base + j + (Superset ? 0 : bit)];
                const auto &source = a[base + j + (Superset ? bit : 0)];
                if constexpr (Inverse)
                    target -= source;
                else
                    target += source;
            }
        }
    }
}
}
// 各集合に部分集合の和を格納する。長さNは空または2冪。時間 O(N log N)、追加領域 O(1)。
template <class T> void subsetZetaTransform(std::vector<T> &a) {
    detail::setTransform<false, false>(a);
}
// 部分集合ゼータ変換を戻す。時間 O(N log N)、追加領域 O(1)。
template <class T> void subsetMobiusTransform(std::vector<T> &a) {
    detail::setTransform<false, true>(a);
}
// 各集合に上位集合の和を格納する。長さNは空または2冪。時間 O(N log N)、追加領域 O(1)。
template <class T> void supersetZetaTransform(std::vector<T> &a) {
    detail::setTransform<true, false>(a);
}
// 上位集合ゼータ変換を戻す。時間 O(N log N)、追加領域 O(1)。
template <class T> void supersetMobiusTransform(std::vector<T> &a) {
    detail::setTransform<true, true>(a);
}
}
