#pragma once
#include <cplib/common.hpp>
#include <type_traits>

namespace cplib {
// 128ビット積で整数偏角を比較する。正のx軸から反時計回り。計算量 O(1)。
// 整数幅: 入力はInt（64ビット）、積は__int128へ昇格。積の差を作らず比較する。
// 全Int範囲（-2^63も含む）で積の絶対値<=2^126、符号込み128ビットで安全。10^9も可。
inline Int argcmp(std::pair<Int, Int> l, std::pair<Int, Int> r) {
    bool a = l.second > 0 || (l.second == 0 && l.first > 0),
         b = r.second > 0 || (r.second == 0 && r.first > 0);
    if (a != b)
        return b < a ? -1 : 1;
    __int128_t x = __int128_t(l.second) * r.first, y = __int128_t(l.first) * r.second;
    return x < y ? -1 : x > y ? 1 : 0;
}

// 整数専用APIへ、浮動小数点やIntに収まらない整数型を暗黙変換させない。
namespace detail {
template <class T> inline constexpr bool argcmpCoordinate =
    std::is_integral_v<T> && std::is_signed_v<T> && sizeof(T) <= sizeof(Int);
}
template <class A, class B, class C, class D>
    requires (!(detail::argcmpCoordinate<A> && detail::argcmpCoordinate<B> &&
                detail::argcmpCoordinate<C> && detail::argcmpCoordinate<D>))
Int argcmp(std::pair<A, B>, std::pair<C, D>) = delete;

// N個の非零ベクトルを正のx軸から反時計回りに安定ソートする。計算量 O(N log N)。
// 整数幅: 入力はInt（64ビット）、積は__int128へ昇格。積の差を作らず比較する。
// 全Int範囲（-2^63も含む）で積の絶対値<=2^126、符号込み128ビットで安全。10^9も可。
inline void argsort(std::vector<std::pair<Int, Int>> &x) {
    std::stable_sort(x.begin(), x.end(), [](auto l, auto r) { return argcmp(l, r) < 0; });
}

// N個の非零ベクトルを正のx軸から反時計回りに偏角ソートした列を返す。計算量 O(N log N)。
// 整数幅: 入力はInt（64ビット）、積は__int128へ昇格。積の差を作らず比較する。
// 全Int範囲（-2^63も含む）で積の絶対値<=2^126、符号込み128ビットで安全。10^9も可。
inline auto argsorted(std::vector<std::pair<Int, Int>> x) {
    argsort(x);
    return x;
}
}
