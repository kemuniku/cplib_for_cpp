#pragma once
#include <cplib/common.hpp>
#include <concepts>
#include <functional>
#include <optional>
#include <string_view>
#include <unordered_map>

namespace cplib {
template <class S> using DigitDpTransition = std::function<std::optional<S>(Int, S, Int)>;
template <class S> using DigitDpDirectTransition = std::function<S(Int, S, Int)>;
template <class S> using DigitDpInitialState = std::function<S(Int, Int)>;
template <class S> using DigitDpOptionalInitialState = std::function<std::optional<S>(Int, Int)>;
namespace detail {
template <class T> struct DigitDpStateType { using type = T; };
template <class T> struct DigitDpStateType<std::optional<T>> { using type = T; };
template <class T> using DigitDpState = typename DigitDpStateType<T>::type;
template <class S, class T> void digitDpAccumulate(std::unordered_map<S, T> &table,
                                                  const S &state, const T &count, const T &zero) {
    auto [it, inserted] = table.try_emplace(state, zero);
    it->second = it->second + count;
}
inline void validateDigits(const std::vector<Int> &digits, Int base) {
    assert(!digits.empty() && base >= 2);
    for (Int digit : digits)
        assert(0 <= digit && digit < base);
}
}
// 0埋めした全桁を遷移させ、上限以下の非負整数を最終状態ごとに数える。
// nextは状態またはoptional<状態>を返す。桁数L、最大状態数K、基数Bとして期待時間 O(LKB)、領域 O(K)。
template <class S, class Next, class T>
    requires (!std::is_invocable_v<S, Int, Int>)
std::unordered_map<S, T> digitDp(const std::vector<Int> &digits, S initialState, Next next,
                                T zero, T one, Int base = 10) {
    detail::validateDigits(digits, base);
    std::array<std::unordered_map<S, T>, 2> dp;
    dp[1].emplace(initialState, one);
    for (Int pos = 0; pos < Int(digits.size()); ++pos) {
        std::array<std::unordered_map<S, T>, 2> target;
        for (Int tight = 0; tight < 2; ++tight) {
            Int limit = tight ? digits[pos] : base - 1;
            for (const auto &[state, count] : dp[tight]) {
                for (Int digit = 0; digit <= limit; ++digit) {
                    std::optional<S> out = next(pos, state, digit);
                    if (out)
                        detail::digitDpAccumulate(target[tight && digit == digits[pos]], *out, count, zero);
                }
            }
        }
        dp = std::move(target);
    }
    for (const auto &[state, count] : dp[1])
        detail::digitDpAccumulate(dp[0], state, count, zero);
    return std::move(dp[0]);
}
// 先頭桁と実際の桁数から状態を作る。桁合わせの0は飛ばし、整数0はinitial(0,1)で扱う。
// 初期化・遷移は状態またはoptional<状態>を返す。期待時間 O(LKB)、領域 O(K)。
template <class Initial, class Next, class T>
    requires std::is_invocable_v<Initial, Int, Int>
auto digitDp(const std::vector<Int> &digits, Initial initial, Next next,
              T zero, T one, Int base = 10) {
    using S = detail::DigitDpState<std::invoke_result_t<Initial, Int, Int>>;
    detail::validateDigits(digits, base);
    std::array<std::unordered_map<S, T>, 2> dp;
    bool unstartedTight = true;
    for (Int pos = 0; pos < Int(digits.size()); ++pos) {
        std::array<std::unordered_map<S, T>, 2> target;
        for (Int tight = 0; tight < 2; ++tight) {
            Int limit = tight ? digits[pos] : base - 1;
            for (const auto &[state, count] : dp[tight]) {
                for (Int digit = 0; digit <= limit; ++digit) {
                    std::optional<S> out = next(pos, state, digit);
                    if (out)
                        detail::digitDpAccumulate(target[tight && digit == digits[pos]], *out, count, zero);
                }
            }
        }
        Int limit = unstartedTight ? digits[pos] : base - 1;
        for (Int digit = 0; digit <= limit; ++digit) {
            if (digit == 0 && pos + 1 < Int(digits.size()))
                continue;
            std::optional<S> out = initial(digit, Int(digits.size()) - pos);
            if (out)
                detail::digitDpAccumulate(target[unstartedTight && digit == digits[pos]], *out, one, zero);
        }
        unstartedTight = unstartedTight && digits[pos] == 0;
        dp = std::move(target);
    }
    for (const auto &[state, count] : dp[1])
        detail::digitDpAccumulate(dp[0], state, count, zero);
    return std::move(dp[0]);
}
// 整数上限をbase進数に直して桁DPする。期待時間 O(LKB)、領域 O(K+L)。
template <std::integral Bound, class Initial, class Next, class T>
auto digitDp(Bound upper, Initial initial, Next next, T zero, T one, Int base = 10) {
    assert(upper >= 0 && base >= 2);
    UInt value = upper;
    std::vector<Int> digits;
    do {
        digits.push_back(value % UInt(base));
        value /= UInt(base);
    } while (value != 0);
    std::reverse(digits.begin(), digits.end());
    return digitDp(digits, initial, next, zero, one, base);
}
// 十進文字列上限で桁DPする。期待時間 O(LK*10)、領域 O(K+L)。
template <class Initial, class Next, class T>
auto digitDp(std::string_view upper, Initial initial, Next next, T zero, T one) {
    std::vector<Int> digits;
    for (char c : upper) {
        assert('0' <= c && c <= '9');
        digits.push_back(c - '0');
    }
    return digitDp(digits, initial, next, zero, one, 10);
}
}
