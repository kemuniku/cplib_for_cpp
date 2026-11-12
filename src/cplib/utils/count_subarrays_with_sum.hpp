#pragma once
#include <cplib/common.hpp>
#include <concepts>
#include <stdexcept>
#include <unordered_map>

namespace cplib {
// 和がxの空でない連続部分列を数える。期待時間 O(N)、領域 O(N)。
// 累積和は符号に応じた64ビット整数で保持し、和・個数の桁あふれはoverflow_error。
template <std::integral T>
Int count_subarrays_with_sum(const std::vector<T> &a, T x) {
    using Sum = std::conditional_t<std::is_signed_v<T>, Int, UInt>;
    Sum prefix = 0, target = x;
    std::unordered_map<Sum, Int> frequency{{0, 1}};
    Int result = 0;
    for (T value : a) {
        Sum next;
        if (__builtin_add_overflow(prefix, Sum(value), &next))
            throw std::overflow_error("累積和が64ビットの範囲外です");
        prefix = next;
        Sum wanted;
        if (!__builtin_sub_overflow(prefix, target, &wanted)) {
            auto it = frequency.find(wanted);
            if (it != frequency.end() && __builtin_add_overflow(result, it->second, &result))
                throw std::overflow_error("個数がIntの範囲外です");
        }
        auto &count = frequency[prefix];
        if (count == std::numeric_limits<Int>::max())
            throw std::overflow_error("頻度がIntの範囲外です");
        ++count;
    }
    return result;
}
}
