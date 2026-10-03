#pragma once
#include <cplib/math/primefactor.hpp>

namespace cplib {
// 約数を列挙する。小さい入力はO(sqrt(x))、大きい入力は素因数分解後O(約数個数)。
inline std::vector<Int> divisor(Int x, bool sorted = true) {
    std::vector<Int> result;
    if (x <= 1000000) {
        for (Int i = 1; i <= x / i; ++i)
            if (x % i == 0) {
                result.push_back(i);
                if (i != x / i)
                    result.push_back(x / i);
            }
    } else {
        auto factor = primefactor_table(x);
        result.push_back(1);
        for (auto [p, count] : factor) {
            std::size_t old = result.size();
            Int mult = 1;
            for (Int k = 1; k <= count; ++k) {
                mult *= p;
                for (std::size_t j = 0; j < old; ++j)
                    result.push_back(result[j] * mult);
            }
        }
    }
    if (sorted)
        std::sort(result.begin(), result.end());
    return result;
}
}
