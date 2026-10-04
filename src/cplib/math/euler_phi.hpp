#pragma once
#include <cplib/common.hpp>

namespace cplib {
// Eulerのφ関数を試し割りで求める。O(sqrt(n))。
inline Int euler_phi(Int n) {
    Int result = n;
    for (Int i = 2; i <= n / i; ++i)
        if (n % i == 0) {
            result -= result / i;
            while (n % i == 0)
                n /= i;
        }
    if (n > 1)
        result -= result / n;
    return result;
}

// 0からnのφを篩で求める。O(n log log n)時間・O(n)領域。
inline std::vector<Int> euler_phi_list(Int n) {
    std::vector<Int> result(n + 1);
    std::iota(result.begin(), result.end(), 0);
    for (Int i = 2; i <= n; ++i)
        if (result[i] == i)
            for (Int j = i; j <= n; j += i)
                result[j] = result[j] / i * (i - 1);
    return result;
}
}
