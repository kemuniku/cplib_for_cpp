#pragma once
#include <cplib/common.hpp>

namespace cplib {
// Σ floor((a*i+b)/m) (0 <= i < n)をO(log m)時間、O(1)空間で返す。
// n, a, b >= 0、m > 0、答えがIntに収まることを前提とする。
inline Int floor_sum(Int nn, Int mm, Int aa, Int bb) {
    assert(nn >= 0 && mm > 0 && aa >= 0 && bb >= 0);
    __int128 n = nn, m = mm, a = aa, b = bb, answer = 0;
    while (true) {
        if (a >= m) {
            answer += n * (n - 1) / 2 * (a / m);
            a %= m;
        }
        if (b >= m) {
            answer += n * (b / m);
            b %= m;
        }
        assert(answer <= std::numeric_limits<Int>::max());
        __int128 y = a * n + b;
        if (y < m)
            break;
        n = y / m;
        b = y % m;
        std::swap(m, a);
    }
    return static_cast<Int>(answer);
}
}
