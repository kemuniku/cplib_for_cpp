#pragma once
#include <cplib/math/isqrt.hpp>

namespace cplib {
struct FloorSumProgression {
    Int x, y, dx, dy, n;
};

// Ak+B mod MをO(sqrt(N))個の等差数列に分解する。O(sqrt(N))時間・領域。
inline std::vector<FloorSumProgression> sqrt_heuristic_for_floor_sum(Int A, Int B, Int N, Int M) {
    A %= M;
    B %= M;
    Int D = isqrt(N), bidx = -1, bval = M;
    std::vector<FloorSumProgression> result;
    for (Int i = 1; i <= D; ++i) {
        Int val = A * i % M;
        val = std::min(val, M - val);
        if (bval > val) {
            bval = val;
            bidx = i;
        }
    }
    if ((bidx * A % M) > M - (bidx * A % M)) {
        for (auto v : sqrt_heuristic_for_floor_sum((M - A) % M, M - 1 - B, N, M))
            result.push_back({v.x, M - 1 - v.y, v.dx, -v.dy, v.n});
        return result;
    }
    for (Int g = 0; g < bidx; ++g) {
        Int a = A * bidx % M, b = (A * g + B) % M, n = (N - g + bidx - 1) / bidx,
            x = (a * (n - 1) + b) / M + 1, l = 0;
        for (Int k = 0; k < x; ++k) {
            Int r = k + 1 == x ? n : (M * (k + 1) - b + a - 1) / a;
            result.push_back({bidx * l + g, (a * l + b) % M, bidx, bidx * A % M, r - l});
            l = r;
        }
    }
    return result;
}
}
