#pragma once
#include <cplib/math/monoid_floor_sum.hpp>

namespace cplib {
// result[j][k] = Σ(i=0..<n) i^j * floor((a*i+b)/m)^k (0<=j<=p, 0<=k<=q)。0^0=1。
// n, a, b, p, q >= 0、m > 0、a*n+bがIntの範囲内であること。
// TはIntから変換できる可換環で、既定値が零の型を指定する。除算は不要。
// 時間O((p+1)(q+1)(p+q+2)log(m+1)log(n+a+b+2))、空間O((p+q+2)^2)。
// 整数型では結果だけでなく途中の環演算も型の範囲内に収まる必要がある。
template <class T>
std::vector<std::vector<T>> generalizedFloorSumTable(Int n, Int m, Int a, Int b, Int p, Int q) {
    assert(n >= 0 && m > 0 && a >= 0 && b >= 0 && p >= 0 && q >= 0);
    assert(n == 0 || a <= (std::numeric_limits<Int>::max() - b) / n);
    using Table = std::vector<std::vector<T>>;
    auto zero = [&]() { return Table(p + 1, std::vector<T>(q + 1)); };
    if (n == 0)
        return zero();
    T one = 1;
    Int degree = std::max(p, q);
    Table binom(degree + 1);
    for (Int j = 0; j <= degree; ++j) {
        binom[j].resize(j + 1);
        binom[j][0] = binom[j][j] = one;
        for (Int k = 1; k < j; ++k)
            binom[j][k] = binom[j - 1][k - 1] + binom[j - 1][k];
    }

    struct Moment {
        T dx{}, dy{};
        Table sums;
    };

    auto combine = [&](const Moment &l, const Moment &r) {
        std::vector<T> xp(p + 1), yp(q + 1);
        xp[0] = yp[0] = one;
        for (Int j = 1; j <= p; ++j)
            xp[j] = xp[j - 1] * l.dx;
        for (Int k = 1; k <= q; ++k)
            yp[k] = yp[k - 1] * l.dy;
        auto shifted = zero();
        for (Int j = 0; j <= p; ++j)
            for (Int s = 0; s <= j; ++s) {
                T c = binom[j][s] * xp[j - s];
                for (Int k = 0; k <= q; ++k)
                    shifted[j][k] = shifted[j][k] + c * r.sums[s][k];
            }
        Moment result{l.dx + r.dx, l.dy + r.dy, zero()};
        for (Int k = 0; k <= q; ++k)
            for (Int t = 0; t <= k; ++t) {
                T c = binom[k][t] * yp[k - t];
                for (Int j = 0; j <= p; ++j)
                    result.sums[j][k] = result.sums[j][k] + c * shifted[j][t];
            }
        for (Int j = 0; j <= p; ++j)
            for (Int k = 0; k <= q; ++k)
                result.sums[j][k] = result.sums[j][k] + l.sums[j][k];
        return result;
    };
    Moment x{one, T{}, zero()}, y{T{}, one, zero()}, e{T{}, T{}, zero()};
    x.sums[0][0] = one;
    return monoidFloorSum(n, m, a, b, x, y, combine, e).sums;
}

// Σ(i=0..<n) i^p * floor((a*i+b)/m)^q を返す。条件・計算量は generalizedFloorSumTable と同じ。
template <class T> T generalizedFloorSum(Int n, Int m, Int a, Int b, Int p, Int q) {
    return generalizedFloorSumTable<T>(n, m, a, b, p, q)[p][q];
}
}
