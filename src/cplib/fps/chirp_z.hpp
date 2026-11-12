#pragma once
#include <cplib/convolution/convolution.hpp>

namespace cplib {
// f(a*r^i), 0<=i<mを返す。可逆なrなら時間 O((n+m)log(n+m))、非可逆なら O(nm)。領域 O(n+m)。
template <class T>
std::vector<T> multipointEvaluationGeometric(const std::vector<T> &f, T a, T r, Int m) {
    assert(m >= 0);
    const Int n = f.size();
    std::vector<T> result(m);
    if (m == 0 || n == 0)
        return result;
    if (n == 1 || a == T(0)) {
        std::fill(result.begin(), result.end(), f[0]);
        return result;
    }
    auto [g, inverse] = inv_gcd(r.val(), T::umod());
    if (g != 1 || m == 1 || r == T(0) || r == T(1) || r == T(-1)) {
        T point = a;
        T atA{}, atMinusA{};
        for (Int j = n; j-- > 0;) {
            atA = atA * a + f[j];
            atMinusA = atMinusA * (-a) + f[j];
        }
        for (Int i = 0; i < m; ++i) {
            if (r == T(0))
                result[i] = i == 0 ? atA : f[0];
            else if (r == T(1) || m == 1)
                result[i] = atA;
            else if (r == T(-1))
                result[i] = i % 2 == 0 ? atA : atMinusA;
            else {
                for (Int j = n; j-- > 0;)
                    result[i] = result[i] * point + f[j];
                point *= r;
            }
        }
        return result;
    }
    std::vector<T> chirp(n + m - 1), inverseChirp(std::max(n, m)), left(n);
    chirp[0] = inverseChirp[0] = 1;
    T power = 1, inversePower = 1;
    for (Int i = 1; i < Int(chirp.size()); ++i) {
        chirp[i] = chirp[i - 1] * power;
        power *= r;
    }
    for (Int i = 1; i < Int(inverseChirp.size()); ++i) {
        inverseChirp[i] = inverseChirp[i - 1] * inversePower;
        inversePower *= T(inverse);
    }
    power = 1;
    for (Int j = 0; j < n; ++j) {
        left[n - 1 - j] = f[j] * power * inverseChirp[j];
        power *= a;
    }
    auto product = convolution(left, chirp);
    for (Int i = 0; i < m; ++i)
        result[i] = product[n - 1 + i] * inverseChirp[i];
    return result;
}
// 等比数列上の多点評価を行う。時間 O((n+m)log(n+m))、非可逆なrでは O(nm)。領域 O(n+m)。
template <class T>
std::vector<T> chirpZ(const std::vector<T> &f, T r, Int m, T a = T(1)) {
    return multipointEvaluationGeometric(f, a, r, m);
}
}
