#pragma once
#include <cplib/fps/berlekamp_massey.hpp>
#include <cplib/fps/bostan_mori.hpp>

namespace cplib {
// 先頭の数項からBerlekamp--Massey法で漸化式を推定し、Bostan--Mori法で第k項（0-indexed）を求める。
// 法は素数であること。元の数列の漸化式の次数がd以下なら、先頭2d項あれば復元できる。
// k < a.size()ならa[k]を返す。空列・全零列は零数列として扱う。
// 時間計算量O(a.size()^2 + M(d) log(k+1))。dは推定次数、M(d)は長さdの畳み込みの計算量。
template <Modint T> T bmbm(const std::vector<T> &a, Int k) {
    assert(k >= 0);
    if (k < Int(a.size()))
        return a[k];
    auto coefficients = berlekampMassey(a);
    if (coefficients.empty())
        return T(0);
    return linearRecurrenceKth(std::vector<T>(a.begin(), a.begin() + coefficients.size()),
                               coefficients, k);
}
}
