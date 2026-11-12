#pragma once
#include <cplib/math/combination.hpp>
#include <cplib/math/isqrt.hpp>
#include <stdexcept>

namespace cplib {
template <class M> class OnlineCombinationPrefixSum {
    Int maxN_, width_;
    UInt modulus_;
    std::vector<M> fact_, factInv_, powersOfTwo_, halfScaledFact_;
    std::vector<std::vector<M>> samples_;
public:
    // n<=maxNの二項係数和を前計算する。素数法pと0<=maxN<pが必要。素数性は検査しない。
    // 時間 O(N sqrt(N+1))、領域 O(N+1)。使用中に動的法を変更しないこと。
    explicit OnlineCombinationPrefixSum(Int maxN) : maxN_(maxN), modulus_(M::umod()) {
        if (maxN < 0 || UInt(maxN) >= modulus_)
            throw std::invalid_argument("0 <= maxN < 法 が必要です");
        width_ = isqrt(maxN) + 1;
        Combination_Type<M> comb(maxN);
        fact_ = std::move(comb.fact);
        factInv_ = std::move(comb.fact_inv);
        powersOfTwo_.resize(maxN + 1);
        powersOfTwo_[0] = 1;
        for (Int n = 1; n <= maxN; ++n)
            powersOfTwo_[n] = powersOfTwo_[n - 1] * 2;
        if (modulus_ == 2)
            return;
        M half = M((modulus_ + 1) / 2);
        std::vector<M> inversePowers(maxN + 1);
        inversePowers[0] = 1;
        for (Int n = 1; n <= maxN; ++n)
            inversePowers[n] = inversePowers[n - 1] * half;
        halfScaledFact_.resize(maxN + 1);
        for (Int n = 0; n <= maxN; ++n)
            halfScaledFact_[n] = fact_[n] * inversePowers[n] * half;
        samples_.resize(maxN / width_ + 1);
        for (Int row = 0; row < Int(samples_.size()); ++row) {
            Int n = row * width_;
            samples_[row].resize(row + 1);
            samples_[row][row] = powersOfTwo_[n];
            M sum = 0;
            Int nextLow = 0, lowColumn = 0, nextHigh = width_ - 1, highColumn = row - 1;
            for (Int k = 0; k <= n / 2; ++k) {
                sum += factInv_[k] * factInv_[n - k];
                if (k == nextLow) {
                    samples_[row][lowColumn++] = sum * fact_[n];
                    nextLow += width_;
                }
                if (highColumn >= 0 && k == nextHigh) {
                    samples_[row][highColumn--] = powersOfTwo_[n] - sum * fact_[n];
                    nextHigh += width_;
                }
            }
            for (M &value : samples_[row])
                value *= inversePowers[n];
        }
    }
    // Σ(i=0..k) C(n,i)を返す。k<0は0、k>=nは2^n。時間 O(sqrt(N+1))、追加領域 O(1)。
    M prefixSum(Int n, Int k) const {
        if (n < 0 || n > maxN_ || M::umod() != modulus_)
            throw std::invalid_argument("nの範囲または構築後の法が不正です");
        if (k < 0)
            return M(0);
        if (k >= n)
            return powersOfTwo_[n];
        if (k == 0)
            return M(1);
        Int row = std::min<Int>((n + width_ / 2) / width_, samples_.size() - 1);
        Int column = std::min((k + width_ / 2) / width_, row);
        Int baseN = row * width_, baseK = column * width_;
        M result = samples_[row][column], sum = 0;
        if (baseN <= n) {
            for (Int x = baseN; x < n; ++x)
                sum += halfScaledFact_[x] * factInv_[x - baseK];
            result -= sum * factInv_[baseK];
        } else {
            for (Int x = std::max(n, baseK); x < baseN; ++x)
                sum += halfScaledFact_[x] * factInv_[x - baseK];
            result += sum * factInv_[baseK];
        }
        result *= powersOfTwo_[n];
        sum = 0;
        if (baseK <= k) {
            for (Int y = baseK + 1; y <= k; ++y)
                sum += factInv_[y] * factInv_[n - y];
            result += sum * fact_[n];
        } else {
            for (Int y = k + 1; y <= std::min(baseK, n); ++y)
                sum += factInv_[y] * factInv_[n - y];
            result -= sum * fact_[n];
        }
        return result;
    }
    // Σ(i=l..<r) C(n,i)を返す。空・逆転区間は0。時間 O(sqrt(N+1))、追加領域 O(1)。
    M rangeSum(Int n, Int l, Int r) const {
        M right = prefixSum(n, r <= 0 ? -1 : r - 1);
        if (l >= r)
            return M(0);
        return right - prefixSum(n, l <= 0 ? -1 : l - 1);
    }
};
// 二項係数和を前計算する。時間 O(N sqrt(N+1))、領域 O(N+1)。
template <class M> auto initOnlineCombinationPrefixSum(Int maxN) {
    return OnlineCombinationPrefixSum<M>(maxN);
}
}
