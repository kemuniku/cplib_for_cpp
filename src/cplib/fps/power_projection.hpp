#pragma once
#include <cplib/fps/formal_power_series.hpp>

namespace cplib {
// [x^(f.size()-1)]f(x)^i g(x), i=0..m。x方向のBostan–Moriを平坦化した巡回積で処理。
template <Modint T>
std::vector<T> powerProjection(const std::vector<T> &f, const std::vector<T> &g, Int m) {
    assert(!f.empty() && m >= 0);
    Int n = f.size() - 1, xStride = 1, yDegree = 1;
    while (xStride < n + 1)
        xStride *= 2;
    std::vector<T> p(xStride), q(xStride);
    for (Int i = 0; i <= n; ++i) {
        if (i < Int(g.size()))
            p[i] = g[i];
        q[i] = -f[i];
    }
    while (n > 0) {
        Int wideStride = 2 * xStride;
        std::vector<T> expandedP(yDegree * wideStride), expandedQ((yDegree + 1) * wideStride);
        for (Int y = 0; y < yDegree; ++y)
            for (Int x = 0; x <= n; ++x) {
                expandedP[y * wideStride + x] = p[y * xStride + x];
                expandedQ[y * wideStride + x] = q[y * xStride + x];
            }
        expandedQ[yDegree * wideStride] = 1;
        auto negativeQ = expandedQ;
        for (Int y = 0; y <= yDegree; ++y)
            for (Int x = 1; x < wideStride; x += 2)
                negativeQ[y * wideStride + x] = -negativeQ[y * wideStride + x];
        Int cycleLength = 2 * yDegree * wideStride;
        auto productP = convolutionCyclicPowerOfTwo(expandedP, negativeQ, cycleLength),
             productQ = convolutionCyclicPowerOfTwo(expandedQ, negativeQ, cycleLength);
        productQ[0] -= T(1);
        Int nextStride = xStride / 2;
        std::vector<T> nextP(2 * yDegree * nextStride), nextQ(2 * yDegree * nextStride);
        for (Int y = 0; y < 2 * yDegree; ++y)
            for (Int x = 0; x <= n / 2; ++x) {
                Int base = y * wideStride + 2 * x;
                nextP[y * nextStride + x] = productP[base + (n & 1)];
                nextQ[y * nextStride + x] = productQ[base];
            }
        p = std::move(nextP);
        q = std::move(nextQ);
        n /= 2;
        xStride = nextStride;
        yDegree *= 2;
    }
    p.resize(yDegree);
    q.resize(yDegree);
    q.push_back(T(1));
    std::reverse(p.begin(), p.end());
    std::reverse(q.begin(), q.end());
    return prefix(p * inv(q, m + 1), m + 1);
}

template <Modint T> auto powerProjection(const std::vector<T> &f, Int m) {
    return powerProjection(f, std::vector<T>{T(1)}, m);
}

template <Modint T> auto powerProjection(const std::vector<T> &f, const std::vector<T> &g) {
    return powerProjection(f, g, Int(f.size()) - 1);
}

template <Modint T> auto powerProjection(const std::vector<T> &f) {
    return powerProjection(f, Int(f.size()) - 1);
}

// [x^i]f(x)^i g(x)。O((m+1)log²(m+2))。
// [x^i] f(x)^i g(x) (i = 0, 1, ..., m) を列挙する。
// fは空でなく、mは非負とする。入力の範囲外の係数は零として扱う。
// f[0] != 0の場合はFPSのpowと同じくm + 1が法以下である必要がある。
// NTTを使える場合O((m + 1) log^2(m + 2))時間。
template <Modint T>
std::vector<T> powerProjectionDiagonal(const std::vector<T> &f, const std::vector<T> &g, Int m) {
    assert(!f.empty() && m >= 0);
    if (f[0].val() == 0) {
        std::vector<T> out(m + 1);
        out[0] = coefficient(g, 0);
        T linear = coefficient(f, 1);
        for (Int i = 1; i <= m; ++i)
            out[i] = out[i - 1] * linear;
        return out;
    }
    Int size = m + 1;
    auto base = prefix(f, size);
    auto weight = prefix(pow(base, m, size) * prefix(g, size), size);
    auto inverse = inv(base, m);
    std::vector<T> shifted(size);
    for (Int i = 1; i <= m; ++i)
        shifted[i] = inverse[i - 1];
    auto out = powerProjection(shifted, weight, m);
    std::reverse(out.begin(), out.end());
    return out;
}

template <Modint T> auto powerProjectionDiagonal(const std::vector<T> &f, Int m) {
    return powerProjectionDiagonal(f, std::vector<T>{T(1)}, m);
}

template <Modint T> auto powerProjectionDiagonal(const std::vector<T> &f, const std::vector<T> &g) {
    return powerProjectionDiagonal(f, g, Int(f.size()) - 1);
}

template <Modint T> auto powerProjectionDiagonal(const std::vector<T> &f) {
    return powerProjectionDiagonal(f, Int(f.size()) - 1);
}
}
