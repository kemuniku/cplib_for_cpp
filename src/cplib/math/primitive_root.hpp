#pragma once
#include <cplib/math/primefactor.hpp>
#include <cplib/math/powmod.hpp>

namespace cplib {
// 素数法の原始根を乱択で探す。素因数分解後、1回の試行は異なる素因数の数を k として O(k log p)。
// 試行回数の期待値は (p-1)/φ(p-1)。素因数分解の計算量は primefactor を参照。
inline Int primitive_root(Int p) {
    assert(isprime(p));
    auto pf = primefactor(p - 1);
    pf.erase(std::unique(pf.begin(), pf.end()), pf.end());
    for (;;) {
        Int a = std::uniform_int_distribution<Int>(1, p - 1)(detail::factor_rng);
        bool ok = true;
        for (Int q : pf)
            if (powmod(a, (p - 1) / q, p) == 1) {
                ok = false;
                break;
            }
        if (ok)
            return a;
    }
}
}
