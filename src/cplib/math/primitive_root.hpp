#pragma once
#include <cplib/math/primefactor.hpp>
#include <cplib/math/powmod.hpp>

namespace cplib {
// 素数法の原始根を乱択で探す。元実装と同じ素因数判定を用いる。
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
