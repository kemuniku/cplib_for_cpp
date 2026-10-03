#pragma once
#include <cplib/modint/modint.hpp>

namespace cplib {
namespace detail {
constexpr Int ntt_powmod(Int a, Int x, Int p) {
    Int ans = 1;
    while (x > 0) {
        if (x & 1)
            ans = ans * a % p;
        a = a * a % p;
        x >>= 1;
    }
    return ans % p;
}

constexpr Int ntt_primitive_root(Int p) {
    std::array<Int, 32> pf{};
    Int size = 0, x = p - 1;
    for (Int i = 2; i * i <= x; ++i)
        if (x % i == 0) {
            pf[size++] = i;
            while (x % i == 0)
                x /= i;
        }
    if (x != 1)
        pf[size++] = x;
    for (Int a = 3; a < p; ++a) {
        bool valid = true;
        for (Int i = 0; i < size; ++i)
            if (ntt_powmod(a, (p - 1) / pf[i], p) == 1) {
                valid = false;
                break;
            }
        if (valid)
            return a;
    }
    return 0;
}

struct NttConfig {
    Int sum_e = 0, m = 0, forth_root = 0, forth_root_inv = 0, primitive_root = 0;
    std::array<Int, 30> rate2{}, rate3{}, irate2{}, irate3{};
};

constexpr NttConfig initNttConfig(Int m) {
    NttConfig c;
    c.m = m;
    c.primitive_root = ntt_primitive_root(m);
    Int e = 0, prod = 1, iprod = 1;
    auto nth = [&](Int n) { return ntt_powmod(c.primitive_root, (m - 1) / n, m); };
    while (e + 2 <= std::countr_zero(UInt(m - 1)) - 1) {
        Int root = nth(Int(1) << (e + 2)), iroot = ntt_powmod(root, m - 2, m);
        c.rate2[e] = root * iprod % m;
        c.irate2[e] = iroot * prod % m;
        prod = prod * root % m;
        iprod = iprod * iroot % m;
        ++e;
    }
    e = 0;
    prod = iprod = 1;
    while (e + 3 <= std::countr_zero(UInt(m - 1))) {
        Int root = nth(Int(1) << (e + 3)), iroot = ntt_powmod(root, m - 2, m);
        c.rate3[e] = root * iprod % m;
        c.irate3[e] = iroot * prod % m;
        prod = prod * root % m;
        iprod = iprod * iroot % m;
        ++e;
    }
    c.sum_e = e + 2;
    c.forth_root = nth(4);
    c.forth_root_inv = ntt_powmod(c.forth_root, m - 2, m);
    return c;
}

inline NttConfig dynamic_ntt_config;
template <class T> struct NttConfigFor;

template <std::uint32_t M, bool Dynamic> struct NttConfigFor<BasicMontgomeryModint<M, Dynamic>> {
    static const NttConfig &get() {
        if constexpr (Dynamic) {
            if (dynamic_ntt_config.m != BasicMontgomeryModint<M, Dynamic>::umod())
                dynamic_ntt_config = initNttConfig(BasicMontgomeryModint<M, Dynamic>::umod());
            return dynamic_ntt_config;
        } else {
            static constexpr auto config = initNttConfig(M);
            return config;
        }
    }
};

template <std::uint32_t M, bool Dynamic> struct NttConfigFor<BasicBarrettModint<M, Dynamic>> {
    static const NttConfig &get() {
        if constexpr (Dynamic) {
            if (dynamic_ntt_config.m != BasicBarrettModint<M, Dynamic>::umod())
                dynamic_ntt_config = initNttConfig(BasicBarrettModint<M, Dynamic>::umod());
            return dynamic_ntt_config;
        } else {
            static constexpr auto config = initNttConfig(M);
            return config;
        }
    }
};
}

// radix-4を基本とする元のNTT。静的法は係数をコンパイル時に構築。O(N log N)。
template <class T>
    requires(BarrettModint<T> || MontgomeryModint<T>)
void ntt(std::vector<T> &f) {
    Int n = f.size();
    if (n <= 1)
        return;
    const auto c = detail::NttConfigFor<T>::get();
    assert(std::has_single_bit(UInt(n)));
    for (Int width = n; width > 1;) {
        if (width == 2) {
            Int offset = width >> 1;
            T root = 1;
            for (Int top = 0; top < n; top += width) {
                for (Int i = top; i < top + offset; ++i) {
                    T c0 = f[i], c1 = f[i + offset] * root;
                    f[i] = c0 + c1;
                    f[i + offset] = c0 - c1;
                }
                root *= T(c.rate2[std::countr_zero(~UInt(top / width))]);
            }
            width >>= 1;
        } else {
            Int offset = width >> 2;
            T root = 1;
            for (Int top = 0; top < n; top += width) {
                T root2 = root * root, root3 = root * root2;
                for (Int i = top; i < top + offset; ++i) {
                    T c0 = f[i], c1 = f[i + offset] * root, c2 = f[i + offset * 2] * root2,
                      c3 = f[i + offset * 3] * root3;
                    T a = c0 + c2, b = c0 - c2, d = c1 + c3, e = (c1 - c3) * T(c.forth_root);
                    f[i] = a + d;
                    f[i + offset] = a - d;
                    f[i + offset * 2] = b + e;
                    f[i + offset * 3] = b - e;
                }
                root *= T(c.rate3[std::countr_zero(~UInt(top / width))]);
            }
            width >>= 2;
        }
    }
}

template <class T>
    requires(BarrettModint<T> || MontgomeryModint<T>)
void intt(std::vector<T> &f) {
    Int n = f.size();
    if (n <= 1)
        return;
    const auto c = detail::NttConfigFor<T>::get();
    assert(std::has_single_bit(UInt(n)));
    for (Int width = (std::countr_zero(UInt(n)) % 2 ? 2 : 4); width <= n; width <<= 2) {
        if (width == 2) {
            Int offset = width >> 1;
            T root = 1;
            for (Int top = 0; top < n; top += width) {
                for (Int i = top; i < top + offset; ++i) {
                    T a = f[i], b = f[i + offset];
                    f[i] = a + b;
                    f[i + offset] = (a - b) * root;
                }
                root *= T(c.irate2[std::countr_zero(~UInt(top / width))]);
            }
        } else {
            Int offset = width >> 2;
            T root = 1;
            for (Int top = 0; top < n; top += width) {
                T root2 = root * root, root3 = root * root2;
                for (Int i = top; i < top + offset; ++i) {
                    T c0 = f[i], c1 = f[i + offset], c2 = f[i + offset * 2], c3 = f[i + offset * 3];
                    T a = c0 + c1, b = c0 - c1, d = c2 + c3, e = (c2 - c3) * T(c.forth_root_inv);
                    f[i] = a + d;
                    f[i + offset] = (b + e) * root;
                    f[i + offset * 2] = (a - d) * root2;
                    f[i + offset * 3] = (b - e) * root3;
                }
                root *= T(c.irate3[std::countr_zero(~UInt(top / width))]);
            }
        }
    }
    T ninv = T(n).inv();
    for (auto &x : f)
        x *= ninv;
}
}
