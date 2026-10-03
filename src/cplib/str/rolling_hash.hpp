#pragma once
#include <cplib/str/private/hash_string_base.hpp>

namespace cplib {
namespace detail {
// この版だけ、元ソースどおり mod と等しい値を 0 に正規化しない。
inline UInt modulo(UInt x) {
    x = (x >> 61) + (x & cplib::detail::hash_string_mod);
    if (x > cplib::detail::hash_string_mod)
        x -= cplib::detail::hash_string_mod;
    return x;
}

inline UInt power(UInt a, UInt n) {
    UInt result = 1;
    while (n) {
        if (n & 1)
            result = modulo(cplib::detail::hash_string_mul(result, a));
        a = modulo(cplib::detail::hash_string_mul(a, a));
        n >>= 1;
    }
    return result;
}

inline UInt base = 0, inverse = 0;
inline bool initialized = false;

inline void initialize(UInt maxa, Int seed) {
    if (initialized)
        return;
    std::mt19937_64 rng(seed == -1 ? std::random_device{}() : UInt(seed));
    UInt k = 0, b = 1;
    while (b <= maxa || std::gcd(cplib::detail::hash_string_mod - 1, k) != 1) {
        k = std::uniform_int_distribution<UInt>(0, cplib::detail::hash_string_mod - 1)(rng);
        b = power(37, k);
    }
    base = b;
    inverse = power(base, cplib::detail::hash_string_mod - 2);
    initialized = true;
}
}

template <class Sequence> class RollingHash {
    Sequence s_;
    std::vector<UInt> accum_, powers_, inverse_powers_;

public:
    explicit RollingHash(Sequence s) : s_(std::move(s)) {
        build();
    }

    void build(UInt maxa = 1000000000, Int seed = -1) {
        detail::initialize(maxa, seed);
        accum_.assign(s_.size() + 1, 0);
        powers_.assign(s_.size() + 1, 1);
        inverse_powers_.assign(s_.size() + 1, 1);
        for (std::size_t i = 0; i < s_.size(); ++i) {
            UInt value;
            if constexpr (std::is_same_v<typename Sequence::value_type, char>)
                value = static_cast<unsigned char>(s_[i]);
            else
                value = UInt(s_[i]);
            accum_[i + 1] =
                detail::modulo(accum_[i] + cplib::detail::hash_string_mul(value, powers_[i]));
            powers_[i + 1] =
                detail::modulo(cplib::detail::hash_string_mul(powers_[i], detail::base));
            inverse_powers_[i + 1] =
                detail::modulo(cplib::detail::hash_string_mul(inverse_powers_[i], detail::inverse));
        }
    }

    template <class A, class B> UInt query(ClosedSlice<A, B> range) const {
        Int l = range.a, r = Int(range.b) + 1;
        assert(0 <= l && l < Int(accum_.size()) && 0 <= r && r < Int(accum_.size()));
        return detail::modulo(cplib::detail::hash_string_mul(
            accum_[r] + cplib::detail::hash_string_mod - accum_[l], inverse_powers_[l]));
    }
};

template <class T> auto initRollingHash(const std::vector<T> &s) {
    return RollingHash<std::vector<T>>(s);
}

inline auto initRollingHash(const std::string &s) {
    return RollingHash<std::string>(s);
}
}
