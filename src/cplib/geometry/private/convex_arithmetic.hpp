#pragma once
#include <cplib/geometry/base.hpp>
#include <cplib/math/bigint.hpp>
#include <cplib/math/fractions.hpp>

namespace cplib::detail {

template <class T>
inline constexpr bool convexInteger = std::is_integral_v<T> || std::is_same_v<T, __int128> ||
                                      std::is_same_v<T, __uint128_t> || std::is_same_v<T, BigInt>;

template <class T> using ConvexWide =
    std::conditional_t<convexInteger<T>, std::conditional_t<(sizeof(T) <= 4), __int128, BigInt>, T>;

template <class T> using ConvexField = std::conditional_t<convexInteger<T>, Fraction<BigInt>, T>;

inline BigInt convexBigInt(__uint128_t value) {
    BigInt result(static_cast<UInt>(value));
    UInt high = static_cast<UInt>(value >> 64);
    if (high != 0)
        result += BigInt(high) << 64;
    return result;
}

inline BigInt convexBigInt(__int128 value) {
    auto unsigned_value = static_cast<__uint128_t>(value);
    return value < 0 ? -convexBigInt(__uint128_t(0) - unsigned_value) : convexBigInt(unsigned_value);
}

template <class T> BigInt convexBigInt(const T &value) {
    return BigInt(value);
}

template <class R, class T, class Id> Point<R> convexPoint(const Point<T, Id> &p) {
    if constexpr ((std::is_same_v<R, Fraction<BigInt>> || std::is_same_v<R, BigInt>) && convexInteger<T>) {
        return {convexBigInt(p.x), convexBigInt(p.y)};
    } else {
        return {R(p.x), R(p.y)};
    }
}

template <class T> Point<T> convexNormal(Point<T> p) {
    return {-p.y, p.x};
}

}
