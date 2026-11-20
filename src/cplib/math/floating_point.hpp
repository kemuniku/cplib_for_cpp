#pragma once
#include <cmath>
#include <limits>
#include <type_traits>

namespace cplib {
// 標準型はそのまま利用する。拡張数値型の対応ヘッダーで特殊化できる。
template <class T> struct NumericLimits : std::numeric_limits<T> {};
template <class T> inline constexpr bool isFloatingPoint =
    std::is_floating_point_v<T> ||
    (NumericLimits<T>::is_specialized && !NumericLimits<T>::is_integer &&
     !NumericLimits<T>::is_exact);
template <class T> concept FloatingPoint = isFloatingPoint<T>;

// 拡張型の同名関数を先に宣言すれば、後続のテンプレートでも選択できる。
// int128.hppにもabsがあるため、stdのGNU拡張版を直接usingせず転送する。
template <class T> auto abs(T value) -> decltype(std::abs(value)) {
    return std::abs(value);
}
using std::fabs;
using std::sqrt;
using std::cbrt;
using std::hypot;
using std::pow;
using std::exp;
using std::exp2;
using std::expm1;
using std::log;
using std::log2;
using std::log10;
using std::log1p;
using std::sin;
using std::cos;
using std::tan;
using std::asin;
using std::acos;
using std::atan;
using std::atan2;
using std::sinh;
using std::cosh;
using std::tanh;
using std::asinh;
using std::acosh;
using std::atanh;
using std::floor;
using std::ceil;
using std::trunc;
using std::round;
using std::rint;
using std::nearbyint;
using std::lround;
using std::llround;
using std::lrint;
using std::llrint;
using std::fmod;
using std::remainder;
using std::remquo;
using std::copysign;
using std::nextafter;
using std::fdim;
using std::fmin;
using std::fmax;
using std::fma;
using std::frexp;
using std::ldexp;
using std::modf;
using std::scalbn;
using std::scalbln;
using std::ilogb;
using std::logb;
using std::erf;
using std::erfc;
using std::tgamma;
using std::lgamma;
using std::isfinite;
using std::isinf;
using std::isnan;
using std::isnormal;
using std::signbit;
using std::fpclassify;
using std::isgreater;
using std::isgreaterequal;
using std::isless;
using std::islessequal;
using std::islessgreater;
using std::isunordered;
}
