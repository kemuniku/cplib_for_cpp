#pragma once
#ifndef __STDC_WANT_IEC_60559_TYPES_EXT__
#define __STDC_WANT_IEC_60559_TYPES_EXT__ 1
#endif
#include <cplib/common.hpp>
#include <cplib/math/floating_point.hpp>
#include <cerrno>
#include <climits>
#include <cstdlib>
#include <istream>
#include <ostream>
#include <stdexcept>
#include <string_view>

// GNU/Linux (glibc 2.26以降) 用。Boost・libquadmathは不要。
// 幾何などのテンプレートより先にincludeする。
// using namespace cplib; の後はsqrt(x), sin(x), cin >> x, cout << xなどが使える。
// std::sqrt等への追加はしない。cplib::sqrtまたは非修飾名で呼ぶ。
namespace cplib {
using float128 = __float128;

template <> inline constexpr bool isFloatingPoint<float128> = true;

// stdを拡張せず、cplib::NumericLimits<float128>でIEEE binary128の特性を公開する。
template <> struct NumericLimits<float128> : std::numeric_limits<long double> {
    static constexpr bool is_specialized = true;
    static constexpr int digits = 113, digits10 = 33, max_digits10 = 36;
    static constexpr int radix = 2;
    static constexpr int min_exponent = -16381, min_exponent10 = -4931;
    static constexpr int max_exponent = 16384, max_exponent10 = 4932;
    static constexpr bool is_signed = true, is_integer = false, is_exact = false;
    static constexpr bool is_bounded = true, is_modulo = false, is_iec559 = true;
    static constexpr bool has_infinity = true, has_quiet_NaN = true, has_signaling_NaN = true;
    static constexpr std::float_denorm_style has_denorm = std::denorm_present;
    static constexpr bool has_denorm_loss = false;
    static constexpr std::float_round_style round_style = std::round_to_nearest;
private:
    static constexpr float128 bits(__uint128_t value) noexcept {
        return std::bit_cast<float128>(value);
    }
public:
    static constexpr float128 (min)() noexcept { return bits(__uint128_t(1) << 112); }
    static constexpr float128 (max)() noexcept { return bits((__uint128_t(0x7fff) << 112) - 1); }
    static constexpr float128 lowest() noexcept { return -(max)(); }
    static constexpr float128 epsilon() noexcept { return bits(__uint128_t(16383 - 112) << 112); }
    static constexpr float128 round_error() noexcept { return float128(0.5); }
    static constexpr float128 infinity() noexcept { return bits(__uint128_t(0x7fff) << 112); }
    static constexpr float128 quiet_NaN() noexcept { return bits(__uint128_t(0xffff) << 111); }
    static constexpr float128 signaling_NaN() noexcept { return bits((__uint128_t(0x7fff) << 112) | 1); }
    static constexpr float128 denorm_min() noexcept { return bits(1); }
};

// 数値型との変換。O(1)。
template <class T> constexpr float128 to_float128(T x) { return static_cast<float128>(x); }
inline double to_float(float128 x) { return static_cast<double>(x); }
inline int cmp(float128 x, float128 y) { return x < y ? -1 : x == y ? 0 : 1; }

namespace detail {
template <class T> inline constexpr bool isNativeFloat128 = std::is_same_v<T, float128>;
template <class... T> inline constexpr bool float128Arguments =
    (isNativeFloat128<T> || ...) && ((std::is_arithmetic_v<T> || isNativeFloat128<T>) && ...);
}

// 制約付きテンプレートにし、標準ライブラリに同じ型の非テンプレート関数があれば優先する。
// 数学関数の精度・定義域・特殊値はglibcのf128関数に従う。各呼び出しO(1)。
#define CPLIB_FLOAT128_UNARY(name) \
    template <class T> requires detail::isNativeFloat128<T> \
    inline float128 name(T x) { return ::name##f128(x); } \
    using ::name##l; \
    template <class T> requires detail::isNativeFloat128<T> \
    inline float128 name##l(T x) { return ::name##f128(x); }

CPLIB_FLOAT128_UNARY(fabs)
CPLIB_FLOAT128_UNARY(sqrt)
CPLIB_FLOAT128_UNARY(cbrt)
CPLIB_FLOAT128_UNARY(exp)
CPLIB_FLOAT128_UNARY(exp2)
CPLIB_FLOAT128_UNARY(expm1)
CPLIB_FLOAT128_UNARY(log)
CPLIB_FLOAT128_UNARY(log2)
CPLIB_FLOAT128_UNARY(log10)
CPLIB_FLOAT128_UNARY(log1p)
CPLIB_FLOAT128_UNARY(sin)
CPLIB_FLOAT128_UNARY(cos)
CPLIB_FLOAT128_UNARY(tan)
CPLIB_FLOAT128_UNARY(asin)
CPLIB_FLOAT128_UNARY(acos)
CPLIB_FLOAT128_UNARY(atan)
CPLIB_FLOAT128_UNARY(sinh)
CPLIB_FLOAT128_UNARY(cosh)
CPLIB_FLOAT128_UNARY(tanh)
CPLIB_FLOAT128_UNARY(asinh)
CPLIB_FLOAT128_UNARY(acosh)
CPLIB_FLOAT128_UNARY(atanh)
CPLIB_FLOAT128_UNARY(floor)
CPLIB_FLOAT128_UNARY(ceil)
CPLIB_FLOAT128_UNARY(trunc)
CPLIB_FLOAT128_UNARY(round)
CPLIB_FLOAT128_UNARY(rint)
CPLIB_FLOAT128_UNARY(nearbyint)
CPLIB_FLOAT128_UNARY(logb)
CPLIB_FLOAT128_UNARY(erf)
CPLIB_FLOAT128_UNARY(erfc)
CPLIB_FLOAT128_UNARY(tgamma)
CPLIB_FLOAT128_UNARY(lgamma)
#undef CPLIB_FLOAT128_UNARY

template <class T> requires detail::isNativeFloat128<T>
inline float128 abs(T x) { return __builtin_fabsf128(x); }

#define CPLIB_FLOAT128_BINARY(name) \
    template <class T, class U> requires detail::float128Arguments<T, U> \
    inline float128 name(T x, U y) { return ::name##f128(float128(x), float128(y)); } \
    using ::name##l; \
    template <class T, class U> requires detail::float128Arguments<T, U> \
    inline float128 name##l(T x, U y) { return ::name##f128(float128(x), float128(y)); }

CPLIB_FLOAT128_BINARY(pow)
CPLIB_FLOAT128_BINARY(hypot)
CPLIB_FLOAT128_BINARY(atan2)
CPLIB_FLOAT128_BINARY(fmod)
CPLIB_FLOAT128_BINARY(remainder)
CPLIB_FLOAT128_BINARY(copysign)
CPLIB_FLOAT128_BINARY(nextafter)
CPLIB_FLOAT128_BINARY(fdim)
CPLIB_FLOAT128_BINARY(fmin)
CPLIB_FLOAT128_BINARY(fmax)
#undef CPLIB_FLOAT128_BINARY

template <class T, class U, class V> requires detail::float128Arguments<T, U, V>
inline float128 hypot(T x, U y, V z) {
    return ::hypotf128(::hypotf128(float128(x), float128(y)), float128(z));
}
template <class T, class U, class V> requires detail::float128Arguments<T, U, V>
inline float128 fma(T x, U y, V z) { return ::fmaf128(float128(x), float128(y), float128(z)); }
using ::fmal;
template <class T, class U, class V> requires detail::float128Arguments<T, U, V>
inline float128 fmal(T x, U y, V z) { return ::fmaf128(float128(x), float128(y), float128(z)); }

#define CPLIB_FLOAT128_INTEGER_RESULT(name, result) \
    template <class T> requires detail::isNativeFloat128<T> \
    inline result name(T x) { return ::name##f128(x); } \
    using ::name##l; \
    template <class T> requires detail::isNativeFloat128<T> \
    inline result name##l(T x) { return ::name##f128(x); }
CPLIB_FLOAT128_INTEGER_RESULT(lround, long)
CPLIB_FLOAT128_INTEGER_RESULT(llround, long long)
CPLIB_FLOAT128_INTEGER_RESULT(lrint, long)
CPLIB_FLOAT128_INTEGER_RESULT(llrint, long long)
CPLIB_FLOAT128_INTEGER_RESULT(ilogb, int)
#undef CPLIB_FLOAT128_INTEGER_RESULT

#define CPLIB_FLOAT128_EXTRA_ARG(name, Arg) \
    template <class T> requires detail::isNativeFloat128<T> \
    inline float128 name(T x, Arg arg) { return ::name##f128(x, arg); } \
    using ::name##l; \
    template <class T> requires detail::isNativeFloat128<T> \
    inline float128 name##l(T x, Arg arg) { return ::name##f128(x, arg); }
CPLIB_FLOAT128_EXTRA_ARG(frexp, int *)
CPLIB_FLOAT128_EXTRA_ARG(ldexp, int)
CPLIB_FLOAT128_EXTRA_ARG(modf, float128 *)
CPLIB_FLOAT128_EXTRA_ARG(scalbn, int)
CPLIB_FLOAT128_EXTRA_ARG(scalbln, long)
#undef CPLIB_FLOAT128_EXTRA_ARG

template <class T, class U> requires detail::float128Arguments<T, U>
inline float128 remquo(T x, U y, int *quotient) { return ::remquof128(float128(x), float128(y), quotient); }
using ::remquol;
template <class T, class U> requires detail::float128Arguments<T, U>
inline float128 remquol(T x, U y, int *quotient) { return ::remquof128(float128(x), float128(y), quotient); }

#define CPLIB_FLOAT128_CLASSIFY(name) \
    template <class T> requires detail::isNativeFloat128<T> \
    inline bool name(T x) { return __builtin_##name(x); }
CPLIB_FLOAT128_CLASSIFY(isfinite)
CPLIB_FLOAT128_CLASSIFY(isinf)
CPLIB_FLOAT128_CLASSIFY(isnan)
CPLIB_FLOAT128_CLASSIFY(isnormal)
CPLIB_FLOAT128_CLASSIFY(signbit)
#undef CPLIB_FLOAT128_CLASSIFY

template <class T> requires detail::isNativeFloat128<T>
inline int fpclassify(T x) { return __builtin_fpclassify(FP_NAN, FP_INFINITE, FP_NORMAL, FP_SUBNORMAL, FP_ZERO, x); }

#define CPLIB_FLOAT128_COMPARE(name) \
    template <class T, class U> requires detail::float128Arguments<T, U> \
    inline bool name(T x, U y) { return __builtin_##name(float128(x), float128(y)); }
CPLIB_FLOAT128_COMPARE(isgreater)
CPLIB_FLOAT128_COMPARE(isgreaterequal)
CPLIB_FLOAT128_COMPARE(isless)
CPLIB_FLOAT128_COMPARE(islessequal)
CPLIB_FLOAT128_COMPARE(islessgreater)
CPLIB_FLOAT128_COMPARE(isunordered)
#undef CPLIB_FLOAT128_COMPARE

// Cロケールの数値文字列を直接128ビットで読む。O(文字数)。不正な文字列はinvalid_argument。
// 末尾の余分な文字は許さない。overflow・ゼロへのunderflowはout_of_range。
inline float128 parseFloat128(const std::string &text) {
    char *end;
    errno = 0;
    float128 result = ::strtof128(text.c_str(), &end);
    if (end == text.c_str() || end != text.c_str() + text.size())
        throw std::invalid_argument("float128の数値文字列ではありません");
    if (errno == ERANGE && (result == 0 || !cplib::isfinite(result)))
        throw std::out_of_range("float128の表現範囲外です");
    return result;
}

namespace detail {
// strfromf128はprintfの符号・幅フラグを受け付けないので、数値の変換だけに使う。
inline std::string formatFloat128(float128 value, char conversion, std::streamsize precision) {
    std::string format = "%";
    if (conversion != 'a' && conversion != 'A')
        format += "." + std::to_string(precision < 0 ? 6 : precision);
    format += conversion;
    int size = ::strfromf128(nullptr, 0, format.c_str(), value);
    if (size < 0)
        throw std::runtime_error("float128の文字列化に失敗しました");
    std::string result(size, '\0');
    ::strfromf128(result.data(), result.size() + 1, format.c_str(), value);
    return result;
}
}

// 有効桁数precisionの十進文字列に変換。既定の36桁で往復変換可能。O(出力文字数)。
template <class T> requires detail::isNativeFloat128<T>
inline std::string to_string(T value, int precision = NumericLimits<float128>::max_digits10) {
    return detail::formatFloat128(value, 'g', precision);
}

// 空白区切りの数値トークンを読む。long doubleを経由しない。O(トークン長)。
// トークン全体を消費し、不正な入力ではfailbitを立てる。
inline std::istream &operator>>(std::istream &in, float128 &value) {
    std::string token;
    if (!(in >> token))
        return in;
    try {
        value = parseFloat128(token);
    } catch (const std::invalid_argument &) {
        in.setstate(std::ios_base::failbit);
    } catch (const std::out_of_range &) {
        in.setstate(std::ios_base::failbit);
    }
    return in;
}

// Cロケール表記で出力する。精度・fixed/scientific/hexfloat・符号・幅・埋め文字に対応。
// std::localeの小数点・桁区切り変更は対象外。時間・追加領域O(出力文字数)。
inline std::ostream &operator<<(std::ostream &out, float128 value) {
    std::ostream::sentry sentry(out);
    if (!sentry)
        return out;
    try {
        auto flags = out.flags();
        auto field = flags & std::ios_base::floatfield;
        char conversion = field == std::ios_base::fixed ? 'f' :
            field == std::ios_base::scientific ? 'e' :
            field == std::ios_base::floatfield ? 'a' : 'g';
        if (flags & std::ios_base::uppercase)
            conversion -= 'a' - 'A';
        auto text = detail::formatFloat128(value, conversion, out.precision());
        if ((flags & std::ios_base::showpoint) && cplib::isfinite(value)) {
            auto end = text.find_first_of(field == std::ios_base::floatfield ? "pP" : "eE");
            if (end == std::string::npos)
                end = text.size();
            int significant = 0;
            for (std::size_t i = 0; i < end; ++i)
                if (text[i] >= '0' && text[i] <= '9' && (significant || text[i] != '0'))
                    ++significant;
            significant = std::max(1, significant);
            if (text.find('.') == std::string::npos)
                text.insert(end++, 1, '.');
            if (field == std::ios_base::fmtflags(0)) {
                auto precision = std::max<std::streamsize>(1, out.precision() < 0 ? 6 : out.precision());
                if (precision > significant)
                    text.insert(end, precision - significant, '0');
            }
        }
        if ((flags & std::ios_base::showpos) && !cplib::signbit(value))
            text.insert(text.begin(), '+');
        auto width = out.width(0);
        if (width > std::streamsize(text.size())) {
            auto count = std::size_t(width) - text.size();
            auto adjust = flags & std::ios_base::adjustfield;
            if (adjust == std::ios_base::left)
                text.append(count, out.fill());
            else {
                std::size_t prefix = 0;
                if (adjust == std::ios_base::internal && (text[0] == '+' || text[0] == '-'))
                    prefix = 1;
                text.insert(prefix, count, out.fill());
            }
        }
        if (out.rdbuf()->sputn(text.data(), std::streamsize(text.size())) != std::streamsize(text.size()))
            out.setstate(std::ios_base::badbit);
    } catch (...) {
        out.setstate(std::ios_base::badbit);
    }
    return out;
}
}

// 組込み型にはcplibのADLが働かないため、通常のcin/cout用にも公開する。
using cplib::operator>>;
using cplib::operator<<;
