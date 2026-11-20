#pragma once
#include <cplib/geometry/base.hpp>
#include <cplib/math/fractions.hpp>
#include <stdexcept>

namespace cplib::detail {
// 円の厳密計算は指定された基底型を維持する。中間値は全てその型に収まること。
// 以下の内部計算量は基底型の算術・gcdを一操作として数える。
template <class T> struct CircleBaseType { using type = T; };
template <class T> struct CircleBaseType<Fraction<T>> { using type = T; };
template <class T> using CircleBase = typename CircleBaseType<T>::type;
template <class T> inline constexpr bool isCircleFraction = false;
template <class T> inline constexpr bool isCircleFraction<Fraction<T>> = true;

template <class T> Int circleSign(const T &value) {
    return value == T(0) ? 0 : value < T(0) ? -1 : 1;
}
template <class T> T circleGcd(T a, T b) {
    if (a < T(0))
        a = -a;
    if (b < T(0))
        b = -b;
    while (b != T(0)) {
        T remainder = a % b;
        a = std::move(b);
        b = std::move(remainder);
    }
    return a;
}
// 幅の見積もり: |分子|,正の分母<2^u と <2^v の分数同士について、
// 加減算はu+v+1、乗除算はu+vビット（いずれも符号を除く）で中間値まで覆える。
// 整数2倍/4倍はそれぞれ+1/+2ビット。内積・外積は2(u+v)+1、二乗ノルムは4u+1。
// gcdと事前約分はこの上限を増やさない。最後の約分結果だけで幅を判断しない。
// circle.hppのhによる十分条件は、この規則を実装の演算順に適用した安全側の上限。
template <class T> struct CircleRational {
    T num, den;
    CircleRational() : num(0), den(1) {
    }
    CircleRational(T value) : num(std::move(value)), den(1) {
    }
    CircleRational(T numerator, T denominator) : num(std::move(numerator)), den(std::move(denominator)) {
        if (den == T(0))
            throw std::invalid_argument("厳密円の分母は非零である必要があります");
        T divisor = circleGcd(num, den);
        num /= divisor;
        den /= divisor;
        if (den < T(0)) {
            num = -num;
            den = -den;
        }
    }
    CircleRational operator-() const {
        return {-num, den};
    }
    friend CircleRational operator+(const CircleRational &a, const CircleRational &b) {
        T divisor = circleGcd(a.den, b.den);
        T left = a.den / divisor, right = b.den / divisor;
        return {a.num * right + b.num * left, left * b.den};
    }
    friend CircleRational operator-(const CircleRational &a, const CircleRational &b) {
        return a + (-b);
    }
    friend CircleRational operator*(const CircleRational &a, const CircleRational &b) {
        T left = circleGcd(a.num, b.den), right = circleGcd(b.num, a.den);
        return {(a.num / left) * (b.num / right), (a.den / right) * (b.den / left)};
    }
    friend CircleRational operator/(const CircleRational &a, const CircleRational &b) {
        return a * CircleRational(b.den, b.num);
    }
    friend CircleRational operator*(const CircleRational &a, Int value) {
        return a * CircleRational(T(value));
    }
    friend bool operator==(const CircleRational &a, const CircleRational &b) {
        return a.num == b.num && a.den == b.den;
    }
    CircleRational &operator+=(const CircleRational &other) {
        return *this = *this + other;
    }
    CircleRational &operator-=(const CircleRational &other) {
        return *this = *this - other;
    }
    CircleRational &operator*=(const CircleRational &other) {
        return *this = *this * other;
    }
    CircleRational &operator*=(Int other) {
        return *this = *this * other;
    }
};
template <class T> auto circleFraction(const T &value) {
    if constexpr (isCircleFraction<T>)
        return CircleRational<CircleBase<T>>(value.num, value.den);
    else
        return CircleRational<T>(value);
}
template <class T> auto circleExactPoint(const Point<T> &point) {
    return Point{circleFraction(point.x), circleFraction(point.y)};
}
template <class B, class T> Point<CircleRational<B>> circlePointLike(const Point<T> &point) {
    auto convert = [](const auto &value) {
        auto fraction = circleFraction(value);
        return CircleRational<B>(B(fraction.num), B(fraction.den));
    };
    return {convert(point.x), convert(point.y)};
}
template <class T> Fraction<T> circleStored(const CircleRational<T> &value) {
    return Fraction<T>(value.num, value.den, false);
}
template <class T> auto circleStoredPoint(const Point<CircleRational<T>> &point) {
    return Point{circleStored(point.x), circleStored(point.y)};
}
template <class T> T circleStoredLike(const CircleRational<CircleBase<T>> &value) {
    if constexpr (isCircleFraction<T>)
        return T(value.num, value.den, false);
    else {
        if (value.den != CircleBase<T>(1))
            throw std::invalid_argument("整数座標の3点として表現できません");
        return T(value.num);
    }
}
template <class T> std::string circleIntegerText(const T &value) {
    if constexpr (requires { value.str(); }) {
        return value.str();
    } else if constexpr (std::is_same_v<T, __int128_t> || std::is_same_v<T, __uint128_t>) {
        bool negative = value < 0;
        __uint128_t magnitude = static_cast<__uint128_t>(value);
        if (negative)
            magnitude = -magnitude;
        std::string result;
        do {
            result.push_back(char('0' + magnitude % 10));
            magnitude /= 10;
        } while (magnitude != 0);
        if (negative)
            result.push_back('-');
        std::reverse(result.begin(), result.end());
        return result;
    } else {
        return std::to_string(value);
    }
}
template <class F = double, class T> auto circleDecimalHead(const T &value) {
    std::string text = circleIntegerText(value);
    std::size_t begin = text[0] == '-' ? 1 : 0;
    std::size_t count = std::min<std::size_t>(text.size() - begin, cplib::NumericLimits<F>::max_digits10);
    // 指定された浮動小数点型で直接読み、doubleを経由させない。
    F head;
    std::istringstream input(text.substr(begin, count));
    input >> head;
    head /= pow(F(10), F(count - 1));
    return std::pair{head, Int(text.size() - begin - 1)};
}
template <class F = double, class T> F circleApprox(const CircleRational<T> &value, bool squareRoot = false) {
    if (value.num == T(0))
        return 0;
    if (squareRoot && value.num < T(0))
        throw std::invalid_argument("負数の平方根は取得できません");
    auto [numerator, numeratorExponent] = circleDecimalHead<F>(value.num);
    auto [denominator, denominatorExponent] = circleDecimalHead<F>(value.den);
    Int exponent = numeratorExponent - denominatorExponent;
    F mantissa = numerator / denominator;
    if (squareRoot) {
        if (exponent % 2 != 0) {
            mantissa *= 10;
            --exponent;
        }
        mantissa = sqrt(mantissa);
        exponent /= 2;
    }
    while (mantissa < 1) {
        mantissa *= 10;
        --exponent;
    }
    while (mantissa >= 10) {
        mantissa /= 10;
        ++exponent;
    }
    constexpr Int minimum = cplib::NumericLimits<F>::min_exponent10;
    constexpr Int maximum = cplib::NumericLimits<F>::max_exponent10;
    if (exponent < minimum - cplib::NumericLimits<F>::max_digits10 - 2 || exponent > maximum)
        throw std::invalid_argument("厳密円の近似出力は指定浮動小数点型の範囲外です");
    // 非正規化数でも、10^exponentを先に0へ丸めない。
    F result = exponent < minimum
        ? (mantissa * pow(F(10), F(minimum))) * pow(F(10), F(exponent - minimum))
        : mantissa * pow(F(10), F(exponent));
    result *= F(circleSign(value.num));
    if (!isfinite(result) || result == 0)
        throw std::invalid_argument("厳密円の近似出力は非零の有限値で表現できません");
    return result;
}
template <class T> Point<double> circlePointApprox(const Point<CircleRational<T>> &point) {
    return {circleApprox(point.x), circleApprox(point.y)};
}
template <class T> auto circleOffsetApprox(const Point<CircleRational<T>> &v,
                                          const CircleRational<T> &squaredFactor) {
    return Point<double>{
        double(circleSign(v.x.num)) * circleApprox(v.x * v.x * squaredFactor, true),
        double(circleSign(v.y.num)) * circleApprox(v.y * v.y * squaredFactor, true)};
}
template <class F, class T> F circleFloatCoordinateAs(const T &value) {
    if constexpr (std::is_arithmetic_v<T> || cplib::isFloatingPoint<T> ||
                  std::is_same_v<T, __int128_t> || std::is_same_v<T, __uint128_t>) {
        return static_cast<F>(value);
    } else if constexpr (isCircleFraction<T>) {
        if (value.den == CircleBase<T>(0))
            throw std::invalid_argument("float円の分数座標は有限である必要があります");
        // 分子・分母を別々にfloatへ変換すると、有限な比でもinf/infになり得る。
        // 最小負整数の符号反転も避け、絶対値を十進表記から扱う。
        CircleRational<CircleBase<T>> raw;
        raw.num = value.num;
        raw.den = value.den;
        return circleApprox<F>(raw) * F(circleSign(value.den));
    } else {
        return circleApprox<F>(circleFraction(value));
    }
}
template <class T> double circleFloatCoordinate(const T &value) {
    return circleFloatCoordinateAs<double>(value);
}
template <class T> Point<double> circleFloatPoint(const Point<T> &point) {
    return {circleFloatCoordinate(point.x), circleFloatCoordinate(point.y)};
}
template <class F, class T> Point<F> circleFloatPointAs(const Point<T> &point) {
    return {circleFloatCoordinateAs<F>(point.x), circleFloatCoordinateAs<F>(point.y)};
}
template <class... T> using CircleFloatType = std::common_type_t<
    double, std::conditional_t<cplib::isFloatingPoint<T>, T, double>...>;
template <class T> Int circleRadicalSign(const CircleRational<T> &x, Int sign,
                                        const CircleRational<T> &radicand) {
    Int xs = circleSign(x.num);
    if (radicand.num == T(0))
        return xs;
    if (xs == 0 || xs == sign)
        return sign;
    Int comparison = circleSign((x * x - radicand).num);
    return comparison == 0 ? 0 : comparison > 0 ? xs : sign;
}
}

namespace cplib {
// 整数・分数の点を同じ基底型のFractionの点へ変換する。O(1)算術操作とgcdの費用。
// 整数幅: 整数入力は分母1にするだけ。B=M、M=10^9ならint可（自動昇格なし）。
template <class T> auto toExactPoint(const Point<T> &point) {
    return detail::circleStoredPoint(detail::circleExactPoint(point));
}
}
