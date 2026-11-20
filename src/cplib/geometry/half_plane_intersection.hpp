#pragma once
#include <cplib/geometry/convex_polygon.hpp>
#include <cplib/math/fractions.hpp>
#include <initializer_list>
#include <stdexcept>

namespace cplib {

// 左閉半平面の共通部分が空、有界、非有界のいずれかを表す。
enum HalfPlaneIntersectionKind { hpiEmpty, hpiBounded, hpiUnbounded };

namespace detail::hpi {

template <class T> struct IsFraction : std::false_type {};
template <class T> struct IsFraction<Fraction<T>> : std::true_type {};

// 浮動小数点では有限値、分数では正の分母を持つ値だけを許す。
template <class T> T checked(T value) {
    static_assert(cplib::isFloatingPoint<T> || IsFraction<T>::value,
                  "half_plane_intersection requires floating point or Fraction coordinates");
    if constexpr (cplib::isFloatingPoint<T>) {
        if (!isfinite(value))
            throw std::invalid_argument("半平面交差の値が有限値ではありません");
    } else if constexpr (IsFraction<T>::value) {
        if (value.den <= 0)
            throw std::invalid_argument("半平面交差には正の分母の有限分数が必要です");
    }
    return value;
}

template <class T> struct Bound {
    // -1 は負の無限大、0 は有限値、1 は正の無限大。
    int infinity = 0;
    T value{};
};

template <class T> Bound<T> finiteBound(T value) {
    return {0, checked(std::move(value))};
}

template <class T> int compare(const Bound<T> &a, const Bound<T> &b) {
    if (a.infinity != b.infinity)
        return a.infinity < b.infinity ? -1 : 1;
    if (a.infinity != 0)
        return 0;
    return a.value < b.value ? -1 : b.value < a.value ? 1 : 0;
}

template <class T> Bound<T> minimum(const Bound<T> &a, const Bound<T> &b) {
    return compare(a, b) <= 0 ? a : b;
}

template <class T> Bound<T> maximum(const Bound<T> &a, const Bound<T> &b) {
    return compare(a, b) >= 0 ? a : b;
}

template <class T> struct Function {
    T slope, intercept;
};

template <class T> struct Envelope {
    std::vector<Function<T>> functions;
    std::vector<Bound<T>> starts;
};

// 一次関数の最大値を表す包絡線を構築する。O(N log N) 演算、O(N) 空間。
template <class T> Envelope<T> makeEnvelope(std::vector<Function<T>> lines) {
    std::sort(lines.begin(), lines.end(), [](const Function<T> &a, const Function<T> &b) {
        if (a.slope != b.slope)
            return a.slope < b.slope;
        return a.intercept < b.intercept;
    });

    Envelope<T> result;
    for (const auto &line : lines) {
        // 同じ傾きなら、切片が最大の関数だけを残す。
        if (!result.functions.empty() && result.functions.back().slope == line.slope) {
            result.functions.pop_back();
            result.starts.pop_back();
        }
        Bound<T> start{-1, T{}};
        while (!result.functions.empty()) {
            const auto &last = result.functions.back();
            start = finiteBound((last.intercept - line.intercept) / (line.slope - last.slope));
            if (compare(start, result.starts.back()) > 0)
                break;
            result.functions.pop_back();
            result.starts.pop_back();
        }
        if (result.functions.empty())
            start = {-1, T{}};
        result.functions.push_back(line);
        result.starts.push_back(std::move(start));
    }
    return result;
}

// 包絡線上の値を二分探索で求める。O(log N) 演算。
template <class T> T valueAt(const Envelope<T> &envelope, const T &x) {
    std::size_t low = 0;
    std::size_t high = envelope.starts.size();
    const auto bound = finiteBound(x);
    while (high - low > 1) {
        const auto middle = low + (high - low) / 2;
        if (compare(envelope.starts[middle], bound) <= 0)
            low = middle;
        else
            high = middle;
    }
    const auto &line = envelope.functions[low];
    return checked(line.slope * x + line.intercept);
}

} // namespace detail::hpi

// vertices は有限の極点のみを持ち、constraints は入力の左閉半平面を保持する。
// dimension は空集合=-1、点=0、線分・半直線・直線=1、面=2。
template <class T> struct HalfPlaneIntersection {
    HalfPlaneIntersectionKind kind = hpiEmpty;
    Int dimension = -1;
    std::vector<Point<T>> vertices;
    std::vector<Line<T>> constraints;

    // 空集合・有界領域を独立した凸多角形に変換する。V を頂点数として時間・領域 O(V)。
    // 空集合は空の多角形にし、点・線分も許す。非有界領域は std::domain_error。
    // 半平面交差が返した正規化済みの頂点順を使い、座標演算は行わない。
    // 変換後の包含判定は ConvexPolygon の比較規則（浮動小数点では EPS 使用）に従う。
    ConvexPolygon<T> toConvexPolygon() const {
        if (kind == hpiUnbounded)
            throw std::domain_error("非有界な半平面交差は凸多角形に変換できません");
        if (kind == hpiEmpty)
            return {};
        return ConvexPolygon<T>(vertices.begin(), vertices.end());
    }

    // 点がすべての左閉半平面を満たすか判定する。O(N) 演算、追加要素数 O(1)。
    // N は制約数。非有限座標や分母が正でない分数は invalid_argument。
    bool contains(Point<T> point) const {
        using detail::hpi::checked;
        checked(point.x);
        checked(point.y);
        if (kind == hpiEmpty)
            return false;
        const T zero(0);
        for (const auto &line : constraints) {
            const T dx = checked(line.t.x - line.s.x);
            const T dy = checked(line.t.y - line.s.y);
            const T side = checked(dx * (point.y - line.s.y) - dy * (point.x - line.s.x));
            if (side < zero)
                return false;
        }
        return true;
    }
};

// 有向直線の左閉半平面の共通部分を求める。O(N log(N+1)) 演算、O(N) 空間。
// T は浮動小数点型または Fraction。計算量では座標型の四則演算・比較を1演算と数える。
// 分数の分母は正で、中間整数演算がオーバーフローしないこと。
// vertices は下側を x 昇順、上側を x 降順に並べ、有界な面なら反時計回り。
// 非有界領域の全体は constraints で表し、制約なしは全平面を返す。
// 入力は有限座標かつ異なる端点。非有限の入力・計算結果は invalid_argument。
// 浮動小数点の比較には EPS を使わず、GEOMETRY_EPS も参照しない。
// 丸めにより退化・平行・近接を誤分類し得るため、厳密な判定には Fraction を使う。
template <class T>
HalfPlaneIntersection<T> half_plane_intersection(std::span<const Line<T>> lines) {
    using namespace detail::hpi;
    static_assert(cplib::isFloatingPoint<T> || IsFraction<T>::value,
                  "half_plane_intersection requires floating point or Fraction coordinates");
    const T zero(0);
    const Bound<T> negativeInfinity{-1, T{}};
    const Bound<T> positiveInfinity{1, T{}};
    HalfPlaneIntersection<T> result;
    result.constraints.assign(lines.begin(), lines.end());

    std::vector<Function<T>> lowerLines, negativeUpperLines;
    auto xmin = negativeInfinity;
    auto xmax = positiveInfinity;
    for (const auto &line : lines) {
        checked(line.s.x);
        checked(line.s.y);
        checked(line.t.x);
        checked(line.t.y);
        const T dx = checked(line.t.x - line.s.x);
        const T dy = checked(line.t.y - line.s.y);
        if (dx == zero && dy == zero)
            throw std::invalid_argument("半平面の端点は異なる必要があります");
        if (dx == zero) {
            if (dy > zero)
                xmax = minimum(xmax, finiteBound(line.s.x));
            else
                xmin = maximum(xmin, finiteBound(line.s.x));
        } else {
            const T slope = checked(dy / dx);
            const T intercept = checked(checked(dx * line.s.y - dy * line.s.x) / dx);
            if (dx > zero)
                lowerLines.push_back({slope, intercept});
            else
                negativeUpperLines.push_back({-slope, -intercept});
        }
    }
    const auto lower = makeEnvelope(std::move(lowerLines));
    const auto upper = makeEnvelope(std::move(negativeUpperLines));
    if (compare(xmin, xmax) > 0)
        return result;

    auto left = xmin;
    auto right = xmax;
    bool hasArea = false;
    bool leftMeet = false;
    bool rightMeet = false;
    if (!lower.functions.empty() && !upper.functions.empty()) {
        // 下側と上側の包絡線を同時にたどり、下側 <= 上側となる x の区間を求める。
        left = positiveInfinity;
        right = negativeInfinity;
        std::size_t i = 0, j = 0;
        auto start = negativeInfinity;
        while (true) {
            const auto nextLower = i + 1 < lower.starts.size() ? lower.starts[i + 1] : positiveInfinity;
            const auto nextUpper = j + 1 < upper.starts.size() ? upper.starts[j + 1] : positiveInfinity;
            const auto finish = minimum(nextLower, nextUpper);
            auto a = maximum(start, xmin);
            auto b = minimum(finish, xmax);
            const T slope = checked(lower.functions[i].slope + upper.functions[j].slope);
            const T intercept = checked(lower.functions[i].intercept + upper.functions[j].intercept);
            auto root = positiveInfinity;
            if (slope != zero)
                root = finiteBound(-intercept / slope);
            if (slope > zero)
                b = minimum(b, root);
            else if (slope < zero)
                a = maximum(a, root);
            else if (intercept > zero) {
                a = positiveInfinity;
                b = negativeInfinity;
            }
            if (compare(a, b) <= 0) {
                const bool identical = slope == zero && intercept == zero;
                const bool aMeet = identical || (slope != zero && compare(a, root) == 0);
                const bool bMeet = identical || (slope != zero && compare(b, root) == 0);
                if (compare(a, left) < 0)
                    leftMeet = aMeet;
                else if (compare(a, left) == 0)
                    leftMeet = leftMeet || aMeet;
                if (compare(b, right) > 0)
                    rightMeet = bMeet;
                else if (compare(b, right) == 0)
                    rightMeet = rightMeet || bMeet;
                left = minimum(left, a);
                right = maximum(right, b);
                if (compare(a, b) < 0 && (slope != zero || intercept < zero))
                    hasArea = true;
            }
            if (finish.infinity == 1)
                break;
            if (compare(nextLower, finish) == 0)
                ++i;
            if (compare(nextUpper, finish) == 0)
                ++j;
            start = finish;
        }
        if (compare(left, right) > 0)
            return result;
    } else if (compare(left, right) < 0) {
        hasArea = true;
    }

    const bool bothSides = !lower.functions.empty() && !upper.functions.empty();
    result.kind = left.infinity == 0 && right.infinity == 0 && bothSides ? hpiBounded : hpiUnbounded;
    if (compare(left, right) == 0)
        result.dimension = bothSides && (leftMeet || rightMeet) ? 0 : 1;
    else
        result.dimension = hasArea ? 2 : 1;

    // 下側を左から右、上側を右から左へ出力し、つなぎ目の重複を除く。
    auto append = [&](Point<T> point) {
        if (result.vertices.empty() || !exact_equal(result.vertices.back(), point))
            result.vertices.push_back(std::move(point));
    };
    if (!lower.functions.empty()) {
        if (left.infinity == 0)
            append({left.value, valueAt(lower, left.value)});
        for (std::size_t k = 1; k < lower.starts.size(); ++k) {
            const auto &x = lower.starts[k];
            if (compare(left, x) < 0 && compare(x, right) < 0)
                append({x.value, valueAt(lower, x.value)});
        }
        if (right.infinity == 0 && compare(left, right) != 0)
            append({right.value, valueAt(lower, right.value)});
    }
    if (!upper.functions.empty()) {
        if (right.infinity == 0) {
            const T y = rightMeet ? valueAt(lower, right.value) : -valueAt(upper, right.value);
            append({right.value, y});
        }
        for (std::size_t k = upper.starts.size(); k-- > 1;) {
            const auto &x = upper.starts[k];
            if (compare(left, x) < 0 && compare(x, right) < 0)
                append({x.value, -valueAt(upper, x.value)});
        }
        if (left.infinity == 0 && compare(left, right) != 0) {
            const T y = leftMeet ? valueAt(lower, left.value) : -valueAt(upper, left.value);
            append({left.value, y});
        }
    }
    if (result.vertices.size() > 1 && exact_equal(result.vertices.front(), result.vertices.back()))
        result.vertices.pop_back();
    return result;
}

// vector の有向直線から左閉半平面の共通部分を求める。O(N log(N+1)) 演算、O(N) 空間。
template <class T, class Allocator>
HalfPlaneIntersection<T> half_plane_intersection(const std::vector<Line<T>, Allocator> &lines) {
    return half_plane_intersection<T>(std::span<const Line<T>>(lines.data(), lines.size()));
}

// array の有向直線から左閉半平面の共通部分を求める。O(N log(N+1)) 演算、O(N) 空間。
template <class T, std::size_t N>
HalfPlaneIntersection<T> half_plane_intersection(const std::array<Line<T>, N> &lines) {
    return half_plane_intersection<T>(std::span<const Line<T>>(lines.data(), lines.size()));
}

// 初期化リストの有向直線から左閉半平面の共通部分を求める。O(N log(N+1)) 演算、O(N) 空間。
template <class T>
HalfPlaneIntersection<T> half_plane_intersection(std::initializer_list<Line<T>> lines) {
    return half_plane_intersection<T>(std::span<const Line<T>>(lines.begin(), lines.size()));
}

// 点が保持された左閉半平面をすべて満たすか判定する。O(N) 演算、追加要素数 O(1)。
template <class T> bool contains(const HalfPlaneIntersection<T> &region, Point<T> point) {
    return region.contains(std::move(point));
}

// 空集合・有界領域を凸多角形に変換する。V を頂点数として時間・領域 O(V)。
// 空集合は空の多角形、非有界領域は std::domain_error。点・線分も変換できる。
template <class T> ConvexPolygon<T> toConvexPolygon(const HalfPlaneIntersection<T> &region) {
    return region.toConvexPolygon();
}

namespace detail::hpi {

// 符号付き整数型（__int128、BigIntを含む）。浮動小数点・Fractionは対象外。
template <class T> concept IntegerCoefficient =
    !std::is_unsigned_v<T> && !std::is_same_v<T, __uint128_t> &&
    requires(T a, T b) { T(a % b); };

// 正の分母を持つ分数を、交差積を作らず比較する。O(log M)回の除算。
// 符号反転が可能な範囲の整数を使う（符号付き最小値は対象外）。
template <IntegerCoefficient T>
int compareRatios(T a, T b, T c, T d) {
    if ((a < T(0)) != (c < T(0)))
        return a < T(0) ? -1 : 1;
    bool reversed = a < T(0);
    if (reversed) {
        a = -a;
        c = -c;
    }
    while (true) {
        const T q = a / b, r = c / d;
        if (q != r) {
            const int order = q < r ? -1 : 1;
            return reversed ? -order : order;
        }
        a %= b;
        c %= d;
        if (a == T(0) || c == T(0)) {
            const int order = a == c ? 0 : a == T(0) ? -1 : 1;
            return reversed ? -order : order;
        }
        std::swap(a, b);
        std::swap(c, d);
        reversed = !reversed;
    }
}

template <IntegerCoefficient T> struct IntegerBound {
    int infinity = 0;
    T num = T(0), den = T(1);
};

template <IntegerCoefficient T> IntegerBound<T> integerBound(T num, T den) {
    if (den < T(0)) {
        num = -num;
        den = -den;
    }
    return {0, std::move(num), std::move(den)};
}

template <IntegerCoefficient T>
int compare(const IntegerBound<T> &a, const IntegerBound<T> &b) {
    if (a.infinity != b.infinity)
        return a.infinity < b.infinity ? -1 : 1;
    if (a.infinity != 0)
        return 0;
    return compareRatios(a.num, a.den, b.num, b.den);
}

template <IntegerCoefficient T>
IntegerBound<T> minimum(const IntegerBound<T> &a, const IntegerBound<T> &b) {
    return compare(a, b) <= 0 ? a : b;
}

template <IntegerCoefficient T>
IntegerBound<T> maximum(const IntegerBound<T> &a, const IntegerBound<T> &b) {
    return compare(a, b) >= 0 ? a : b;
}

// (p*x+q)/d を表す。d>0で、係数は入力の大きさを保つ。
template <IntegerCoefficient T> struct IntegerFunction {
    T p, q, d;
};

template <IntegerCoefficient T> struct IntegerEnvelope {
    std::vector<IntegerFunction<T>> functions;
    std::vector<IntegerBound<T>> starts;
};

template <IntegerCoefficient T>
IntegerEnvelope<T> makeIntegerEnvelope(std::vector<IntegerFunction<T>> lines) {
    std::sort(lines.begin(), lines.end(), [](const auto &a, const auto &b) {
        const T left = a.p * b.d, right = b.p * a.d;
        if (left != right)
            return left < right;
        return a.q * b.d < b.q * a.d;
    });
    IntegerEnvelope<T> result;
    for (const auto &line : lines) {
        if (!result.functions.empty()) {
            const auto &last = result.functions.back();
            if (last.p * line.d == line.p * last.d) {
                result.functions.pop_back();
                result.starts.pop_back();
            }
        }
        IntegerBound<T> start{-1};
        while (!result.functions.empty()) {
            const auto &last = result.functions.back();
            start = integerBound<T>(last.q * line.d - line.q * last.d,
                                    line.p * last.d - last.p * line.d);
            if (compare(start, result.starts.back()) > 0)
                break;
            result.functions.pop_back();
            result.starts.pop_back();
        }
        if (result.functions.empty())
            start = {-1};
        result.functions.push_back(line);
        result.starts.push_back(std::move(start));
    }
    return result;
}

template <IntegerCoefficient T>
Fraction<T> integerValueAt(const IntegerFunction<T> &line, const IntegerBound<T> &x) {
    return Fraction<T>(line.p * x.num + line.q * x.den, line.d * x.den);
}

template <IntegerCoefficient T>
Fraction<T> integerValueAt(const IntegerEnvelope<T> &envelope, const IntegerBound<T> &x) {
    std::size_t low = 0, high = envelope.starts.size();
    while (high - low > 1) {
        const auto middle = low + (high - low) / 2;
        if (compare(envelope.starts[middle], x) <= 0)
            low = middle;
        else
            high = middle;
    }
    return integerValueAt(envelope.functions[low], x);
}

} // namespace detail::hpi

// 整数係数 ax+by+c<=0 の共通部分。頂点は約分済みの有理数で保持する。
// kind、dimension、頂点順は HalfPlaneIntersection と同じ。constraints は入力係数。
template <detail::hpi::IntegerCoefficient T> struct IntegerHalfPlaneIntersection {
    HalfPlaneIntersectionKind kind = hpiEmpty;
    Int dimension = -1;
    std::vector<Point<Fraction<T>>> vertices;
    std::vector<std::array<T, 3>> constraints;

    // 座標演算をせず凸多角形に変換する。時間・領域O(V)。非有界ならdomain_error。
    // 整数幅: 格納済み分子・分母のコピーのみ。追加幅不要。
    ConvexPolygon<Fraction<T>> toConvexPolygon() const {
        if (kind == hpiUnbounded)
            throw std::domain_error("非有界な半平面交差は凸多角形に変換できません");
        if (kind == hpiEmpty)
            return {};
        return ConvexPolygon<Fraction<T>>(vertices.begin(), vertices.end());
    }

    // 整数点の包含。O(N)演算。係数と座標の積・和がTに収まること。
    // 整数幅: A=max(1,|a|,|b|),C=max(1,|c|),Q=max(1,|x|,|y|)としてB=2*A*Q+C。
    // A=C=Q=10^9なら62ビット、long longで十分。領域構築時の幅とは別条件。
    bool contains(const Point<T> &point) const {
        if (kind == hpiEmpty)
            return false;
        for (const auto &[a, b, c] : constraints)
            if (a * point.x + b * point.y + c > T(0))
                return false;
        return true;
    }

    // 有理点の包含。O(N log M)演算。分母は正で、下記の積・和がTに収まること。
    // 分母同士は掛けず、a*x.num と b*y.num+c*y.den を作って分数比較する。
    // 整数幅: 係数上限A,C、分子絶対値上限P>=1、正の分母上限D>=1としてB=A*P+C*D。
    // A=C=P=D=10^9なら62ビット、long longで十分。約分後だけでなく入力の分子・分母で見積もる。
    bool contains(const Point<Fraction<T>> &point) const {
        if (point.x.den <= T(0) || point.y.den <= T(0))
            throw std::invalid_argument("包含判定には正の分母の有限分数が必要です");
        if (kind == hpiEmpty)
            return false;
        for (const auto &[a, b, c] : constraints)
            if (detail::hpi::compareRatios<T>(a * point.x.num, point.x.den,
                    -(b * point.y.num + c * point.y.den), point.y.den) > 0)
                return false;
        return true;
    }
};

// 整数係数 {a,b,c} の閉半平面 ax+by+c<=0 の共通部分を厳密に求める。
// Tは符号付き整数（__int128、BigInt可）。自動昇格・浮動小数点演算は行わない。
// a=b=0はc<=0なら制約なし、c>0なら空集合。無入力は全平面。
// A=max(1,|a|,|b|), C=max(1,|c|)とすると、中間値は絶対値
// max(4*A*A*C, 2*A*A*A)以下。これと符号反転がTに収まること。
// 分数比較は交差積を作らず除算で行い、汎用Fractionの四則演算は使わない。
// 中間値上限をMとしてO(N log(N+1)+N log(M+1))整数演算、O(N)空間。
// 返した頂点に対する後続の面積・凸多角形演算のオーバーフロー条件は別途必要。
template <detail::hpi::IntegerCoefficient T>
// 整数幅: B=max(4*A^2*C,2*A^3)としてW(B)ビット（base.hppの定義）。
// 全係数の絶対値<=10^9ならB<=4*10^27、93ビット、T=__int128で十分（long long不可）。
// 座標上限MのLineをtoCoefficientsした場合はA<=2M,C<=4M^2なのでB<=64M^4。
// M=10^9なら127ビット。係数を作る型と、この関数へ渡す型の両方に注意。
IntegerHalfPlaneIntersection<T> half_plane_intersection(std::span<const std::array<T, 3>> lines) {
    using namespace detail::hpi;
    using Bound = IntegerBound<T>;
    const T zero(0);
    const Bound negativeInfinity{-1}, positiveInfinity{1};
    IntegerHalfPlaneIntersection<T> result;
    result.constraints.assign(lines.begin(), lines.end());
    std::vector<IntegerFunction<T>> lowerLines, negativeUpperLines;
    auto xmin = negativeInfinity, xmax = positiveInfinity;
    for (const auto &[a, b, c] : lines) {
        if (b == zero) {
            if (a > zero)
                xmax = minimum(xmax, integerBound<T>(-c, a));
            else if (a < zero)
                xmin = maximum(xmin, integerBound<T>(-c, a));
            else if (c > zero)
                return result;
        } else if (b < zero)
            lowerLines.push_back({a, c, -b});
        else
            negativeUpperLines.push_back({a, c, b});
    }
    if (compare(xmin, xmax) > 0)
        return result;
    const auto lower = makeIntegerEnvelope(std::move(lowerLines));
    const auto upper = makeIntegerEnvelope(std::move(negativeUpperLines));
    auto left = xmin, right = xmax;
    bool hasArea = false, leftMeet = false, rightMeet = false;
    if (!lower.functions.empty() && !upper.functions.empty()) {
        left = positiveInfinity;
        right = negativeInfinity;
        std::size_t i = 0, j = 0;
        auto start = negativeInfinity;
        while (true) {
            const auto nextLower = i + 1 < lower.starts.size() ? lower.starts[i + 1] : positiveInfinity;
            const auto nextUpper = j + 1 < upper.starts.size() ? upper.starts[j + 1] : positiveInfinity;
            const auto finish = minimum(nextLower, nextUpper);
            auto a = maximum(start, xmin), b = minimum(finish, xmax);
            const auto &lo = lower.functions[i];
            const auto &hi = upper.functions[j];
            const T slope = lo.p * hi.d + hi.p * lo.d;
            const T intercept = lo.q * hi.d + hi.q * lo.d;
            auto root = positiveInfinity;
            if (slope != zero)
                root = integerBound<T>(-intercept, slope);
            if (slope > zero)
                b = minimum(b, root);
            else if (slope < zero)
                a = maximum(a, root);
            else if (intercept > zero) {
                a = positiveInfinity;
                b = negativeInfinity;
            }
            if (compare(a, b) <= 0) {
                const bool identical = slope == zero && intercept == zero;
                const bool aMeet = identical || (slope != zero && compare(a, root) == 0);
                const bool bMeet = identical || (slope != zero && compare(b, root) == 0);
                const int leftOrder = compare(a, left), rightOrder = compare(b, right);
                if (leftOrder < 0)
                    leftMeet = aMeet;
                else if (leftOrder == 0)
                    leftMeet = leftMeet || aMeet;
                if (rightOrder > 0)
                    rightMeet = bMeet;
                else if (rightOrder == 0)
                    rightMeet = rightMeet || bMeet;
                left = minimum(left, a);
                right = maximum(right, b);
                if (compare(a, b) < 0 && (slope != zero || intercept < zero))
                    hasArea = true;
            }
            if (finish.infinity == 1)
                break;
            if (compare(nextLower, finish) == 0)
                ++i;
            if (compare(nextUpper, finish) == 0)
                ++j;
            start = finish;
        }
        if (compare(left, right) > 0)
            return result;
    } else if (compare(left, right) < 0)
        hasArea = true;

    const bool bothSides = !lower.functions.empty() && !upper.functions.empty();
    result.kind = left.infinity == 0 && right.infinity == 0 && bothSides ? hpiBounded : hpiUnbounded;
    if (compare(left, right) == 0)
        result.dimension = bothSides && (leftMeet || rightMeet) ? 0 : 1;
    else
        result.dimension = hasArea ? 2 : 1;

    // 全座標をここで約分する。同値比較にも分子・分母の交差積を使わない。
    auto samePoint = [](const auto &a, const auto &b) {
        return a.x.num == b.x.num && a.x.den == b.x.den &&
               a.y.num == b.y.num && a.y.den == b.y.den;
    };
    auto append = [&](const Bound &x, Fraction<T> y) {
        Point<Fraction<T>> point{Fraction<T>(x.num, x.den), std::move(y)};
        if (result.vertices.empty() || !samePoint(result.vertices.back(), point))
            result.vertices.push_back(std::move(point));
    };
    if (!lower.functions.empty()) {
        if (left.infinity == 0)
            append(left, integerValueAt(lower, left));
        for (std::size_t k = 1; k < lower.starts.size(); ++k) {
            const auto &x = lower.starts[k];
            if (compare(left, x) < 0 && compare(x, right) < 0)
                append(x, integerValueAt(lower.functions[k], x));
        }
        if (right.infinity == 0 && compare(left, right) != 0)
            append(right, integerValueAt(lower, right));
    }
    if (!upper.functions.empty()) {
        if (right.infinity == 0)
            append(right, rightMeet ? integerValueAt(lower, right) : -integerValueAt(upper, right));
        for (std::size_t k = upper.starts.size(); k-- > 1;) {
            const auto &x = upper.starts[k];
            if (compare(left, x) < 0 && compare(x, right) < 0)
                append(x, -integerValueAt(upper.functions[k], x));
        }
        if (left.infinity == 0 && compare(left, right) != 0)
            append(left, leftMeet ? integerValueAt(lower, left) : -integerValueAt(upper, left));
    }
    if (result.vertices.size() > 1 && samePoint(result.vertices.front(), result.vertices.back()))
        result.vertices.pop_back();
    return result;
}

// 整数係数版のvector、array、初期化リスト入口。計算量・前提はspan版と同じ。
template <detail::hpi::IntegerCoefficient T, class Allocator>
// 整数幅: span版と同じB=max(4*A^2*C,2*A^3)。全係数10^9以下なら93ビット、__int128。
IntegerHalfPlaneIntersection<T> half_plane_intersection(const std::vector<std::array<T, 3>, Allocator> &lines) {
    return half_plane_intersection<T>(std::span<const std::array<T, 3>>(lines.data(), lines.size()));
}

template <detail::hpi::IntegerCoefficient T, std::size_t N>
// 整数幅: span版と同じB=max(4*A^2*C,2*A^3)。全係数10^9以下なら93ビット、__int128。
IntegerHalfPlaneIntersection<T> half_plane_intersection(const std::array<std::array<T, 3>, N> &lines) {
    return half_plane_intersection<T>(std::span<const std::array<T, 3>>(lines.data(), lines.size()));
}

template <detail::hpi::IntegerCoefficient T>
// 整数幅: span版と同じB=max(4*A^2*C,2*A^3)。全係数10^9以下なら93ビット、__int128。
IntegerHalfPlaneIntersection<T> half_plane_intersection(std::initializer_list<std::array<T, 3>> lines) {
    return half_plane_intersection<T>(std::span<const std::array<T, 3>>(lines.begin(), lines.size()));
}

// 包含・凸多角形変換の自由関数。計算量・前提は対応するメンバー関数と同じ。
template <detail::hpi::IntegerCoefficient T>
// 整数幅: メンバー版と同じB=2*A*Q+C。A=C=Q=10^9なら62ビット、long long。
bool contains(const IntegerHalfPlaneIntersection<T> &region, const Point<T> &point) {
    return region.contains(point);
}

template <detail::hpi::IntegerCoefficient T>
// 整数幅: メンバー版と同じB=A*P+C*D。全上限10^9なら62ビット、long long。
bool contains(const IntegerHalfPlaneIntersection<T> &region, const Point<Fraction<T>> &point) {
    return region.contains(point);
}

template <detail::hpi::IntegerCoefficient T>
// 整数幅: 格納済み分子・分母のコピーのみ。追加幅不要。
ConvexPolygon<Fraction<T>> toConvexPolygon(const IntegerHalfPlaneIntersection<T> &region) {
    return region.toConvexPolygon();
}

} // namespace cplib
