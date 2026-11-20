#pragma once
#include <cplib/geometry/circle_integer_impl.hpp>

namespace cplib {
// 円の計算量は基底型の算術・gcdを一操作として数える。多倍長演算と十進変換の費用は桁数に依存する。
// gcdは中間値の絶対値上限Vに対して O(log(V+1)) 回の除算を要する。
// 整数・Fractionは型を自動昇格せず、全中間値が指定基底型に収まることを前提とする。
// 厳密円の整数幅:
// Mは「保持する3点points()」の座標上限。元の構築引数の上限とは異なる。
// initCircle(中心,半径)は保持座標<=2M、initCircle(中心,通過点)は<=3M、
// initDiameterCircleは<=4M・尺度2になる（この行のMは構築引数の座標・半径上限）。
// 以下のhは、全入力とcircleExactData内の中間整数、および算出した中心・半径二乗の
// 分子の絶対値と正の分母がすべて2^h未満となる幅（符号を除く）。hは2以上とする。
// 整数の保持座標・問い合わせ座標がM以下、尺度の絶対値が1または2なら、
// k=ceil(log2(M+1)), h=6k+11で十分。例えばM=10^9ではh=191。
// 一般の整数尺度S>=1ではs=ceil(log2(S+1))として
// h=max(6k+11,4k+2s+8,2k+s+4,2s)で十分。問い合わせ入力もこの幅に収める。
// 根拠: 差<=2M、中心オフセットの分子<=32M^3・分母<=16M^2。
// その二乗和の分子<=2048M^6・分母<=256M^4で、尺度による分母のS倍/S^2倍も考慮する。
// 同一行列式に由来する分母を毎回約分するため、二乗和もこの上限で覆える。
// h式は十分条件であり最小幅ではない。共通分母・約分を使えばさらに小さくできる場合がある。
// Fraction入力にはMからhへの式を使わず、circle_integer_impl.hppの演算別規則でhを求める。
// 以下で「hに対し...ビット」と書く幅は符号込み。M=10^9の一般の3点円の
// 中心・交点・接線計算は、この保証では__int128でも足りずBigInt等の多倍長整数を使う。
// 近似値を返す関数も、近似へ変換する前の厳密計算には同じ幅が必要。自動昇格はしない。
enum CirclePointLocation { circleOutside, circleBoundary, circleInside };
enum CircleResultKind { circleFinite, circleInfinite };
template <class F> struct CircleIntersectionsT {
    CircleResultKind kind = circleFinite;
    std::vector<Point<F>> points;
};
template <class F> struct CircleTangentT {
    // 第1引数の円上の接点。
    Point<F> first;
    // common_tangentsでは第2引数の円上の接点、tangent_linesでは渡した点p。
    Point<F> second;
    // 接線の単位方向ベクトル。接点が一致する場合にも使える。
    Point<F> direction;
};
template <class F> struct CircleTangentsT {
    CircleResultKind kind = circleFinite;
    std::vector<CircleTangentT<F>> tangents;
};
using CircleIntersections = CircleIntersectionsT<double>;
using CircleTangent = CircleTangentT<double>;
using CircleTangents = CircleTangentsT<double>;
template <class T, class... S> using CircleResultScalar = std::common_type_t<
    std::conditional_t<cplib::isFloatingPoint<T>, T, double>,
    std::conditional_t<cplib::isFloatingPoint<S>, S, float>...>;
template <class T = double> class Circle;

namespace detail {
template <class F> inline void circleCheckPoint(Point<F> p) {
    if (!isfinite(p.x) || !isfinite(p.y))
        throw std::invalid_argument("円の座標と中間値は有限である必要があります");
}
inline void circleCheckTolerance(long double tolerance) {
    if (!isfinite(tolerance) || tolerance < 0 || tolerance >= 1)
        throw std::invalid_argument("相対許容誤差は0以上1未満である必要があります");
}
template <class F> inline F circleLength(Point<F> p) {
    circleCheckPoint(p);
    F result = hypot(p.x, p.y);
    if (!isfinite(result))
        throw std::invalid_argument("円の距離は使用する浮動小数点型で表現可能である必要があります");
    return result;
}
template <class F> inline Point<F> circlePerpendicular(Point<F> p) {
    return {-p.y, p.x};
}
template <class C> auto circleExactData(const C &circle) {
    auto stored = circle.points();
    auto a = circleExactPoint(stored[0]), b = circleExactPoint(stored[1]);
    auto c = circleExactPoint(stored[2]);
    auto scale = circleFraction(circle.point_scale());
    using R = decltype(scale);
    if (a == b && a == c)
        return std::pair{a * (R(1) / scale), R(0)};
    auto u = b - a, v = c - a;
    auto determinant = cplib::cross(u, v) * 2;
    Point<R> offset{
        (v.y * cplib::norm(u) - u.y * cplib::norm(v)) / determinant,
        (u.x * cplib::norm(v) - v.x * cplib::norm(u)) / determinant};
    return std::pair{(a + offset) * (R(1) / scale), cplib::norm(offset) / (scale * scale)};
}
template <class Derived> class CircleOperations {
public:
    // 点の円盤内外・境界を分類する。厳密版はO(1)算術操作、float版は O(1)。
    // 整数幅: 同名の自由関数classifyの条件と同じ（各オーバーロードのコメントを参照）。
    template <class S> CirclePointLocation classify(Point<S> p, long double tolerance = 1e-10L) const;
    // 境界を含む円盤の点包含を判定する。厳密版はO(1)算術操作、float版は O(1)。
    // 整数幅: 同名の自由関数containsの条件と同じ（各オーバーロードのコメントを参照）。
    template <class S> bool contains(Point<S> p, long double tolerance = 1e-10L) const;
    // 円周上かを判定する。厳密版はO(1)算術操作、float版は O(1)。
    // 整数幅: 同名の自由関数on_circleの条件と同じ（各オーバーロードのコメントを参照）。
    template <class S> bool on_circle(Point<S> p, long double tolerance = 1e-10L) const;
    // 円・直線・閉線分との交点を返す。O(1)算術操作と近似変換の費用。
    // 整数幅: 同名の自由関数cross_pointsの条件と同じ（各オーバーロードのコメントを参照）。
    template <class Object> auto cross_points(const Object &object, long double tolerance = 1e-10L) const;
    // 交点数を返す。無限個は-1。厳密版はO(1)算術操作、float版は O(1)。
    // 整数幅: 同名の自由関数intersection_countの条件と同じ（各オーバーロードのコメントを参照）。
    template <class Object> Int intersection_count(const Object &object, long double tolerance = 1e-10L) const;
    // 点からの接線を返す。O(1)算術操作と近似変換の費用。
    // 整数幅: 同名の自由関数tangent_linesの条件と同じ（各オーバーロードのコメントを参照）。
    template <class S> auto tangent_lines(Point<S> p, long double tolerance = 1e-10L) const;
    // 点からの接線数を返す。無限個は-1。O(1)算術操作。
    // 整数幅: 同名の自由関数tangent_countの条件と同じ（各オーバーロードのコメントを参照）。
    template <class S> Int tangent_count(Point<S> p, long double tolerance = 1e-10L) const;
    // 交点をfloat64で近似して返す。O(1)算術操作と変換費用。
    // 整数幅: 同名の自由関数cross_points_approxの条件と同じ（各オーバーロードのコメントを参照）。
    template <class Object> CircleIntersections cross_points_approx(const Object &object) const;
    // 点からの接線をfloat64で近似して返す。O(1)算術操作と変換費用。
    // 整数幅: 同名の自由関数tangent_lines_approxの条件と同じ（各オーバーロードのコメントを参照）。
    template <class S> CircleTangents tangent_lines_approx(Point<S> p) const;
    // 共通接線を返す。O(1)算術操作と近似変換の費用。
    // 整数幅: 同名の自由関数common_tangentsの条件と同じ（各オーバーロードのコメントを参照）。
    auto common_tangents(const Derived &other, long double tolerance = 1e-10L) const;
    // 共通接線数を返す。無限個は-1。O(1)算術操作。
    // 整数幅: 同名の自由関数common_tangent_countの条件と同じ（各オーバーロードのコメントを参照）。
    Int common_tangent_count(const Derived &other, long double tolerance = 1e-10L) const;
};
}

template <class T> class Circle : public detail::CircleOperations<Circle<T>> {
    static_assert(!cplib::isFloatingPoint<T>);
    using B = detail::CircleBase<T>;
    static_assert(!std::is_unsigned_v<B> && !std::is_same_v<B, __uint128_t>);
    std::array<Point<T>, 3> points_{};
    T scale_ = T(1);
    bool initialized_ = false;
public:
    // 整数幅: 0と1の初期化のみ。追加幅不要。
    // 未初期化の厳密円を作る。使用前に構築済みの円を代入する。O(1)算術操作。
    Circle() = default;
    // 共通尺度付きの3点を保持する。全点一致以外の共線や尺度0はinvalid_argument。O(1)算術操作。
    // 整数幅: 保持座標上限M、尺度絶対値SとしてB=max(8M^2,S)。M=10^9,S<=2なら64ビット、long long。
    // これは構築時の共線判定だけの幅。中心取得や包含判定では別の幅が必要。
    explicit Circle(std::array<Point<T>, 3> points, T scale = T(1))
        : points_(std::move(points)), scale_(std::move(scale)), initialized_(true) {
        if (detail::circleFraction(scale_).num == B(0))
            throw std::invalid_argument("円の尺度は非零である必要があります");
        auto a = detail::circleExactPoint(points_[0]), b = detail::circleExactPoint(points_[1]);
        auto c = detail::circleExactPoint(points_[2]);
        if (cplib::cross(b - a, c - a).num == B(0) && !(a == b && a == c))
            throw std::invalid_argument("円の3点は非共線または全点一致である必要があります");
    }
    // ID付きの保持点から座標を取り出して構築する。尺度と計算量は上のコンストラクタと同じ。
    template <class Id> requires (!std::is_void_v<Id>)
    explicit Circle(const std::array<Point<T, Id>, 3> &points, T scale = T(1))
        : Circle(std::array<Point<T>, 3>{points[0], points[1], points[2]}, std::move(scale)) {}
    // 保持する3点の分子座標を返す。実際の座標はpoint_scaleで割る。O(1)回の値コピー。
    // 整数幅: コピーのみ。保持座標がTに収まればよい。
    auto points() const {
        if (!initialized_)
            throw std::invalid_argument("厳密円は初期化されていません");
        return points_;
    }
    // 3点の共通尺度を返す。通常構築は1、直径円は2。O(1)回の値コピー。
    // 整数幅: コピーのみ。保持尺度がTに収まればよい。
    T point_scale() const {
        if (!initialized_)
            throw std::invalid_argument("厳密円は初期化されていません");
        return scale_;
    }
    // 共通尺度で割った実際の3点をFractionで返す。O(1)算術操作。
    // 整数幅: 保持座標上限M、尺度絶対値S>=1としてB=max(M,S)。M=10^9,S<=2ならint可。
    auto points_exact() const {
        auto stored = points();
        auto scale = detail::circleFraction(scale_);
        using R = decltype(scale);
        std::array<Point<Fraction<B>>, 3> result;
        for (std::size_t i = 0; i < 3; ++i)
            result[i] = detail::circleStoredPoint(detail::circleExactPoint(stored[i]) * (R(1) / scale));
        return result;
    }
    // 指定基底型のFractionで中心を厳密に返す。O(1)算術操作。
    // 整数幅: 冒頭のhに対しh+1ビット。保持座標M=10^9,尺度<=2なら192ビットで十分、BigInt等を使う。
    // 中心だけの取得でも内部で半径二乗まで計算する。
    auto center_exact() const {
        return detail::circleStoredPoint(detail::circleExactData(*this).first);
    }
    // 指定基底型のFractionで半径二乗を厳密に返す。O(1)算術操作。
    // 整数幅: 冒頭のhに対しh+1ビット。保持座標M=10^9,尺度<=2なら192ビットで十分、BigInt等を使う。
    // 中心だけの取得でも内部で半径二乗まで計算する。
    auto radius_squared_exact() const {
        return detail::circleStored(detail::circleExactData(*this).second);
    }
    // 中心を厳密に返す。O(1)算術操作。
    // 整数幅: 冒頭のhに対しh+1ビット。保持座標M=10^9,尺度<=2なら192ビットで十分、BigInt等を使う。
    // 中心だけの取得でも内部で半径二乗まで計算する。
    auto center() const {
        return center_exact();
    }
    // 中心をfloat64へ近似する。O(1)算術操作と十進変換の費用。
    // 整数幅: 近似変換前にh+1ビット。保持座標M=10^9,尺度<=2なら192ビットで十分、BigInt等。
    Point<double> center_approx() const {
        return detail::circlePointApprox(detail::circleExactData(*this).first);
    }
    // 半径をfloat64へ近似する。O(1)算術操作と十進変換の費用。
    // 整数幅: 近似変換前にh+1ビット。保持座標M=10^9,尺度<=2なら192ビットで十分、BigInt等。
    double radius_approx() const {
        return detail::circleApprox(detail::circleExactData(*this).second, true);
    }
    // 半径をfloat64へ近似する。O(1)算術操作と十進変換の費用。
    // 整数幅: 近似変換前にh+1ビット。保持座標M=10^9,尺度<=2なら192ビットで十分、BigInt等。
    double radius() const {
        return radius_approx();
    }
};
template <class T> requires cplib::isFloatingPoint<T>
class Circle<T> : public detail::CircleOperations<Circle<T>> {
    Point<T> center_{};
    T radius_ = 0;
public:
    // 原点を中心とする点円を作る。計算量 O(1)。
    Circle() = default;
    // 有限中心と非負の有限半径で構築する。不正値はinvalid_argument。計算量 O(1)。
    Circle(Point<T> center, T radius) : center_(center), radius_(radius) {
        detail::circleCheckPoint(center);
        if (!isfinite(radius) || radius < 0)
            throw std::invalid_argument("円の半径は非負かつ有限である必要があります");
    }
    // 中心を返す。計算量 O(1)。
    Point<T> center() const {
        return center_;
    }
    // 半径を返す。計算量 O(1)。
    T radius() const {
        return radius_;
    }
    // 中心をfloat64で返す。計算量 O(1)。
    Point<double> center_approx() const {
        return {double(center_.x), double(center_.y)};
    }
    // 半径をfloat64で返す。計算量 O(1)。
    double radius_approx() const {
        return radius_;
    }
};

// 同じ座標型の円の一致を厳密に判定する。O(1)算術操作。
// 整数幅: 中心と半径二乗を計算後、約分済み分子・分母を比較する。h+1ビットで十分。
template <class T> bool operator==(const Circle<T> &a, const Circle<T> &b) {
    if constexpr (cplib::isFloatingPoint<T>)
        return exact_equal(a.center(), b.center()) && a.radius() == b.radius();
    else {
        auto x = detail::circleExactData(a), y = detail::circleExactData(b);
        return x.first == y.first && x.second == y.second;
    }
}
// 中心と半径で構築する。floatがあればdouble以上、long doubleがあればlong doubleを使う。O(1)算術操作と変換費用。
// 整数幅: M=max(1,|中心座標|,半径)としてB=2M^2（保持座標<=2M、共線判定の外積は2r^2）。
// M=10^9なら62ビット、long longで構築可。その後の操作は保持3点から別途見積もる。
template <class T, class R, class Id = void> auto initCircle(const Point<T, Id> &p, const R &radius) {
    if constexpr (cplib::isFloatingPoint<T> || cplib::isFloatingPoint<R>) {
        using F = detail::CircleFloatType<T, R>;
        return Circle<F>(detail::circleFloatPointAs<F>(p), detail::circleFloatCoordinateAs<F>(radius));
    } else {
        using B = detail::CircleBase<T>;
        auto center = detail::circleExactPoint(p);
        auto supplied = detail::circleFraction(radius);
        detail::CircleRational<B> r(B(supplied.num), B(supplied.den));
        if (r.num < B(0))
            throw std::invalid_argument("円の半径は非負である必要があります");
        Point<T> a{detail::circleStoredLike<T>(center.x + r), p.y};
        Point<T> b{p.x, detail::circleStoredLike<T>(center.y + r)};
        Point<T> c{detail::circleStoredLike<T>(center.x - r), p.y};
        return Circle<T>({a, b, c});
    }
}
// x, y, rの順で中心と半径を読み込む。構築の計算量はinitCircleと同じで、入力費用が加わる。
// 整数幅: 中心・半径版initCircleと同じB=2M^2。M=10^9ならlong longで構築可。
template <class T> std::istream &operator>>(std::istream &in, Circle<T> &circle) {
    Point<T> center;
    T radius;
    if (in >> center >> radius) {
        if constexpr (cplib::isFloatingPoint<T>)
            circle = Circle<T>(center, radius);
        else
            circle = initCircle(center, radius);
    }
    return in;
}
// 中心と通過点で構築する。非floatは同じ座標型を使う。O(1)算術操作と変換費用。
// 整数幅: 構築引数の上限Mに対しB=16M^2。保持座標<=3M、共線判定は2*norm(q-p)。
// M=10^9なら65ビット、T=__int128で構築可。後続の操作の幅は別条件。
template <class T, class S, class FirstId = void, class SecondId = void>
auto initCircle(const Point<T, FirstId> &p, const Point<S, SecondId> &q) {
    if constexpr (cplib::isFloatingPoint<T> || cplib::isFloatingPoint<S>) {
        using F = detail::CircleFloatType<T, S>;
        auto center = detail::circleFloatPointAs<F>(p);
        return Circle<F>(center, detail::circleLength(detail::circleFloatPointAs<F>(q) - center));
    } else {
        static_assert(std::is_same_v<T, S>);
        auto center = detail::circleExactPoint(p), through = detail::circleExactPoint(q);
        auto v = through - center;
        auto b = center + Point{-v.y, v.x}, c = center - v;
        return Circle<T>({q, Point<T>{detail::circleStoredLike<T>(b.x), detail::circleStoredLike<T>(b.y)},
                         Point<T>{detail::circleStoredLike<T>(c.x), detail::circleStoredLike<T>(c.y)}});
    }
}
// 3点から円を作る。非floatは元の型と3点を保持する。O(1)算術操作と変換費用。
template <class T, class S, class U>
// 整数幅: B=8M^2（共線判定）。M=10^9なら64ビット、long longで構築可。後続の幅は別条件。
auto initCircle(const Point<T> &a, const Point<S> &b, const Point<U> &c) {
    if constexpr (cplib::isFloatingPoint<T> || cplib::isFloatingPoint<S> || cplib::isFloatingPoint<U>) {
        using F = detail::CircleFloatType<T, S, U>;
        auto p = detail::circleFloatPointAs<long double>(a), q = detail::circleFloatPointAs<long double>(b);
        auto r = detail::circleFloatPointAs<long double>(c);
        detail::circleCheckPoint(p);
        detail::circleCheckPoint(q);
        detail::circleCheckPoint(r);
        auto u = q - p, v = r - p;
        long double scale = std::max({abs(u.x), abs(u.y), abs(v.x), abs(v.y)});
        if (scale == 0)
            return Circle<F>(Point<F>{F(p.x), F(p.y)}, 0);
        auto un = u / scale, vn = v / scale;
        long double rawDeterminant = cplib::cross(u, v);
        long double determinant = isfinite(rawDeterminant) && rawDeterminant != 0
            ? (rawDeterminant / scale) / scale * 2 : 2 * cplib::cross(un, vn);
        if (determinant == 0)
            throw std::invalid_argument("円の3点は非共線または全点一致である必要があります");
        Point<long double> offset{(vn.y * cplib::norm(un) - un.y * cplib::norm(vn)) / determinant * scale,
                                  (un.x * cplib::norm(vn) - vn.x * cplib::norm(un)) / determinant * scale};
        auto center = p + offset;
        auto distance = [&](Point<long double> point) { return detail::circleLength(point - center); };
        return Circle<F>(Point<F>{F(center.x), F(center.y)},
                         F(std::max({distance(p), distance(q), distance(r)})));
    } else {
        static_assert(std::is_same_v<T, S> && std::is_same_v<T, U>);
        return Circle<T>({a, b, c});
    }
}
// 2点を直径の両端とする円を作る。非floatは尺度2の3点で保持する。O(1)算術操作と変換費用。
// 整数幅: 元の2点の座標上限Mに対しB=16M^2。保持座標<=4M、尺度2、外積は2*norm(a-b)。
// M=10^9なら65ビット、T=__int128で構築可。後続の幅は保持座標上限から見積もる。
template <class T> auto initDiameterCircle(const Point<T> &a, const Point<T> &b) {
    if constexpr (cplib::isFloatingPoint<T>) {
        using F = detail::CircleFloatType<T>;
        auto p = detail::circleFloatPointAs<F>(a), q = detail::circleFloatPointAs<F>(b);
        auto center = p + (q - p) * F(0.5);
        return Circle<F>(center, std::max(detail::circleLength(p - center), detail::circleLength(q - center)));
    } else {
        auto p = detail::circleExactPoint(a), q = detail::circleExactPoint(b);
        auto sum = p + q, v = p - q;
        auto middle = sum + Point{-v.y, v.x};
        auto first = p * 2, last = q * 2;
        return Circle<T>({Point<T>{detail::circleStoredLike<T>(first.x), detail::circleStoredLike<T>(first.y)},
                          Point<T>{detail::circleStoredLike<T>(middle.x), detail::circleStoredLike<T>(middle.y)},
                          Point<T>{detail::circleStoredLike<T>(last.x), detail::circleStoredLike<T>(last.y)}}, T(2));
    }
}
// 円をfloat64へ明示的に近似変換する。O(1)算術操作と十進変換の費用。
// 整数幅: 厳密円ではh+1ビット。保持座標M=10^9,尺度<=2なら192ビットで十分、BigInt等。
template <class T> auto toFloatCircle(const Circle<T> &circle) {
    return Circle<double>(circle.center_approx(), circle.radius_approx());
}
// 中心を返す。非floatはFractionの点となる。O(1)算術操作。
// 整数幅: 厳密円ではh+1ビット。保持座標M=10^9,尺度<=2なら192ビットで十分、BigInt等。
template <class T> auto center(const Circle<T> &circle) {
    return circle.center();
}
// 半径を返す。非floatはfloat64に近似する。O(1)算術操作と変換費用。
// 整数幅: 厳密円ではh+1ビット。保持座標M=10^9,尺度<=2なら192ビットで十分、BigInt等。
template <class T> auto radius(const Circle<T> &circle) {
    return circle.radius();
}
// 中心を厳密に返す。O(1)算術操作。
// 整数幅: 厳密円ではh+1ビット。保持座標M=10^9,尺度<=2なら192ビットで十分、BigInt等。
template <class T> auto center_exact(const Circle<T> &circle) {
    return circle.center_exact();
}
// 半径二乗を厳密に返す。O(1)算術操作。
// 整数幅: 厳密円ではh+1ビット。保持座標M=10^9,尺度<=2なら192ビットで十分、BigInt等。
template <class T> auto radius_squared_exact(const Circle<T> &circle) {
    return circle.radius_squared_exact();
}
// 保持する3点を返す。O(1)回の値コピー。
// 整数幅: 保持データのコピーのみ。追加幅不要。
template <class T> auto points(const Circle<T> &circle) {
    return circle.points();
}
// 保持3点の共通尺度を返す。O(1)回の値コピー。
// 整数幅: 保持データのコピーのみ。追加幅不要。
template <class T> auto point_scale(const Circle<T> &circle) {
    return circle.point_scale();
}
// 尺度で割った3点を厳密に返す。O(1)算術操作。
// 整数幅: 整数の保持座標上限M、尺度絶対値S>=1に対しB=max(M,S)。M=10^9,S<=2ならint可。
template <class T> auto points_exact(const Circle<T> &circle) {
    return circle.points_exact();
}
// 中心をfloat64へ近似する。O(1)算術操作と変換費用。
// 整数幅: 厳密円ではh+1ビット。保持座標M=10^9,尺度<=2なら192ビットで十分、BigInt等。
template <class T> auto center_approx(const Circle<T> &circle) {
    return circle.center_approx();
}
// 半径をfloat64へ近似する。O(1)算術操作と変換費用。
// 整数幅: 厳密円ではh+1ビット。保持座標M=10^9,尺度<=2なら192ビットで十分、BigInt等。
template <class T> double radius_approx(const Circle<T> &circle) {
    return circle.radius_approx();
}

// 同じ座標型の2円盤の共通部分面積を非負の近似値で返す。
// 離れている場合・外接・点円は0、包含・内接・同一円は小さい円盤の面積。
// 戻り値はCircleResultScalar<T>、内部計算はlong double。整数・Fractionも近似する。
// 距離と結果は使用する浮動小数点型で有限に表現できること。
// O(1)算術操作と近似変換費用、追加領域 O(1)。厳密円の中間値は指定基底型に収まること。
template <class T>
// 整数幅: 厳密な中心差を取ってから近似するので2h+2ビットで十分。M=10^9ではBigInt等。
CircleResultScalar<T> intersection_area(const Circle<T> &a, const Circle<T> &b) {
    long double distance, small, large;
    if constexpr (cplib::isFloatingPoint<T>) {
        distance = detail::circleLength(detail::circleFloatPointAs<long double>(b.center()) -
                                        detail::circleFloatPointAs<long double>(a.center()));
        small = a.radius();
        large = b.radius();
    } else {
        auto x = detail::circleExactData(a), y = detail::circleExactData(b);
        auto difference = y.first - x.first;
        // 大きな座標の円でも、中心差を厳密に取ってから近似する。
        distance = detail::circleLength(Point<long double>{detail::circleApprox<long double>(difference.x),
                                                           detail::circleApprox<long double>(difference.y)});
        small = detail::circleApprox<long double>(x.second, true);
        large = detail::circleApprox<long double>(y.second, true);
    }
    if (small > large)
        std::swap(small, large);
    if (small == 0 || distance - large >= small)
        return CircleResultScalar<T>(0);

    long double result;
    if (large - distance >= small) {
        result = (std::acos(-1.0L) * small) * small;
    } else {
        long double scale = std::max(distance, large);
        long double d = distance / scale, r = small / scale, R = large / scale;
        // Heronの公式による共通弦の半分の長さ。接する付近の差を先に取る。
        long double gap = large - small;
        long double inner = small < large / 2
            ? ((distance - large) + small) / distance
            : (distance - gap) / distance;
        long double height = sqrt((small - (distance - large)) / scale) *
                             sqrt(d + R + r) *
                             sqrt(inner) *
                             sqrt(1 + gap / distance) / 2;
        // 小さい円の中心から共通弦までの符号付き距離（scaleで正規化）。
        long double chord = small < large / 2
            ? (((distance - large) / distance) * (d + R) + (r / d) * r) / 2
            : (d - (gap / distance) * (r + R)) / 2;
        auto segment = [](long double angle, long double radius) {
            if (angle < 0.01L) {
                // angle - sin(angle)*cos(angle)の桁落ちを避ける。
                long double squared = angle * angle;
                long double coefficient = 2.0L / 3 + squared * (-2.0L / 15 +
                    squared * (4.0L / 315 - squared * (2.0L / 2835)));
                long double length = angle * radius;
                return ((coefficient * length) * angle) * length;
            }
            return ((angle - std::sin(angle) * std::cos(angle)) * radius) * radius;
        };
        result = segment(atan2(height, chord), small) +
                 segment(atan2(height, d - chord), large);
    }
    auto out = static_cast<CircleResultScalar<T>>(result);
    if (!isfinite(out))
        throw std::invalid_argument("2円の共通部分面積は有限で表現できる必要があります");
    return out;
}

// 点の内外・境界を分類する。floatは相対許容誤差、非floatは厳密。O(1)算術操作。
template <class T, class S>
// 整数幅: 整数3点+整数問い合わせではH=max(1,|保持座標|,|尺度*問い合わせ座標|)、
// S=|尺度|としてB=max(S,96H^4)。
// 3点をそのまま保持する尺度1の円・問い合わせが全て10^9以下なら128ビット、T=__int128で十分。
// 生成された保持座標や尺度が大きい場合はHで再計算する。戻り値は小さくても内部は四乗。
// Fractionを含む経路は中心・二乗半径を用い、冒頭のhに対して9h+7ビットで十分。
CirclePointLocation classify(const Circle<T> &circle, const Point<S> &point, long double tolerance = 1e-10L) {
    static_assert(cplib::isFloatingPoint<T> || !cplib::isFloatingPoint<S>,
                  "厳密円へのfloat入力にはtoFloatCircleで明示的に変換してください");
    if constexpr (cplib::isFloatingPoint<T>) {
        using F = CircleResultScalar<T, S>;
        auto center = detail::circleFloatPointAs<F>(circle.center());
        F radius = circle.radius();
        detail::circleCheckTolerance(tolerance);
        F distance = detail::circleLength(detail::circleFloatPointAs<F>(point) - center);
        if (abs(distance - radius) <= tolerance * std::max(distance, radius))
            return circleBoundary;
        return distance < radius ? circleInside : circleOutside;
    } else if constexpr (!detail::isCircleFraction<T> && !detail::isCircleFraction<S>) {
        using B = detail::CircleBase<T>;
        auto a = circle.points();
        B x = B(point.x) * circle.point_scale(), y = B(point.y) * circle.point_scale();
        if (exact_equal(a[0], a[1]) && exact_equal(a[0], a[2]))
            return x == a[0].x && y == a[0].y ? circleBoundary : circleOutside;
        Int orientation = detail::circleSign(
            (a[1].x - a[0].x) * (a[2].y - a[0].y) - (a[1].y - a[0].y) * (a[2].x - a[0].x));
        B ax = a[0].x - x, ay = a[0].y - y;
        B bx = a[1].x - x, by = a[1].y - y;
        B cx = a[2].x - x, cy = a[2].y - y;
        B determinant = (ax * ax + ay * ay) * (bx * cy - by * cx)
                      - (bx * bx + by * by) * (ax * cy - ay * cx)
                      + (cx * cx + cy * cy) * (ax * by - ay * bx);
        Int sign = detail::circleSign(determinant) * orientation;
        return sign > 0 ? circleInside : sign < 0 ? circleOutside : circleBoundary;
    } else {
        using B = detail::CircleBase<T>;
        auto [center, radiusSquared] = detail::circleExactData(circle);
        auto difference = cplib::norm(detail::circlePointLike<B>(point) - center) - radiusSquared;
        Int sign = detail::circleSign(difference.num);
        return sign > 0 ? circleOutside : sign < 0 ? circleInside : circleBoundary;
    }
}
// 境界を含む円盤の点包含を判定する。O(1)算術操作。
template <class T, class S>
// 整数幅: classifyと同じB=max(S,96H^4)。尺度1の保持3点・問い合わせが10^9以下なら__int128。
// Fractionを含む場合は9h+7ビットで十分。
bool contains(const Circle<T> &circle, const Point<S> &point, long double tolerance = 1e-10L) {
    return classify(circle, point, tolerance) != circleOutside;
}
// 円周上かを判定する。O(1)算術操作。
template <class T, class S>
// 整数幅: classifyと同じB=max(S,96H^4)。尺度1の保持3点・問い合わせが10^9以下なら__int128。
// Fractionを含む場合は9h+7ビットで十分。
bool on_circle(const Circle<T> &circle, const Point<S> &point, long double tolerance = 1e-10L) {
    return classify(circle, point, tolerance) == circleBoundary;
}
namespace detail {
template <class T, class S> auto circleLineData(const Circle<T> &circle, const Line<S> &line) {
    using B = CircleBase<T>;
    auto [center, radiusSquared] = circleExactData(circle);
    auto start = circlePointLike<B>(line.s), direction = circlePointLike<B>(line.t) - start;
    auto a = cplib::norm(direction);
    if (a.num == B(0))
        throw std::invalid_argument("直線は非退化である必要があります");
    auto difference = start - center;
    auto b = cplib::dot(difference, direction) * 2;
    auto discriminant = b * b - a * (cplib::norm(difference) - radiusSquared) * 4;
    using Vector = decltype(start);
    using Scalar = decltype(a);
    struct Data {
        Vector s, v;
        Scalar a, b, discriminant;
    };
    return Data{start, direction, a, b, discriminant};
}
template <class T> auto circlePairData(const Circle<T> &a, const Circle<T> &b) {
    auto x = circleExactData(a), y = circleExactData(b);
    auto direction = y.first - x.first;
    auto distance = cplib::norm(direction);
    using R = decltype(distance);
    R factor(0), height(0);
    if (distance.num != CircleBase<T>(0)) {
        factor = (distance + x.second - y.second) / (distance * 2);
        height = x.second - distance * factor * factor;
    }
    struct Data {
        decltype(direction) center, v;
        R distanceSquared, factor, heightSquared;
    };
    return Data{x.first, direction, distance, factor, height};
}
}
// 円周と直線の交点を返す。非floatは本数を先に厳密判定する。O(1)算術操作と近似変換費用。
template <class T, class S>
// 整数幅: 冒頭のhに対し53h+41ビットで十分（判別式に加え、交点近似用の方向成分の二乗と分数積を含む）。
// 約分を無視した安全側の十分条件。M=10^9の一般の整数3点円ではBigInt等を使う。
auto cross_points(const Circle<T> &circle, const Line<S> &line, long double tolerance = 1e-10L) {
    static_assert(cplib::isFloatingPoint<T> || !cplib::isFloatingPoint<S>,
                  "厳密円へのfloat入力にはtoFloatCircleで明示的に変換してください");
    CircleIntersectionsT<CircleResultScalar<T, S>> result;
    if constexpr (cplib::isFloatingPoint<T>) {
        using F = CircleResultScalar<T, S>;
        auto center = detail::circleFloatPointAs<F>(circle.center());
        F radius = circle.radius();
        detail::circleCheckTolerance(tolerance);
        auto start = detail::circleFloatPointAs<F>(line.s), finish = detail::circleFloatPointAs<F>(line.t);
        detail::circleCheckPoint(start);
        detail::circleCheckPoint(finish);
        auto v = finish - start;
        F length = detail::circleLength(v);
        if (length == 0)
            throw std::invalid_argument("直線は非退化である必要があります");
        auto unit = v / length, w = center - start;
        detail::circleCheckPoint(w);
        F signedDistance = cplib::cross(unit, w), distance = abs(signedDistance);
        if (!isfinite(distance))
            throw std::invalid_argument("直線距離は有限である必要があります");
        F eps = tolerance * std::max(radius, distance);
        if (distance > radius && distance - radius > eps)
            return result;
        auto foot = center - detail::circlePerpendicular(unit) * signedDistance;
        if (abs(distance - radius) <= eps) {
            result.points = {foot};
        } else {
            F ratio = distance / radius;
            F height = radius * sqrt(std::max(F(0), (1 - ratio) * (1 + ratio)));
            result.points = {foot - unit * height, foot + unit * height};
        }
    } else {
        auto data = detail::circleLineData(circle, line);
        Int sign = detail::circleSign(data.discriminant.num);
        if (sign < 0)
            return result;
        auto foot = detail::circlePointApprox(data.s - data.v * (data.b / (data.a * 2)));
        if (sign == 0) {
            result.points = {foot};
        } else {
            auto offset = detail::circleOffsetApprox(data.v, data.discriminant / (data.a * data.a * 4));
            result.points = {foot - offset, foot + offset};
        }
    }
    for (auto point : result.points)
        detail::circleCheckPoint(point);
    return result;
}
// 円周と閉線分の交点を返す。非floatは端点の採否も厳密に判定する。O(1)算術操作と近似変換費用。
template <class T, class S>
// 整数幅: 冒頭のhに対し65h+54ビットで十分（直線との交点計算と、端点との大小判定で平方根を消す二乗を含む）。
// 約分を無視した安全側の十分条件。M=10^9の一般の整数3点円ではBigInt等を使う。
auto cross_points(const Circle<T> &circle, const Segment<S> &segment, long double tolerance = 1e-10L) {
    static_assert(cplib::isFloatingPoint<T> || !cplib::isFloatingPoint<S>,
                  "厳密円へのfloat入力にはtoFloatCircleで明示的に変換してください");
    CircleIntersectionsT<CircleResultScalar<T, S>> result;
    if constexpr (cplib::isFloatingPoint<T>) {
        using F = CircleResultScalar<T, S>;
        detail::circleCheckTolerance(tolerance);
        auto start = detail::circleFloatPointAs<F>(segment.s), finish = detail::circleFloatPointAs<F>(segment.t);
        auto direction = finish - start;
        F length = detail::circleLength(direction);
        if (length == 0) {
            if (on_circle(circle, segment.s, tolerance))
                result.points.push_back(start);
            return result;
        }
        auto unit = direction / length;
        auto candidates = cross_points(circle, Line<F>{start, finish}, tolerance);
        for (auto p : candidates.points) {
            F parameter = cplib::dot(p - start, unit);
            if (parameter >= -tolerance * length && parameter - length <= tolerance * length)
                result.points.push_back(p);
        }
    } else {
        if (exact_equal(segment.s, segment.t)) {
            if (on_circle(circle, segment.s))
                result.points.push_back(detail::circlePointApprox(detail::circleExactPoint(segment.s)));
            return result;
        }
        auto data = detail::circleLineData(circle, Line<S>{segment.s, segment.t});
        Int sign = detail::circleSign(data.discriminant.num);
        if (sign < 0)
            return result;
        std::vector<Int> accepted;
        Int index = 0;
        for (Int direction : {-1, 1}) {
            if (direction == 1 && sign == 0)
                continue;
            if (detail::circleRadicalSign(-data.b, direction, data.discriminant) >= 0 &&
                detail::circleRadicalSign(-data.b - data.a * 2, direction, data.discriminant) <= 0)
                accepted.push_back(index);
            ++index;
        }
        if (accepted.empty())
            return result;
        auto candidates = cross_points(circle, Line<S>{segment.s, segment.t});
        for (Int i : accepted)
            result.points.push_back(candidates.points[i]);
    }
    return result;
}
// 2円周の交点を返す。同一の正半径円は無限個。O(1)算術操作と近似変換費用。
template <class T>
// 整数幅: 冒頭のhに対し57h+40ビットで十分（中心差・交点の足・高さ二乗・近似用の分数積を含む）。
// 約分を無視した安全側の十分条件。M=10^9の一般の整数3点円ではBigInt等を使う。
auto cross_points(const Circle<T> &a, const Circle<T> &b, long double tolerance = 1e-10L) {
    CircleIntersectionsT<CircleResultScalar<T>> result;
    if constexpr (cplib::isFloatingPoint<T>) {
        using F = T;
        detail::circleCheckTolerance(tolerance);
        auto v = b.center() - a.center();
        F distance = detail::circleLength(v);
        if (distance == 0) {
            if (a.radius() == b.radius()) {
                if (a.radius() == 0)
                    result.points = {a.center()};
                else
                    result.kind = circleInfinite;
            }
            return result;
        }
        F scale = std::max({distance, a.radius(), b.radius()});
        F dn = distance / scale, ra = a.radius() / scale, rb = b.radius() / scale;
        if (dn == 0)
            throw std::invalid_argument("円の尺度比は使用する浮動小数点型で表現可能である必要があります");
        F sum = ra + rb, difference = abs(ra - rb);
        if (dn > sum + tolerance * std::max(dn, sum) ||
            dn < difference - tolerance * std::max(dn, difference))
            return result;
        auto unit = v / distance;
        F x = (dn + (ra - rb) * (ra + rb) / dn) * 0.5 * scale;
        auto foot = a.center() + unit * x;
        if (abs(dn - sum) <= tolerance * std::max(dn, sum) ||
            abs(dn - difference) <= tolerance * std::max(dn, difference)) {
            result.points = {foot};
        } else {
            F height = scale * F(0.5) * sqrt(std::max(F(0), (sum + dn) * (sum - dn)))
                * sqrt(std::max(F(0), ((dn + difference) / dn) * ((dn - difference) / dn)));
            auto offset = detail::circlePerpendicular(unit) * height;
            result.points = {foot - offset, foot + offset};
        }
    } else {
        auto data = detail::circlePairData(a, b);
        if (data.distanceSquared.num == detail::CircleBase<T>(0)) {
            if (a == b) {
                if (a.radius_squared_exact().num == detail::CircleBase<T>(0))
                    result.points = {a.center_approx()};
                else
                    result.kind = circleInfinite;
            }
            return result;
        }
        Int sign = detail::circleSign(data.heightSquared.num);
        if (sign < 0)
            return result;
        auto foot = detail::circlePointApprox(data.center + data.v * data.factor);
        if (sign == 0) {
            result.points = {foot};
        } else {
            auto offset = detail::circleOffsetApprox(Point{-data.v.y, data.v.x},
                data.heightSquared / data.distanceSquared);
            result.points = {foot - offset, foot + offset};
        }
    }
    for (auto p : result.points)
        detail::circleCheckPoint(p);
    return result;
}
// 直線との交点数0〜2を返す。非floatは厳密。O(1)算術操作。
template <class T, class S>
// 整数幅: 冒頭のhに対し33h+27ビットで十分（分数の判別式b*b-4*a*cを含む）。
// 約分を無視した安全側の十分条件。M=10^9の一般の整数3点円ではBigInt等を使う。
Int intersection_count(const Circle<T> &circle, const Line<S> &line, long double tolerance = 1e-10L) {
    static_assert(cplib::isFloatingPoint<T> || !cplib::isFloatingPoint<S>,
                  "厳密円へのfloat入力にはtoFloatCircleで明示的に変換してください");
    if constexpr (cplib::isFloatingPoint<T>)
        return cross_points(circle, line, tolerance).points.size();
    else {
        Int sign = detail::circleSign(detail::circleLineData(circle, line).discriminant.num);
        return sign < 0 ? 0 : sign == 0 ? 1 : 2;
    }
}
// 閉線分との交点数0〜2を返す。端点も非floatでは厳密に扱う。O(1)算術操作。
template <class T, class S>
// 整数幅: 冒頭のhに対し65h+54ビットで十分（判別式と、端点との大小判定で平方根を消す二乗を含む）。
// 約分を無視した安全側の十分条件。M=10^9の一般の整数3点円ではBigInt等を使う。
Int intersection_count(const Circle<T> &circle, const Segment<S> &segment, long double tolerance = 1e-10L) {
    static_assert(cplib::isFloatingPoint<T> || !cplib::isFloatingPoint<S>,
                  "厳密円へのfloat入力にはtoFloatCircleで明示的に変換してください");
    if constexpr (cplib::isFloatingPoint<T>) {
        return cross_points(circle, segment, tolerance).points.size();
    } else {
        if (exact_equal(segment.s, segment.t))
            return on_circle(circle, segment.s);
        auto data = detail::circleLineData(circle, Line<S>{segment.s, segment.t});
        Int sign = detail::circleSign(data.discriminant.num), count = 0;
        if (sign < 0)
            return 0;
        for (Int direction : {-1, 1}) {
            if (direction == 1 && sign == 0)
                continue;
            if (detail::circleRadicalSign(-data.b, direction, data.discriminant) >= 0 &&
                detail::circleRadicalSign(-data.b - data.a * 2, direction, data.discriminant) <= 0)
                ++count;
        }
        return count;
    }
}
// 2円周の交点数を返す。同一正半径円は-1、同一点円は1。O(1)算術操作。
template <class T>
// 整数幅: 冒頭のhに対し45h+33ビットで十分（中心差の二乗・分数の比率・高さ二乗を含む）。
// 約分を無視した安全側の十分条件。M=10^9の一般の整数3点円ではBigInt等を使う。
Int intersection_count(const Circle<T> &a, const Circle<T> &b, long double tolerance = 1e-10L) {
    if constexpr (cplib::isFloatingPoint<T>) {
        auto result = cross_points(a, b, tolerance);
        return result.kind == circleInfinite ? -1 : Int(result.points.size());
    } else {
        auto data = detail::circlePairData(a, b);
        if (data.distanceSquared.num == detail::CircleBase<T>(0)) {
            if (a != b)
                return 0;
            return a.radius_squared_exact().num == detail::CircleBase<T>(0) ? 1 : -1;
        }
        Int sign = detail::circleSign(data.heightSquared.num);
        return sign < 0 ? 0 : sign == 0 ? 1 : 2;
    }
}
// 点からの接線数を返す。点円自身は-1、それ以外の点から点円へは1。O(1)算術操作。
template <class T, class S>
// 整数幅: max(classifyの必要幅,h+1)。半径0の確認で中心・半径二乗も計算する。
// 整数classifyだけで使える__int128では保証できない。M=10^9の一般の3点円ではBigInt等。
Int tangent_count(const Circle<T> &circle, const Point<S> &p, long double tolerance = 1e-10L) {
    auto location = classify(circle, p, tolerance);
    bool zero;
    if constexpr (cplib::isFloatingPoint<T>)
        zero = circle.radius() == 0;
    else
        zero = circle.radius_squared_exact().num == detail::CircleBase<T>(0);
    if (zero)
        return location == circleBoundary ? -1 : 1;
    return location == circleOutside ? 2 : location == circleBoundary ? 1 : 0;
}
// 共通接線数0〜4を返す。同一円は-1。点円も扱う。非floatは厳密。O(1)算術操作。
template <class T>
// 整数幅: 冒頭のhに対し22h+18ビットで十分（(中心間距離二乗-r1二乗-r2二乗)の二乗と分数の差を含む）。
// 約分を無視した安全側の十分条件。M=10^9の一般の整数3点円ではBigInt等を使う。
Int common_tangent_count(const Circle<T> &a, const Circle<T> &b, long double tolerance = 1e-10L) {
    if constexpr (cplib::isFloatingPoint<T>) {
        using F = T;
        detail::circleCheckTolerance(tolerance);
        F d = detail::circleLength(b.center() - a.center());
        if (d == 0)
            return a.radius() == b.radius() ? -1 : 0;
        if (a.radius() == 0 && b.radius() == 0)
            return 1;
        F scale = std::max({d, a.radius(), b.radius()});
        F dn = d / scale, ra = a.radius() / scale, rb = b.radius() / scale;
        Int count = 0;
        for (F side : {F(1), F(-1)}) {
            if (side < 0 && (a.radius() == 0 || b.radius() == 0))
                continue;
            F delta = abs(ra - side * rb), eps = tolerance * std::max(dn, delta);
            if (delta > dn + eps)
                continue;
            count += abs(delta - dn) <= eps ? 1 : 2;
        }
        return count;
    } else {
        auto x = detail::circleExactData(a), y = detail::circleExactData(b);
        auto distance = cplib::norm(y.first - x.first);
        using B = detail::CircleBase<T>;
        if (distance.num == B(0))
            return a == b ? -1 : 0;
        bool firstZero = x.second.num == B(0), secondZero = y.second.num == B(0);
        if (firstZero && secondZero)
            return 1;
        auto difference = distance - x.second - y.second;
        if (firstZero || secondZero) {
            Int sign = detail::circleSign(difference.num);
            return sign < 0 ? 0 : sign == 0 ? 1 : 2;
        }
        Int comparison = detail::circleSign((difference * difference - x.second * y.second * 4).num);
        if (comparison < 0)
            return 2;
        return comparison == 0 ? (difference.num > B(0) ? 3 : 1) : (difference.num > B(0) ? 4 : 0);
    }
}
// 共通接線を接点と単位方向で返す。存在と重複は非floatでは厳密判定する。
// O(1)算術操作と近似変換費用。同一円ではkind=circleInfinite。
template <class T>
// 整数幅: 冒頭のhに対し22h+18ビットで十分（接線の存在判定の二乗と、近似用の分数積・商を含む）。
// 約分を無視した安全側の十分条件。M=10^9の一般の整数3点円ではBigInt等を使う。
auto common_tangents(const Circle<T> &a, const Circle<T> &b, long double tolerance = 1e-10L) {
    CircleTangentsT<CircleResultScalar<T>> result;
    if constexpr (cplib::isFloatingPoint<T>) {
        using F = T;
        detail::circleCheckTolerance(tolerance);
        auto v = b.center() - a.center();
        F distance = detail::circleLength(v);
        if (distance == 0) {
            if (a.radius() == b.radius())
                result.kind = circleInfinite;
            return result;
        }
        auto unit = v / distance;
        F scale = std::max({distance, a.radius(), b.radius()});
        F dn = distance / scale;
        if (dn == 0)
            throw std::invalid_argument("円の尺度比は使用する浮動小数点型で表現可能である必要があります");
        for (F side : {F(1), F(-1)}) {
            if (side < 0 && (a.radius() == 0 || b.radius() == 0))
                continue;
            F delta = a.radius() / scale - side * (b.radius() / scale);
            if (abs(delta) > dn + tolerance * std::max(dn, abs(delta)))
                continue;
            F ratio = std::clamp(delta / dn, F(-1), F(1));
            F height = sqrt(std::max(F(0), (1 - ratio) * (1 + ratio)));
            bool tangent = abs(abs(delta) - dn) <= tolerance * std::max(dn, abs(delta));
            for (F sign : {F(-1), F(1)}) {
                if (sign > 0 && (tangent || (a.radius() == 0 && b.radius() == 0)))
                    continue;
                F normalRatio = tangent ? (delta < 0 ? -1.0 : 1.0) : ratio;
                auto normal = unit * normalRatio + detail::circlePerpendicular(unit) * (tangent ? 0 : height * sign);
                auto first = a.center() + normal * a.radius();
                auto second = b.center() + normal * (side * b.radius());
                result.tangents.push_back({first, second, detail::circlePerpendicular(normal)});
            }
        }
    } else {
        using B = detail::CircleBase<T>;
        auto x = detail::circleExactData(a), y = detail::circleExactData(b);
        auto v = y.first - x.first;
        auto distance = cplib::norm(v);
        using R = decltype(distance);
        if (distance.num == B(0)) {
            if (a == b)
                result.kind = circleInfinite;
            return result;
        }
        bool firstZero = x.second.num == B(0), secondZero = y.second.num == B(0);
        struct Family { Int side, existence, contactSign; };
        std::vector<Family> families;
        if (!(firstZero && secondZero)) {
            auto difference = distance - x.second - y.second;
            auto radicand = x.second * y.second * 4;
            for (Int side : {1, -1}) {
                if (side < 0 && (firstZero || secondZero))
                    continue;
                Int existence = detail::circleRadicalSign(difference, side, radicand);
                if (existence < 0)
                    continue;
                Int contactSign = side > 0 && (x.second - y.second).num < B(0) ? -1 : 1;
                families.push_back({side, existence, contactSign});
            }
            if (families.empty())
                return result;
        }
        auto unit = detail::circleOffsetApprox(v, R(1) / distance);
        auto firstCenter = detail::circlePointApprox(x.first), secondCenter = detail::circlePointApprox(y.first);
        if (firstZero && secondZero) {
            result.tangents.push_back({firstCenter, secondCenter, unit});
            return result;
        }
        double firstRadius = detail::circleApprox(x.second, true);
        double secondRadius = detail::circleApprox(y.second, true);
        for (auto family : families) {
            double ratio = detail::circleApprox(x.second / distance, true)
                         - family.side * detail::circleApprox(y.second / distance, true);
            double cosine = family.existence == 0 ? double(family.contactSign) : std::clamp(ratio, -1.0, 1.0);
            double sine = family.existence == 0 ? 0 : sqrt(std::max(0.0, (1 - cosine) * (1 + cosine)));
            for (double sign : {-1.0, 1.0}) {
                if (sign > 0 && family.existence == 0)
                    continue;
                auto normal = unit * cosine + detail::circlePerpendicular(unit) * (sine * sign);
                auto first = firstCenter + normal * firstRadius;
                auto second = secondCenter + normal * (family.side * secondRadius);
                result.tangents.push_back({first, second, detail::circlePerpendicular(normal)});
            }
        }
    }
    for (auto tangent : result.tangents) {
        detail::circleCheckPoint(tangent.first);
        detail::circleCheckPoint(tangent.second);
        detail::circleCheckPoint(tangent.direction);
    }
    return result;
}
// 点からの接線を返す。点円自身では無限個。O(1)算術操作と近似変換費用。
template <class T, class S>
// 整数幅: 冒頭のhに対し30h+19ビットで十分（接線判定と、接点・単位方向を近似する前の分数積・商を含む）。
// 約分を無視した安全側の十分条件。M=10^9の一般の整数3点円ではBigInt等を使う。
auto tangent_lines(const Circle<T> &circle, const Point<S> &point, long double tolerance = 1e-10L) {
    if constexpr (cplib::isFloatingPoint<T>) {
        using F = CircleResultScalar<T, S>;
        Circle<F> promoted(detail::circleFloatPointAs<F>(circle.center()), F(circle.radius()));
        return common_tangents(promoted, Circle<F>(detail::circleFloatPointAs<F>(point), 0), tolerance);
    } else {
        using B = detail::CircleBase<T>;
        Int count = tangent_count(circle, point);
        CircleTangentsT<CircleResultScalar<T>> result;
        if (count == -1) {
            result.kind = circleInfinite;
            return result;
        }
        if (count == 0)
            return result;
        auto [center, radiusSquared] = detail::circleExactData(circle);
        auto q = detail::circlePointLike<B>(point), v = q - center;
        auto distanceSquared = cplib::norm(v);
        using R = decltype(distanceSquared);
        if (radiusSquared.num == B(0)) {
            result.tangents.push_back({detail::circlePointApprox(center), detail::circlePointApprox(q),
                detail::circleOffsetApprox(v, R(1) / distanceSquared)});
            return result;
        }
        auto foot = center + v * (radiusSquared / distanceSquared);
        auto squaredFactor = radiusSquared * (distanceSquared - radiusSquared) / (distanceSquared * distanceSquared);
        auto perpendicular = Point{-v.y, v.x};
        auto offset = detail::circleOffsetApprox(perpendicular, squaredFactor);
        auto normalBase = detail::circleOffsetApprox(v, radiusSquared / (distanceSquared * distanceSquared));
        auto normalOffset = detail::circleOffsetApprox(perpendicular,
            (distanceSquared - radiusSquared) / (distanceSquared * distanceSquared));
        for (double sign : {-1.0, 1.0}) {
            if (count == 1 && sign > 0)
                continue;
            auto normal = normalBase + normalOffset * sign;
            auto contact = detail::circlePointApprox(foot) + offset * sign;
            detail::circleCheckPoint(contact);
            result.tangents.push_back({contact, detail::circlePointApprox(q), detail::circlePerpendicular(normal)});
        }
        return result;
    }
}
// 円・直線・閉線分との交点をfloat64へ近似して返す。O(1)算術操作と変換費用。
template <class T, class Object>
// 整数幅: cross_pointsと同じ。直線53h+41、線分65h+54、円57h+40ビットで十分。M=10^9ではBigInt等。
CircleIntersections cross_points_approx(const Circle<T> &circle, const Object &object) {
    auto source = cross_points(circle, object);
    CircleIntersections result;
    result.kind = source.kind;
    for (auto p : source.points)
        result.points.push_back({double(p.x), double(p.y)});
    return result;
}
// 整数幅: tangent_linesと同じ30h+19ビットで十分。M=10^9の一般の整数3点円ではBigInt等。
template <class T, class S> CircleTangents tangent_lines_approx(const Circle<T> &circle, const Point<S> &p) {
    auto source = tangent_lines(circle, p);
    CircleTangents result;
    result.kind = source.kind;
    for (auto t : source.tangents)
        result.tangents.push_back({{double(t.first.x), double(t.first.y)},
            {double(t.second.x), double(t.second.y)}, {double(t.direction.x), double(t.direction.y)}});
    return result;
}
}

namespace cplib::detail {
template <class Derived>
template <class S>
CirclePointLocation CircleOperations<Derived>::classify(Point<S> p, long double tolerance) const {
    return cplib::classify(static_cast<const Derived &>(*this), p, tolerance);
}
template <class Derived>
template <class S>
bool CircleOperations<Derived>::contains(Point<S> p, long double tolerance) const {
    return cplib::contains(static_cast<const Derived &>(*this), p, tolerance);
}
template <class Derived>
template <class S>
bool CircleOperations<Derived>::on_circle(Point<S> p, long double tolerance) const {
    return cplib::on_circle(static_cast<const Derived &>(*this), p, tolerance);
}
template <class Derived>
template <class Object>
auto CircleOperations<Derived>::cross_points(const Object &object, long double tolerance) const {
    return cplib::cross_points(static_cast<const Derived &>(*this), object, tolerance);
}
template <class Derived>
template <class Object>
Int CircleOperations<Derived>::intersection_count(const Object &object, long double tolerance) const {
    return cplib::intersection_count(static_cast<const Derived &>(*this), object, tolerance);
}
template <class Derived>
template <class S>
auto CircleOperations<Derived>::tangent_lines(Point<S> p, long double tolerance) const {
    return cplib::tangent_lines(static_cast<const Derived &>(*this), p, tolerance);
}
template <class Derived>
template <class S>
Int CircleOperations<Derived>::tangent_count(Point<S> p, long double tolerance) const {
    return cplib::tangent_count(static_cast<const Derived &>(*this), p, tolerance);
}
template <class Derived>
auto CircleOperations<Derived>::common_tangents(const Derived &other, long double tolerance) const {
    return cplib::common_tangents(static_cast<const Derived &>(*this), other, tolerance);
}
template <class Derived>
Int CircleOperations<Derived>::common_tangent_count(const Derived &other, long double tolerance) const {
    return cplib::common_tangent_count(static_cast<const Derived &>(*this), other, tolerance);
}
}

namespace cplib::detail {
template <class Derived>
template <class Object>
CircleIntersections CircleOperations<Derived>::cross_points_approx(const Object &object) const {
    return cplib::cross_points_approx(static_cast<const Derived &>(*this), object);
}
template <class Derived>
template <class S>
CircleTangents CircleOperations<Derived>::tangent_lines_approx(Point<S> p) const {
    return cplib::tangent_lines_approx(static_cast<const Derived &>(*this), p);
}
}
