#pragma once
#include <cplib/common.hpp>
#include <cplib/math/floating_point.hpp>
#include <cmath>
#include <functional>
#include <istream>
#include <iomanip>
#include <sstream>
#include <type_traits>

namespace cplib {
template <class T, class Id = void> struct Point;

inline long double GEOMETRY_EPS = 1e-10L;

// 整数版の幅の表記（geometry全体で共通）:
// M=max(1, 全入力座標の絶対値)、W(B)=1+ceil(log2(B+1))（符号込みのビット数）。
// 「整数幅: B=...」は全中間値の絶対値の安全側上限で、W(B)ビットあれば十分。
// 最小の必要幅とは限らない。Mには問い合わせ点・直線や線分の端点も含める。
// M<2^kなら、B=c*M^dに対しd*k+ceil(log2(c))+1ビットが十分条件となる。
// 明記した関数以外はライブラリ側で計算用の型を広げない。int/long long/__int128を32/64/128ビットとすると、
// B<=2^31-1 / 2^63-1 / 2^127-1で選ぶ。Intは64ビット。
// std::is_integral_vで__int128を判定するAPIはGNU拡張が有効な環境（-std=gnu++20等）を前提とする。
// 例: M=10^9ではW(M)=31、W(2M)=32、W(2M^2)=62、W(8M^2)=64。
// 結果をlong longに代入するだけではintの中間計算は広がらない。入力のPoint<T>から広げる。
// 各幅は符号付き整数用。Fraction<T>には適用せず、約分前の分子・分母の演算も別途見積もる。
// 符号反転・abs・除算で最小負整数を扱う場合も、反転後が型に収まること。

// 浮動小数点の場合は許容誤差付き、それ以外は厳密に等値を判定する。計算量 O(1)。
// 絶対誤差と両辺の絶対値に対する相対誤差を使う。
// 整数幅: 値の比較のみ。入力を保持できればよい（M=10^9ならint可）。
template <class T, class S> bool geometry_eq(const T &x, const S &y) {
    if constexpr (cplib::isFloatingPoint<T> || cplib::isFloatingPoint<S>) {
        using Common = std::common_type_t<T, S>;
        Common a = x, b = y;
        if (a == b)
            return true;
        if (!isfinite(a) || !isfinite(b))
            return false;
        auto diff = abs(a - b);
        return diff < GEOMETRY_EPS || diff < GEOMETRY_EPS *
            std::max(abs(a), abs(b));
    } else
        return x == y;
}

// 許容誤差付き等値判定の否定を返す。計算量 O(1)。
// 整数幅: 値の比較のみ。入力を保持できればよい（M=10^9ならint可）。
template <class T, class S> bool geometry_neq(const T &x, const S &y) {
    return !geometry_eq(x, y);
}

// 許容誤差を考慮して x >= y を判定する。計算量 O(1)。
// 整数幅: 値の比較のみ。入力を保持できればよい（M=10^9ならint可）。
template <class T, class S> bool geometry_ge(const T &x, const S &y) {
    return geometry_eq(x, y) || x > y;
}

// 許容誤差を考慮して x <= y を判定する。計算量 O(1)。
// 整数幅: 値の比較のみ。入力を保持できればよい（M=10^9ならint可）。
template <class T, class S> bool geometry_le(const T &x, const S &y) {
    return geometry_eq(x, y) || x < y;
}

// 許容誤差を考慮して x > y を判定する。計算量 O(1)。
// 整数幅: 値の比較のみ。入力を保持できればよい（M=10^9ならint可）。
template <class T, class S> bool geometry_gt(const T &x, const S &y) {
    return !geometry_le(x, y);
}

// 許容誤差を考慮して x < y を判定する。計算量 O(1)。
// 整数幅: 値の比較のみ。入力を保持できればよい（M=10^9ならint可）。
template <class T, class S> bool geometry_lt(const T &x, const S &y) {
    return !geometry_ge(x, y);
}

template <class T> struct PointKey {
    T x, y;
    // ハッシュ用キーを厳密に比較する。計算量 O(1)。
    bool operator==(const PointKey &) const = default;
};

// 座標型の演算を定数時間とすると、座標の加減算、スカラー倍・除算、内積、辞書順比較はいずれも O(1)。
// ストリーム出力は座標の出力に要する時間に比例する。
// 整数幅: 比較・コピー・入出力は追加幅不要。各演算の幅は下記。
template <class T> struct Point<T, void> {
    T x{}, y{};

    // 整数幅: B=M（コピー）。
    Point operator+() const {
        return *this;
    }

    // 整数幅: B=M（符号反転）。最小負整数は不可。
    Point operator-() const {
        return {-x, -y};
    }

    // 整数幅: B=2M。M=10^9ならint可。
    Point &operator+=(Point q) {
        x += q.x;
        y += q.y;
        return *this;
    }

    // 整数幅: B=2M。M=10^9ならint可。
    Point &operator-=(Point q) {
        x -= q.x;
        y -= q.y;
        return *this;
    }

    // 整数幅: C=max(1,|c|)としてB=M*C。M=C=10^9ならlong long。結果はTに格納する。
    template <class S> Point &operator*=(S c) {
        if constexpr (cplib::isFloatingPoint<T> && !std::is_arithmetic_v<S>) {
            x *= static_cast<T>(c);
            y *= static_cast<T>(c);
        } else {
            x *= c;
            y *= c;
        }
        return *this;
    }

    // 整数幅: 整数c!=0ならB=M。M=10^9ならint可。整数除算で0方向へ丸める。
    template <class S> Point &operator/=(S c) {
        if constexpr (cplib::isFloatingPoint<T> && !std::is_arithmetic_v<S>) {
            x /= static_cast<T>(c);
            y /= static_cast<T>(c);
        } else {
            x /= c;
            y /= c;
        }
        return *this;
    }

    // 整数幅: B=2M。M=10^9ならint可（operator+=と同じ）。
    friend Point operator+(Point p, Point q) {
        return p += q;
    }

    // 整数幅: B=2M。M=10^9ならint可（operator-=と同じ）。
    friend Point operator-(Point p, Point q) {
        return p -= q;
    }

    template <class S>
        requires(!std::is_base_of_v<Point, S>)
    // 整数幅: B=M*max(1,|c|)。M=|c|=10^9ならlong long。
    friend Point operator*(Point p, S c) {
        return p *= c;
    }

    template <class S>
        requires(!std::is_base_of_v<Point, S>)
    // 整数幅: B=M*max(1,|c|)。M=|c|=10^9ならlong long。
    friend Point operator*(S c, Point p) {
        return p *= c;
    }

    // 整数幅: 整数c!=0ならB=M。M=10^9ならint可。
    template <class S> friend Point operator/(Point p, S c) {
        return p /= c;
    }

    // 整数幅: 内積なのでB=2M^2。M=10^9なら62ビット、long long。
    friend T operator*(Point p, Point q) {
        return p.x * q.x + p.y * q.y;
    }

    friend bool operator<(Point p, Point q) {
        return p.x == q.x ? p.y < q.y : p.x < q.x;
    }

    friend bool operator>(Point p, Point q) {
        return q < p;
    }

    friend bool operator==(Point p, Point q) {
        return geometry_eq(p.x, q.x) && geometry_eq(p.y, q.y);
    }

    friend bool operator<=(Point p, Point q) {
        return !(p > q);
    }

    friend bool operator>=(Point p, Point q) {
        return !(p < q);
    }

    friend std::istream &operator>>(std::istream &in, Point &p) {
        return in >> p.x >> p.y;
    }

    friend std::ostream &operator<<(std::ostream &out, Point p) {
        return out << p.x << ' ' << p.y;
    }

    // EPSに依存しないハッシュ用キーに変換する。計算量 O(1)。
    // 整数幅: 追加の座標演算なし。B=M、M=10^9ならint可。
    PointKey<T> toPointKey() const {
        return {x, y};
    }

    // 点を文字列に変換する。計算量は座標の出力に要する時間に比例する。
    // 整数幅: 追加の座標演算なし。B=M、M=10^9ならint可。
    std::string to_string() const {
        std::ostringstream out;
        if constexpr (cplib::isFloatingPoint<T>)
            out << std::setprecision(cplib::NumericLimits<T>::max_digits10);
        out << *this;
        return out.str();
    }
};

// IDは入力点の付加情報。座標の比較・入出力には使わず、座標演算で新しく作る点には付かない。
template <class T, class Id> struct Point : Point<T> {
    Id id{};

    Point() = default;
    Point(T x, T y, Id id = Id{}) : Point<T>{std::move(x), std::move(y)}, id(std::move(id)) {}

    Point &operator+=(Point<T> q) {
        Point<T>::operator+=(std::move(q));
        return *this;
    }
    Point &operator-=(Point<T> q) {
        Point<T>::operator-=(std::move(q));
        return *this;
    }
    template <class S> Point &operator*=(S c) {
        Point<T>::operator*=(std::move(c));
        return *this;
    }
    template <class S> Point &operator/=(S c) {
        Point<T>::operator/=(std::move(c));
        return *this;
    }
};

template <class T> Point(T) -> Point<T>;
template <class T> Point(T, T) -> Point<T>;
template <class T, class Id> Point(T, T, Id) -> Point<T, Id>;

// 座標を厳密に比較する。計算量 O(1)。
// 整数幅: 追加の座標演算なし。B=M、M=10^9ならint可。
template <class T> bool exact_equal(Point<T> p, Point<T> q) {
    return p.x == q.x && p.y == q.y;
}
// 座標の近似一致を判定する。計算量 O(1)。
// 整数幅: 追加の座標演算なし。B=M、M=10^9ならint可。
template <class T> bool almost_equal(Point<T> p, Point<T> q) {
    return geometry_eq(p.x, q.x) && geometry_eq(p.y, q.y);
}

template <class T> PointKey<T> toPointKey(Point<T> p) {
    return p.toPointKey();
}

// 点を (x, y) で初期化する。計算量 O(1)。
// 整数幅: 追加の座標演算なし。B=M、M=10^9ならint可。
template <class T> auto initPoint(T x, T y) {
    return Point<T>{x, y};
}

// 座標とIDを指定して点を初期化する。座標のみを使う演算の計算量はPoint<T>と同じ。
template <class T, class Id> auto initPoint(T x, T y, Id id) {
    return Point<T, Id>{std::move(x), std::move(y), std::move(id)};
}

// 点を(p.first, p.second)で初期化する。計算量 O(1)。
// 整数幅: 追加の座標演算なし。B=M、M=10^9ならint可。
template <class T> auto initPoint(std::pair<T, T> p) {
    return Point<T>{p.first, p.second};
}

// 2点p,qの内積 (p.x * q.x + p.y * q.y) を返す。計算量 O(1)。
// 整数幅: B=2M^2。M=10^9なら62ビット、T=long longで十分（int不可）。
template <class T> T dot(Point<T> p, Point<T> q) {
    return p * q;
}

// 2点p,qの外積 (p.x * q.y - p.y * q.x) を返す。計算量 O(1)。
// 整数幅: B=2M^2。M=10^9なら62ビット、T=long longで十分（int不可）。
template <class T> T cross(Point<T> p, Point<T> q) {
    return p.x * q.y - p.y * q.x;
}

// pのノルムの2乗 (p.x * p.x + p.y * p.y) を返す。計算量 O(1)。
// 整数幅: B=2M^2。M=10^9なら62ビット、T=long longで十分（int不可）。
template <class T> T norm(Point<T> p) {
    return dot(p, p);
}

// 点を辞書順に比較し、-1、0、1を返す。計算量 O(1)。
// 整数幅: 追加の座標演算なし。B=M、M=10^9ならint可。
template <class T> Int cmp(Point<T> p, Point<T> q) {
    return p < q ? -1 : q < p ? 1 : 0;
}

template <class T> std::string to_string(Point<T> p) {
    return p.to_string();
}

// 直線を2点で表す。ストリーム出力は座標の出力に要する時間に比例する。
template <class T> struct Line {
    Point<T> s, t;

    friend std::ostream &operator<<(std::ostream &out, const Line &l) {
        return out << '(' << l.s << ", " << l.t << ')';
    }

    // 直線の接ベクトル（直線の方向のベクトル）を返す。計算量 O(1)。
    // 整数幅: 座標差を取るのでB=2M。M=10^9なら32ビット、int可。
    Point<T> toVector() const {
        return t - s;
    }

    // ax+by+c=0 の係数を返す。s->t の左閉半平面は ax+by+c<=0。O(1)演算。
    // 正規化・型の昇格は行わない。差・積・和が T に収まること。
    // 整数幅: B=4M^2（差は2M、差と座標の積の和）。M=10^9なら63ビット、long long。
    std::array<T, 3> toCoefficients() const {
        T a = t.y - s.y;
        T b = s.x - t.x;
        return {a, b, -a * s.x - b * s.y};
    }

    // ax+by+c=0 を表す整数係数 (a,b,c) で、gcd(a,b,c)=1 を満たす辞書順最小のものを返す。
    // 計算量 O(log(M+1))。
    // 端点は異なること。整数幅は toCoefficients と同じ B=4M^2。M=10^9 なら long long。
    std::array<T, 3> toNormalizedCoefficients() const
        requires(std::is_integral_v<T> && std::is_signed_v<T>) {
        auto coefficients = toCoefficients();
        auto &[a, b, c] = coefficients;
        assert(a != 0 || b != 0);
        // c は a, b の整数線形結合なので gcd(a,b,c)=gcd(a,b)。
        T g = std::gcd(a, b);
        if (a > 0 || (a == 0 && b > 0))
            g = -g;
        a /= g;
        b /= g;
        c /= g;
        return coefficients;
    }

    // 直線を文字列に変換する。計算量は座標の出力に要する時間に比例する。
    // 整数幅: 追加の座標演算なし。B=M、M=10^9ならint可。
    std::string to_string() const {
        std::ostringstream out;
        if constexpr (cplib::isFloatingPoint<T>)
            out << std::setprecision(cplib::NumericLimits<T>::max_digits10);
        out << *this;
        return out.str();
    }
};

// 異なる2点 (s, t) を通る直線を初期化する。計算量 O(1)。
// 整数幅: 追加の座標演算なし。B=M、M=10^9ならint可。
template <class T> auto initLine(Point<T> s, Point<T> t) {
    assert(s != t);
    return Line<T>{s, t};
}

// 直線ax + by + c = 0を初期化する。整数型には使用できない。計算量 O(1)。
template <class T>
    requires(!std::is_integral_v<T>)
Line<T> initLine(T a, T b, T c) {
    assert(geometry_neq(a, T(0)) || geometry_neq(b, T(0)));
    if (geometry_eq(b, T(0)))
        return {{-c / a, T(0)}, {-c / a, T(1)}};
    return {{T(0), -c / b}, {T(1), (-a - c) / b}};
}

template <class T> Point<T> toVector(const Line<T> &l) {
    return l.toVector();
}

template <class T> std::array<T, 3> toCoefficients(const Line<T> &l) {
    return l.toCoefficients();
}

template <class T>
    requires(std::is_integral_v<T> && std::is_signed_v<T>)
std::array<T, 3> toNormalizedCoefficients(const Line<T> &l) {
    return l.toNormalizedCoefficients();
}

template <class T> std::string to_string(const Line<T> &l) {
    return l.to_string();
}

template <class T> struct Segment {
    Point<T> s, t;

    // 線分を通る直線に変換する。計算量 O(1)。
    // 整数幅: 座標のコピーと比較のみ。B=M、M=10^9ならint可。
    Line<T> toLine() const {
        return initLine(s, t);
    }

    // 線分の始点から終点へのベクトルを返す。計算量 O(1)。
    // 整数幅: 座標差を取るのでB=2M。M=10^9なら32ビット、int可。
    auto toVector() const {
        return t - s;
    }

    operator Line<T>() const {
        return toLine();
    }
};

// 異なる2点 (s, t) を結ぶ線分を初期化する。計算量 O(1)。
// 整数幅: 追加の座標演算なし。B=M、M=10^9ならint可。
template <class T> auto initSegment(Point<T> s, Point<T> t) {
    assert(s != t);
    return Segment<T>{s, t};
}

template <class T> auto toLine(const Segment<T> &s) {
    return s.toLine();
}

template <class T> auto toVector(const Segment<T> &s) {
    return s.toVector();
}
}

namespace std {
template <class T, class Id> struct hash<cplib::Point<T, Id>> {
    // 両座標のハッシュから点のハッシュを求める。計算量は座標のハッシュ計算に依存する。
    size_t operator()(cplib::Point<T> p) const requires (!cplib::isFloatingPoint<T>) {
        size_t a = hash<T>{}(p.x), b = hash<T>{}(p.y);
        return a ^ (b + 0x9e3779b97f4a7c15ULL + (a << 6) + (a >> 2));
    }
};
}

namespace std {
template <class T> struct hash<cplib::PointKey<T>> {
    // 各座標の厳密一致に整合するハッシュを返す。計算量は座標のハッシュ計算に依存する。
    size_t operator()(cplib::PointKey<T> p) const {
        size_t a = hash<T>{}(p.x), b = hash<T>{}(p.y);
        return a ^ (b + 0x9e3779b97f4a7c15ULL + (a << 6) + (a >> 2));
    }
};
}
