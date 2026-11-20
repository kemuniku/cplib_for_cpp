#pragma once
#include <cplib/common.hpp>
#include <functional>
#include <sstream>

namespace cplib {
inline Int FRACTION_REDUCE_LIMIT = 1000000000;

// 以下の計算量は T の四則演算・比較を O(1)、値の絶対値の上限を V とする。
template <class T> struct Fraction {
    T num{}, den{};
    // 既定値を作る。計算量 O(1)。
    Fraction() = default;

    // 整数 n を分数にする。計算量 O(1)。
    Fraction(T n) : num(n), den(1) {
    }

    // 分子と分母から分数を作る。約分時 O(log(V+1))、省略時 O(1)。
    Fraction(T n, T d, bool do_reduce = true) : num(n), den(d) {
        if (do_reduce)
            reduce();
    }

    // 0/0かどうかを返す。O(1)。
    bool isNaN() const {
        return den == 0 && num == 0;
    }

    // 分子分母を約分し、分母の符号を正規化する。O(log max(|num|,|den|))。
    void reduce() {
        if (isNaN())
            return;
        T a = num < T(0) ? -num : num;
        T b = den < T(0) ? -den : den;
        while (b != T(0)) {
            T remainder = a % b;
            a = std::move(b);
            b = std::move(remainder);
        }
        num /= a;
        den /= a;
        if (den < 0) {
            den = -den;
            num = -num;
        }
    }

    // 閾値を超えた場合などに約分する。最悪 O(log(V+1))。
    void check_and_reduce() {
        if (den < 0 || num > FRACTION_REDUCE_LIMIT || den > FRACTION_REDUCE_LIMIT)
            reduce();
    }

    // 符号を正規化した逆数を返す。O(1)。
    Fraction inv() const {
        return num < 0 ? Fraction(-den, -num, false) : Fraction(den, num, false);
    }

    // 符号を反転した分数を返す。計算量 O(1)。
    Fraction operator-() const {
        return Fraction(-num, den, false);
    }

    // 分数を加算する。最悪 O(log(V+1))。
    Fraction &operator+=(const Fraction &y) {
        if (isNaN() || y.isNaN()) {
            num = den = 0;
            return *this;
        }
        if (den == 0 && y.den == 0) {
            if ((num > 0) != (y.num > 0))
                num = den = 0;
            return *this;
        }
        if (den == 0 || y.den == 0) {
            if (den != 0)
                *this = y;
            return *this;
        }
        num = num * y.den + y.num * den;
        den *= y.den;
        check_and_reduce();
        return *this;
    }

    // 分数を減算する。最悪 O(log(V+1))。
    Fraction &operator-=(const Fraction &y) {
        return *this += -y;
    }

    // 分数を乗算する。最悪 O(log(V+1))。
    Fraction &operator*=(const Fraction &y) {
        if (isNaN() || y.isNaN()) {
            num = den = 0;
            return *this;
        }
        num *= y.num;
        den *= y.den;
        check_and_reduce();
        return *this;
    }

    // 分数を除算する。最悪 O(log(V+1))。
    Fraction &operator/=(const Fraction &y) {
        if (isNaN() || y.isNaN()) {
            num = den = 0;
            return *this;
        }
        T n = num * y.den, d = den * y.num;
        num = n;
        den = d;
        check_and_reduce();
        return *this;
    }

    // 分数の和を返す。最悪 O(log(V+1))。
    friend Fraction operator+(Fraction x, const Fraction &y) {
        return x += y;
    }

    // 分数の差を返す。最悪 O(log(V+1))。
    friend Fraction operator-(Fraction x, const Fraction &y) {
        return x -= y;
    }

    // 分数の積を返す。最悪 O(log(V+1))。
    friend Fraction operator*(Fraction x, const Fraction &y) {
        return x *= y;
    }

    // 分数の商を返す。最悪 O(log(V+1))。
    friend Fraction operator/(Fraction x, const Fraction &y) {
        return x /= y;
    }

    // 左辺が右辺より小さいか判定する。計算量 O(1)。
    friend bool operator<(const Fraction &x, const Fraction &y) {
        if (x.isNaN() || y.isNaN())
            return false;
        if (x.den == 0 && y.den == 0)
            return x.num < y.num;
        return x.num * y.den < y.num * x.den;
    }

    // 左辺が右辺より大きいか判定する。計算量 O(1)。
    friend bool operator>(const Fraction &x, const Fraction &y) {
        return y < x;
    }

    // 分数が等しいか判定する。計算量 O(1)。
    friend bool operator==(const Fraction &x, const Fraction &y) {
        if (x.isNaN() || y.isNaN())
            return false;
        if (x.den == 0 && y.den == 0)
            return (x.num > 0) == (y.num > 0);
        return x.num * y.den == y.num * x.den;
    }

    // 左辺が右辺以下か判定する。計算量 O(1)。
    friend bool operator<=(const Fraction &x, const Fraction &y) {
        return !(x > y);
    }

    // 左辺が右辺以上か判定する。計算量 O(1)。
    friend bool operator>=(const Fraction &x, const Fraction &y) {
        return !(x < y);
    }

    // 約分した分数をストリームに出力する。約分 O(log(V+1)) と出力文字数に比例する時間。
    friend std::ostream &operator<<(std::ostream &out, Fraction x) {
        x.reduce();
        return out << x.num << '/' << x.den;
    }
};

// 分子と分母から分数を作る。計算量は Fraction(n, d, reduce) と同じ。
template <class T> Fraction<T> initFraction(T n, T d, bool reduce = true) {
    return {n, d, reduce};
}

// 整数を分数にする。計算量 O(1)。
template <class T> Fraction<T> initFraction(T n) {
    return {n};
}

// 0/0 かどうかを判定する。計算量 O(1)。
template <class T> bool isNaN(const Fraction<T> &x) {
    return x.isNaN();
}

// 分数を約分する。計算量は x.reduce() と同じ。
template <class T> void reduce(Fraction<T> &x) {
    x.reduce();
}

// 分数の逆数を返す。計算量 O(1)。
template <class T> Fraction<T> inv(const Fraction<T> &x) {
    return x.inv();
}

// 分数の絶対値を返す。計算量 O(1)。
template <class T> Fraction<T> abs(Fraction<T> x) {
    if (x.num < 0)
        x.num = -x.num;
    return x;
}

// 分数同士を比較し、小さい場合は -1、等しい場合は 0、大きい場合は 1 を返す。計算量 O(1)。
template <class T> int cmp(const Fraction<T> &x, const Fraction<T> &y) {
    return x < y ? -1 : x == y ? 0 : 1;
}

// 分数と整数を比較し、小さい場合は -1、等しい場合は 0、大きい場合は 1 を返す。計算量 O(1)。
template <class T> int cmp(const Fraction<T> &x, T y) {
    return cmp(x, Fraction<T>(y));
}

// 整数と分数を比較し、小さい場合は -1、等しい場合は 0、大きい場合は 1 を返す。計算量 O(1)。
template <class T> int cmp(T x, const Fraction<T> &y) {
    return cmp(Fraction<T>(x), y);
}

// 分数を double に変換する。計算量 O(1)。
template <class T> double toFloat(const Fraction<T> &x) {
    return static_cast<double>(x.num) / static_cast<double>(x.den);
}

// 約分した分数の文字列を返す。計算量はストリーム出力と同じ。
template <class T> std::string to_string(const Fraction<T> &x) {
    std::ostringstream out;
    out << x;
    return out.str();
}

// 分数の非負整数乗を求める。O(log n)回の分数演算。
template <class T> Fraction<T> pow(Fraction<T> x, Int n) {
    Fraction<T> result(T(1));
    while (n > 0) {
        if (n & 1)
            result *= x;
        x *= x;
        n >>= 1;
    }
    return result;
}

// 約分して分数のハッシュ値を返す。計算量は約分と T のハッシュ計算の和。
template <class T> std::size_t hashValue(Fraction<T> x) {
    x.reduce();
    auto a = std::hash<T>{}(x.num), b = std::hash<T>{}(x.den);
    return a ^ (b + 0x9e3779b97f4a7c15ULL + (a << 6) + (a >> 2));
}
}

namespace std {
template <class T> struct hash<cplib::Fraction<T>> {
    // 約分した分数のハッシュ値を返す。計算量は cplib::hashValue(x) と同じ。
    size_t operator()(const cplib::Fraction<T> &x) const {
        return cplib::hashValue(x);
    }
};
}
