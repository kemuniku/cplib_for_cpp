#pragma once
#include <cplib/common.hpp>
#include <functional>
#include <sstream>

namespace cplib {
inline Int FRACTION_REDUCE_LIMIT = 1000000000;

template <class T> struct Fraction {
    T num{}, den{};
    Fraction() = default;

    Fraction(T n) : num(n), den(1) {
    }

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
        using std::abs;
        auto g = std::gcd(abs(num), abs(den));
        num /= g;
        den /= g;
        if (den < 0) {
            den = -den;
            num = -num;
        }
    }

    void check_and_reduce() {
        if (den < 0 || num > FRACTION_REDUCE_LIMIT || den > FRACTION_REDUCE_LIMIT)
            reduce();
    }

    // 符号を正規化した逆数を返す。O(1)。
    Fraction inv() const {
        return num < 0 ? Fraction(-den, -num, false) : Fraction(den, num, false);
    }

    Fraction operator-() const {
        return Fraction(-num, den, false);
    }

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

    Fraction &operator-=(const Fraction &y) {
        return *this += -y;
    }

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

    friend Fraction operator+(Fraction x, const Fraction &y) {
        return x += y;
    }

    friend Fraction operator-(Fraction x, const Fraction &y) {
        return x -= y;
    }

    friend Fraction operator*(Fraction x, const Fraction &y) {
        return x *= y;
    }

    friend Fraction operator/(Fraction x, const Fraction &y) {
        return x /= y;
    }

    friend bool operator<(const Fraction &x, const Fraction &y) {
        if (x.isNaN() || y.isNaN())
            return false;
        if (x.den == 0 && y.den == 0)
            return x.num < y.num;
        return x.num * y.den < y.num * x.den;
    }

    friend bool operator>(const Fraction &x, const Fraction &y) {
        return y < x;
    }

    friend bool operator==(const Fraction &x, const Fraction &y) {
        if (x.isNaN() || y.isNaN())
            return false;
        if (x.den == 0 && y.den == 0)
            return (x.num > 0) == (y.num > 0);
        return x.num * y.den == y.num * x.den;
    }

    friend bool operator<=(const Fraction &x, const Fraction &y) {
        return !(x > y);
    }

    friend bool operator>=(const Fraction &x, const Fraction &y) {
        return !(x < y);
    }

    // 約分した分数を文字列にする。
    friend std::ostream &operator<<(std::ostream &out, Fraction x) {
        x.reduce();
        return out << x.num << '/' << x.den;
    }
};

template <class T> Fraction<T> initFraction(T n, T d, bool reduce = true) {
    return {n, d, reduce};
}

template <class T> Fraction<T> initFraction(T n) {
    return {n};
}

template <class T> bool isNaN(const Fraction<T> &x) {
    return x.isNaN();
}

template <class T> void reduce(Fraction<T> &x) {
    x.reduce();
}

template <class T> Fraction<T> inv(const Fraction<T> &x) {
    return x.inv();
}

template <class T> Fraction<T> abs(Fraction<T> x) {
    if (x.num < 0)
        x.num = -x.num;
    return x;
}

// 分数同士を比較し、小さい場合は -1、等しい場合は 0、大きい場合は 1 を返す。
template <class T> int cmp(const Fraction<T> &x, const Fraction<T> &y) {
    return x < y ? -1 : x == y ? 0 : 1;
}

// 分数と整数を比較し、小さい場合は -1、等しい場合は 0、大きい場合は 1 を返す。
template <class T> int cmp(const Fraction<T> &x, T y) {
    return cmp(x, Fraction<T>(y));
}

// 整数と分数を比較し、小さい場合は -1、等しい場合は 0、大きい場合は 1 を返す。
template <class T> int cmp(T x, const Fraction<T> &y) {
    return cmp(Fraction<T>(x), y);
}

template <class T> double toFloat(const Fraction<T> &x) {
    return static_cast<double>(x.num) / static_cast<double>(x.den);
}

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

template <class T> std::size_t hash(Fraction<T> x) {
    x.reduce();
    auto a = std::hash<T>{}(x.num), b = std::hash<T>{}(x.den);
    return a ^ (b + 0x9e3779b97f4a7c15ULL + (a << 6) + (a >> 2));
}
}

namespace std {
template <class T> struct hash<cplib::Fraction<T>> {
    size_t operator()(const cplib::Fraction<T> &x) const {
        return cplib::hash(x);
    }
};
}
