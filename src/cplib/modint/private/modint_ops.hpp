#pragma once
#include <cplib/common.hpp>
#include <ostream>

namespace cplib::detail {
// 共通の算術API。還元処理・内部表現は各版に持たせる。
template <class Derived> struct ModintOps {
    using cplib_modint_tag = void;

    // 剰余の和を返す。計算量 O(1)。
    friend constexpr Derived operator+(Derived a, Derived b) {
        return a += b;
    }

    // 剰余の差を返す。計算量 O(1)。
    friend constexpr Derived operator-(Derived a, Derived b) {
        return a -= b;
    }

    // 剰余の積を返す。計算量 O(1)。
    friend constexpr Derived operator*(Derived a, Derived b) {
        return a *= b;
    }

    // b の逆元を掛ける。計算量 O(log mod)。
    friend constexpr Derived operator/(Derived a, Derived b) {
        return a /= b;
    }

    // 非負整数乗を返す。計算量 O(log(n+1))。
    constexpr Derived pow(Int n) const {
        Derived result = 1, a = static_cast<const Derived &>(*this);
        while (n > 0) {
            if (n & 1)
                result *= a;
            a *= a;
            n >>= 1;
        }
        return result;
    }

    // 拡張ユークリッド法で逆元を返す。計算量 O(log mod)。
    constexpr Derived inv() const {
        Int x = static_cast<const Derived &>(*this).val(), y = Derived::mod(), u = 1, v = 0;
        assert(x != 0);
        while (y > 0) {
            Int t = x / y;
            x -= t * y;
            u -= t * v;
            std::swap(x, y);
            std::swap(u, v);
        }
        return Derived(u);
    }

    // b の逆元を掛けて代入する。計算量 O(log mod)。
    constexpr Derived &operator/=(Derived b) {
        auto &a = static_cast<Derived &>(*this);
        return a *= b.inv();
    }

    // 正規化された値をストリームに出力する。計算量 O(log mod)。
    friend std::ostream &operator<<(std::ostream &s, const Derived &x) {
        return s << x.val();
    }
};

template <std::integral I> constexpr UInt normalizeModint(I a, std::uint32_t m) {
    assert(m);
    if (UInt(a) < m)
        return UInt(a);
    if constexpr (std::is_unsigned_v<I>)
        return UInt(a) % m;
    else {
        Int r = Int(a) % Int(m);
        return UInt(r < 0 ? r + m : r);
    }
}
}

namespace cplib {
template <class T>
concept Modint = requires { typename T::cplib_modint_tag; };

// 整数から modint を作る。計算量 O(1)。
template <Modint M, class I> constexpr M init(I a) {
    return M(a);
}

// 正規化された剰余を返す。計算量 O(1)。
template <Modint M> constexpr Int val(M a) {
    return a.val();
}

// 符号なしの法を返す。計算量 O(1)。
template <Modint M> constexpr auto umod() {
    return M::umod();
}

// 符号なしの法を返す。計算量 O(1)。
template <Modint M> constexpr auto umod(M) {
    return M::umod();
}

// 符号付きの法を返す。計算量 O(1)。
template <Modint M> constexpr auto mod() {
    return M::mod();
}

// 符号付きの法を返す。計算量 O(1)。
template <Modint M> constexpr auto mod(M) {
    return M::mod();
}

// 法を返す。計算量 O(1)。
template <Modint M> constexpr auto get_M() {
    return M::umod();
}

// 還元用パラメータを返す。計算量 O(1)。
template <Modint M> constexpr auto get_param() {
    return M::get_param();
}

// 逆元を返す。計算量 O(log mod)。
template <Modint M> constexpr M inv(M a) {
    return a.inv();
}

// 非負整数乗を返す。計算量 O(log(n+1))。
template <Modint M> constexpr M pow(M a, Int n) {
    return a.pow(n);
}

// 剰余を十進文字列に変換する。計算量 O(log mod)。
template <Modint M> std::string to_string(M a) {
    return std::to_string(a.val());
}

// 動的法を設定する。計算量 O(1)。
template <Modint M, std::integral I> void setMod(I m) {
    M::setMod(m);
}
}
