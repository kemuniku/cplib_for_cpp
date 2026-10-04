#pragma once
#include <cplib/matrix/field_matrix_ops.hpp>
#include <sstream>

namespace cplib {
// 固定長版は連続したarray<H*W,T>を保持し、寸法を型に含める。
template <Int H, Int W, class T> class StaticMatrix {
    static_assert(H >= 0 && W >= 0);
    std::array<T, H * W> arr_{};

public:
    using value_type = T;
    StaticMatrix() = default;

    explicit StaticMatrix(T value) {
        arr_.fill(value);
    }

    explicit StaticMatrix(const std::array<std::array<T, W>, H> &rows) {
        for (Int i = 0; i < H; ++i)
            for (Int j = 0; j < W; ++j)
                (*this)(i, j) = rows[i][j];
    }

    constexpr Int h() const {
        return H;
    }

    constexpr Int w() const {
        return W;
    }

    T &operator()(Int row, Int col) {
        return arr_[row * W + col];
    }

    const T &operator()(Int row, Int col) const {
        return arr_[row * W + col];
    }

    struct Row {
        StaticMatrix *owner;
        Int row;

        T &operator[](Int col) {
            return (*owner)(row, col);
        }

        operator std::array<T, W>() const {
            std::array<T, W> out;
            for (Int j = 0; j < W; ++j)
                out[j] = (*owner)(row, j);
            return out;
        }

        Row &operator=(const std::array<T, W> &values) {
            for (Int j = 0; j < W; ++j)
                (*owner)(row, j) = values[j];
            return *this;
        }

        Row &operator=(const Row &other) {
            return *this = std::array<T, W>(other);
        }
    };

    Row operator[](Int row) {
        return {this, row};
    }

    std::array<T, W> operator[](Int row) const {
        std::array<T, W> out;
        for (Int j = 0; j < W; ++j)
            out[j] = (*this)(row, j);
        return out;
    }

    bool operator==(const StaticMatrix &) const = default;

    StaticMatrix operator-() const {
        StaticMatrix out;
        for (Int i = 0; i < H * W; ++i)
            out.arr_[i] = -arr_[i];
        return out;
    }

    template <Int C> auto operator*(const StaticMatrix<W, C, T> &b) const {
        StaticMatrix<H, C, T> out;
        for (Int i = 0; i < H; ++i)
            for (Int j = 0; j < C; ++j)
                for (Int k = 0; k < W; ++k)
                    out(i, j) += (*this)(i, k) * b(k, j);
        return out;
    }

    StaticMatrix &operator*=(const StaticMatrix<W, W, T> &b) {
        *this = (*this) * b;
        return *this;
    }

    StaticMatrix &operator*=(T value) {
        for (T &x : arr_)
            x *= value;
        return *this;
    }

    friend StaticMatrix operator*(StaticMatrix a, T b) {
        a *= b;
        return a;
    }

    friend StaticMatrix operator*(T a, StaticMatrix b) {
        b *= a;
        return b;
    }

    StaticMatrix &operator+=(const StaticMatrix &b) {
        assert(h() == b.h() && w() == b.w());
        for (Int i = 0; i < h(); ++i)
            for (Int j = 0; j < w(); ++j) {
                T a = (*this)(i, j);
                T rhs = b(i, j);
                (*this)(i, j) = a + rhs;
            }
        return *this;
    }

    StaticMatrix &operator+=(T b) {
        for (Int i = 0; i < h(); ++i)
            for (Int j = 0; j < w(); ++j) {
                T a = (*this)(i, j);
                (*this)(i, j) = a + b;
            }
        return *this;
    }

    friend StaticMatrix operator+(StaticMatrix a, const StaticMatrix &b) {
        a += b;
        return a;
    }

    friend StaticMatrix operator+(StaticMatrix a, T b) {
        a += b;
        return a;
    }

    friend StaticMatrix operator+(T a, StaticMatrix matrix) {
        for (Int i = 0; i < matrix.h(); ++i)
            for (Int j = 0; j < matrix.w(); ++j) {
                T b = matrix(i, j);
                matrix(i, j) = a + b;
            }
        return matrix;
    }

    StaticMatrix &operator-=(const StaticMatrix &b) {
        assert(h() == b.h() && w() == b.w());
        for (Int i = 0; i < h(); ++i)
            for (Int j = 0; j < w(); ++j) {
                T a = (*this)(i, j);
                T rhs = b(i, j);
                (*this)(i, j) = a - rhs;
            }
        return *this;
    }

    StaticMatrix &operator-=(T b) {
        for (Int i = 0; i < h(); ++i)
            for (Int j = 0; j < w(); ++j) {
                T a = (*this)(i, j);
                (*this)(i, j) = a - b;
            }
        return *this;
    }

    friend StaticMatrix operator-(StaticMatrix a, const StaticMatrix &b) {
        a -= b;
        return a;
    }

    friend StaticMatrix operator-(StaticMatrix a, T b) {
        a -= b;
        return a;
    }

    friend StaticMatrix operator-(T a, StaticMatrix matrix) {
        for (Int i = 0; i < matrix.h(); ++i)
            for (Int j = 0; j < matrix.w(); ++j) {
                T b = matrix(i, j);
                matrix(i, j) = a - b;
            }
        return matrix;
    }

    StaticMatrix &operator&=(const StaticMatrix &b)
        requires(std::integral<T>)
    {
        assert(h() == b.h() && w() == b.w());
        for (Int i = 0; i < h(); ++i)
            for (Int j = 0; j < w(); ++j) {
                T a = (*this)(i, j);
                T rhs = b(i, j);
                (*this)(i, j) = a & rhs;
            }
        return *this;
    }

    StaticMatrix &operator&=(T b)
        requires(std::integral<T>)
    {
        for (Int i = 0; i < h(); ++i)
            for (Int j = 0; j < w(); ++j) {
                T a = (*this)(i, j);
                (*this)(i, j) = a & b;
            }
        return *this;
    }

    friend StaticMatrix operator&(StaticMatrix a, const StaticMatrix &b)
        requires(std::integral<T>)
    {
        a &= b;
        return a;
    }

    friend StaticMatrix operator&(StaticMatrix a, T b)
        requires(std::integral<T>)
    {
        a &= b;
        return a;
    }

    friend StaticMatrix operator&(T a, StaticMatrix matrix)
        requires(std::integral<T>)
    {
        for (Int i = 0; i < matrix.h(); ++i)
            for (Int j = 0; j < matrix.w(); ++j) {
                T b = matrix(i, j);
                matrix(i, j) = a & b;
            }
        return matrix;
    }

    StaticMatrix &operator|=(const StaticMatrix &b)
        requires(std::integral<T>)
    {
        assert(h() == b.h() && w() == b.w());
        for (Int i = 0; i < h(); ++i)
            for (Int j = 0; j < w(); ++j) {
                T a = (*this)(i, j);
                T rhs = b(i, j);
                (*this)(i, j) = a | rhs;
            }
        return *this;
    }

    StaticMatrix &operator|=(T b)
        requires(std::integral<T>)
    {
        for (Int i = 0; i < h(); ++i)
            for (Int j = 0; j < w(); ++j) {
                T a = (*this)(i, j);
                (*this)(i, j) = a | b;
            }
        return *this;
    }

    friend StaticMatrix operator|(StaticMatrix a, const StaticMatrix &b)
        requires(std::integral<T>)
    {
        a |= b;
        return a;
    }

    friend StaticMatrix operator|(StaticMatrix a, T b)
        requires(std::integral<T>)
    {
        a |= b;
        return a;
    }

    friend StaticMatrix operator|(T a, StaticMatrix matrix)
        requires(std::integral<T>)
    {
        for (Int i = 0; i < matrix.h(); ++i)
            for (Int j = 0; j < matrix.w(); ++j) {
                T b = matrix(i, j);
                matrix(i, j) = a | b;
            }
        return matrix;
    }

    StaticMatrix &operator^=(const StaticMatrix &b)
        requires(std::integral<T>)
    {
        assert(h() == b.h() && w() == b.w());
        for (Int i = 0; i < h(); ++i)
            for (Int j = 0; j < w(); ++j) {
                T a = (*this)(i, j);
                T rhs = b(i, j);
                (*this)(i, j) = a ^ rhs;
            }
        return *this;
    }

    StaticMatrix &operator^=(T b)
        requires(std::integral<T>)
    {
        for (Int i = 0; i < h(); ++i)
            for (Int j = 0; j < w(); ++j) {
                T a = (*this)(i, j);
                (*this)(i, j) = a ^ b;
            }
        return *this;
    }

    friend StaticMatrix operator^(StaticMatrix a, const StaticMatrix &b)
        requires(std::integral<T>)
    {
        a ^= b;
        return a;
    }

    friend StaticMatrix operator^(StaticMatrix a, T b)
        requires(std::integral<T>)
    {
        a ^= b;
        return a;
    }

    friend StaticMatrix operator^(T a, StaticMatrix matrix)
        requires(std::integral<T>)
    {
        for (Int i = 0; i < matrix.h(); ++i)
            for (Int j = 0; j < matrix.w(); ++j) {
                T b = matrix(i, j);
                matrix(i, j) = a ^ b;
            }
        return matrix;
    }

    StaticMatrix &operator<<=(const StaticMatrix &b)
        requires(std::integral<T>)
    {
        assert(h() == b.h() && w() == b.w());
        for (Int i = 0; i < h(); ++i)
            for (Int j = 0; j < w(); ++j) {
                T a = (*this)(i, j);
                T rhs = b(i, j);
                (*this)(i, j) = T(std::make_unsigned_t<T>(a) << rhs);
            }
        return *this;
    }

    StaticMatrix &operator<<=(T b)
        requires(std::integral<T>)
    {
        for (Int i = 0; i < h(); ++i)
            for (Int j = 0; j < w(); ++j) {
                T a = (*this)(i, j);
                (*this)(i, j) = T(std::make_unsigned_t<T>(a) << b);
            }
        return *this;
    }

    friend StaticMatrix operator<<(StaticMatrix a, const StaticMatrix &b)
        requires(std::integral<T>)
    {
        a <<= b;
        return a;
    }

    friend StaticMatrix operator<<(StaticMatrix a, T b)
        requires(std::integral<T>)
    {
        a <<= b;
        return a;
    }

    friend StaticMatrix operator<<(T a, StaticMatrix matrix)
        requires(std::integral<T>)
    {
        for (Int i = 0; i < matrix.h(); ++i)
            for (Int j = 0; j < matrix.w(); ++j) {
                T b = matrix(i, j);
                matrix(i, j) = T(std::make_unsigned_t<T>(a) << b);
            }
        return matrix;
    }

    StaticMatrix &operator>>=(const StaticMatrix &b)
        requires(std::integral<T>)
    {
        assert(h() == b.h() && w() == b.w());
        for (Int i = 0; i < h(); ++i)
            for (Int j = 0; j < w(); ++j) {
                T a = (*this)(i, j);
                T rhs = b(i, j);
                (*this)(i, j) = a >> rhs;
            }
        return *this;
    }

    StaticMatrix &operator>>=(T b)
        requires(std::integral<T>)
    {
        for (Int i = 0; i < h(); ++i)
            for (Int j = 0; j < w(); ++j) {
                T a = (*this)(i, j);
                (*this)(i, j) = a >> b;
            }
        return *this;
    }

    friend StaticMatrix operator>>(StaticMatrix a, const StaticMatrix &b)
        requires(std::integral<T>)
    {
        a >>= b;
        return a;
    }

    friend StaticMatrix operator>>(StaticMatrix a, T b)
        requires(std::integral<T>)
    {
        a >>= b;
        return a;
    }

    friend StaticMatrix operator>>(T a, StaticMatrix matrix)
        requires(std::integral<T>)
    {
        for (Int i = 0; i < matrix.h(); ++i)
            for (Int j = 0; j < matrix.w(); ++j) {
                T b = matrix(i, j);
                matrix(i, j) = a >> b;
            }
        return matrix;
    }

    StaticMatrix &operator/=(const StaticMatrix &b)
        requires(std::integral<T>)
    {
        assert(h() == b.h() && w() == b.w());
        for (Int i = 0; i < h(); ++i)
            for (Int j = 0; j < w(); ++j) {
                T a = (*this)(i, j);
                T rhs = b(i, j);
                (*this)(i, j) = a / rhs;
            }
        return *this;
    }

    StaticMatrix &operator/=(T b)
        requires(std::integral<T>)
    {
        for (Int i = 0; i < h(); ++i)
            for (Int j = 0; j < w(); ++j) {
                T a = (*this)(i, j);
                (*this)(i, j) = a / b;
            }
        return *this;
    }

    friend StaticMatrix operator/(StaticMatrix a, const StaticMatrix &b)
        requires(std::integral<T>)
    {
        a /= b;
        return a;
    }

    friend StaticMatrix operator/(StaticMatrix a, T b)
        requires(std::integral<T>)
    {
        a /= b;
        return a;
    }

    friend StaticMatrix operator/(T a, StaticMatrix matrix)
        requires(std::integral<T>)
    {
        for (Int i = 0; i < matrix.h(); ++i)
            for (Int j = 0; j < matrix.w(); ++j) {
                T b = matrix(i, j);
                matrix(i, j) = a / b;
            }
        return matrix;
    }

    StaticMatrix &operator%=(const StaticMatrix &b)
        requires(std::integral<T>)
    {
        assert(h() == b.h() && w() == b.w());
        for (Int i = 0; i < h(); ++i)
            for (Int j = 0; j < w(); ++j) {
                T a = (*this)(i, j);
                T rhs = b(i, j);
                (*this)(i, j) = a % rhs;
            }
        return *this;
    }

    StaticMatrix &operator%=(T b)
        requires(std::integral<T>)
    {
        for (Int i = 0; i < h(); ++i)
            for (Int j = 0; j < w(); ++j) {
                T a = (*this)(i, j);
                (*this)(i, j) = a % b;
            }
        return *this;
    }

    friend StaticMatrix operator%(StaticMatrix a, const StaticMatrix &b)
        requires(std::integral<T>)
    {
        a %= b;
        return a;
    }

    friend StaticMatrix operator%(StaticMatrix a, T b)
        requires(std::integral<T>)
    {
        a %= b;
        return a;
    }

    friend StaticMatrix operator%(T a, StaticMatrix matrix)
        requires(std::integral<T>)
    {
        for (Int i = 0; i < matrix.h(); ++i)
            for (Int j = 0; j < matrix.w(); ++j) {
                T b = matrix(i, j);
                matrix(i, j) = a % b;
            }
        return matrix;
    }

    static StaticMatrix identity(Int n = H, T one = T(1), T zero = T(0)) {
        assert(H == W && n == H);
        StaticMatrix out(zero);
        for (Int i = 0; i < H; ++i)
            out(i, i) = one;
        return out;
    }

    StaticMatrix pow(Int exponent) const
        requires(H == W)
    {
        StaticMatrix result = identity(), base = *this;
        while (exponent > 0) {
            if (exponent & 1)
                result *= base;
            base *= base;
            exponent >>= 1;
        }
        return result;
    }

    T sum() const {
        T result = T(0);
        for (T x : arr_)
            result += x;
        return result;
    }

    std::string str() const {
        std::ostringstream out;
        for (Int i = 0; i < H; ++i) {
            if (i)
                out << '\n';
            for (Int j = 0; j < W; ++j) {
                if (j)
                    out << ' ';
                if constexpr (std::is_same_v<T, bool>)
                    out << ((*this)(i, j) ? "true" : "false");
                else
                    out << (*this)(i, j);
            }
        }
        return out.str();
    }

    std::size_t hash() const {
        std::size_t seed = H * W;
        for (T value : arr_)
            seed ^= std::hash<T>{}(value) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
        return seed;
    }

    // 左上height行width列の階数を求める。O(h*w*min(h,w))。
    Int rank(Int height = H, Int width = W) const {
        return fieldRank(matrixRows(*this, height, width), width);
    }

    // 左上n×nの行列式を求める。空行列は1。O(n^3)。
    T determinant(Int n = H) const {
        assert(0 <= n && n <= std::min(H, W));
        return fieldDeterminant(matrixRows(*this, n, n));
    }

    // 対称な左上n×n（nは偶数）のhafnianを求める。O(n^2*2^(n/2))。
    T hafnian(Int n = H) const {
        assert(0 <= n && n <= std::min(H, W));
        return fieldHafnian(matrixRows(*this, n, n));
    }

    // 左上height行width列でAx=bの特殊解と核の基底を返す。解なしはstd::nullopt。
    template <class B> auto solveLinearSystem(const B &b, Int height = H, Int width = W) const {
        return fieldSolve(matrixRows(*this, height, width), width, b);
    }

    // 左上n×nの逆行列を返す。範囲外は零、特異行列はstd::nullopt。O(n^3)。
    std::optional<StaticMatrix> inverse(Int n = H) const {
        assert(0 <= n && n <= std::min(H, W));
        auto rows = fieldAdjugateInverse(matrixRows(*this, n, n), false);
        if (!rows)
            return std::nullopt;
        StaticMatrix out;
        for (Int i = 0; i < n; ++i)
            for (Int j = 0; j < n; ++j)
                out(i, j) = (*rows)[i][j];
        return out;
    }

    // 左上n×nの余因子行列を返す。範囲外は零。O(n^3)。
    StaticMatrix adjugate(Int n = H) const {
        assert(0 <= n && n <= std::min(H, W));
        auto rows = fieldAdjugateInverse(matrixRows(*this, n, n), true);
        StaticMatrix out;
        for (Int i = 0; i < n; ++i)
            for (Int j = 0; j < n; ++j)
                out(i, j) = (*rows)[i][j];
        return out;
    }
};

template <Int H, Int W, class T> auto initMatrix(T value) {
    return StaticMatrix<H, W, T>(value);
}

template <std::size_t H, std::size_t W, class T>
auto initMatrix(const std::array<std::array<T, W>, H> &rows) {
    return StaticMatrix<H, W, T>(rows);
}

template <std::size_t H, std::size_t W, class T>
auto toMatrix(const std::array<std::array<T, W>, H> &rows) {
    return initMatrix(rows);
}

template <Int H, Int W, class T> auto identity_matrix(Int n = H, T one = T(1), T zero = T(0)) {
    return StaticMatrix<H, W, T>::identity(n, one, zero);
}

template <Int H, Int W, class T> Int h(const StaticMatrix<H, W, T> &m) {
    return m.h();
}

template <Int H, Int W, class T> Int w(const StaticMatrix<H, W, T> &m) {
    return m.w();
}

template <Int H, Int W, class T> auto hash(const StaticMatrix<H, W, T> &m) {
    return m.hash();
}

template <Int H, Int W, class T> auto pow(const StaticMatrix<H, W, T> &m, Int n) {
    return m.pow(n);
}

template <Int H, Int W, class T> T sum(const StaticMatrix<H, W, T> &m) {
    return m.sum();
}

template <Int H, Int W, class T>
Int rank(const StaticMatrix<H, W, T> &m, Int height = H, Int width = W) {
    return m.rank(height, width);
}

template <Int H, Int W, class T> T determinant(const StaticMatrix<H, W, T> &m, Int n = H) {
    return m.determinant(n);
}

template <Int H, Int W, class T> T hafnian(const StaticMatrix<H, W, T> &m, Int n = H) {
    return m.hafnian(n);
}

template <Int H, Int W, class T> auto inverse(const StaticMatrix<H, W, T> &m, Int n = H) {
    return m.inverse(n);
}

template <Int H, Int W, class T> auto adjugate(const StaticMatrix<H, W, T> &m, Int n = H) {
    return m.adjugate(n);
}

template <Int H, Int W, class T> auto to_string(const StaticMatrix<H, W, T> &m) {
    return m.str();
}

template <Int H, Int W, class T, class B>
auto solveLinearSystem(const StaticMatrix<H, W, T> &m, const B &b, Int height = H, Int width = W) {
    return m.solveLinearSystem(b, height, width);
}
}

namespace std {
template <cplib::Int H, cplib::Int W, class T> struct hash<cplib::StaticMatrix<H, W, T>> {
    size_t operator()(const cplib::StaticMatrix<H, W, T> &m) const {
        return m.hash();
    }
};
}
