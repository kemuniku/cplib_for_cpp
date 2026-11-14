#pragma once
#include <cplib/matrix/field_matrix_ops.hpp>
#include <sstream>

namespace cplib {
// 固定長版は連続したarray<H*W,T>を保持し、寸法を型に含める。
// 要素演算をO(1)として、要素ごとの算術・ビット演算と比較は O(HW+1)。
// 行列積 H×W と W×C は O(HWC+HC+1)。寸法・要素取得と非constの行ビュー取得は O(1)。
// constの行取得、行の配列への変換・行への代入は O(W+1)。
template <Int H, Int W, class T> class StaticMatrix {
    static_assert(H >= 0 && W >= 0);
    std::array<T, H * W> arr_{};

public:
    using value_type = T;
    // 零初期化した行列を作る。計算量 O(HW+1)。
    StaticMatrix() = default;

    // 指定値で行列を初期化する。計算量 O(HW+1)。
    explicit StaticMatrix(T value) {
        arr_.fill(value);
    }

    // 配列をコピーして行列を作る。計算量 O(HW+1)。
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
            std::array<T, W> output_value;
            for (Int j = 0; j < W; ++j)
                output_value[j] = (*owner)(row, j);
            return output_value;
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
        std::array<T, W> output_value;
        for (Int j = 0; j < W; ++j)
            output_value[j] = (*this)(row, j);
        return output_value;
    }

    bool operator==(const StaticMatrix &) const = default;

    StaticMatrix operator-() const {
        StaticMatrix output_value;
        for (Int i = 0; i < H * W; ++i)
            output_value.arr_[i] = -arr_[i];
        return output_value;
    }

    template <Int C> auto operator*(const StaticMatrix<W, C, T> &b) const {
        StaticMatrix<H, C, T> output_value;
        for (Int i = 0; i < H; ++i)
            for (Int j = 0; j < C; ++j)
                for (Int k = 0; k < W; ++k)
                    output_value(i, j) += (*this)(i, k) * b(k, j);
        return output_value;
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

    // 正方単位行列を作る。計算量 O(H²+1)。
    static StaticMatrix identity(Int n = H, T one = T(1), T zero = T(0)) {
        assert(H == W && n == H);
        StaticMatrix output_value(zero);
        for (Int i = 0; i < H; ++i)
            output_value(i, i) = one;
        return output_value;
    }

    // H次正方行列を非負のk乗にする。計算量 O(H³ log(k+1)+H²+1)。
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

    // 全要素の和を返す。計算量 O(HW+1)。
    T sum() const {
        T result = T(0);
        for (T x : arr_)
            result += x;
        return result;
    }

    // 文字列に変換する。要素の出力時間の合計と O(HW+H+1) が必要。
    std::string str() const {
        std::ostringstream output_value;
        for (Int i = 0; i < H; ++i) {
            if (i)
                output_value << '\n';
            for (Int j = 0; j < W; ++j) {
                if (j)
                    output_value << ' ';
                if constexpr (std::is_same_v<T, bool>)
                    output_value << ((*this)(i, j) ? "true" : "false");
                else
                    output_value << (*this)(i, j);
            }
        }
        return output_value.str();
    }

    // ハッシュを返す。要素のハッシュを O(1) として O(HW+1)。
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
    // height=h、width=wとして計算量 O((h+1)(w+1)²)。
    template <class B> auto solveLinearSystem(const B &b, Int height = H, Int width = W) const {
        return fieldSolve(matrixRows(*this, height, width), width, b);
    }

    // 左上n×nの逆行列を返す。範囲外は零、特異行列はstd::nullopt。O(n^3)。
    std::optional<StaticMatrix> inverse(Int n = H) const {
        assert(0 <= n && n <= std::min(H, W));
        auto rows = fieldAdjugateInverse(matrixRows(*this, n, n), false);
        if (!rows)
            return std::nullopt;
        StaticMatrix output_value;
        for (Int i = 0; i < n; ++i)
            for (Int j = 0; j < n; ++j)
                output_value(i, j) = (*rows)[i][j];
        return output_value;
    }

    // 左上n×nの余因子行列を返す。範囲外は零。O(n^3)。
    StaticMatrix adjugate(Int n = H) const {
        assert(0 <= n && n <= std::min(H, W));
        auto rows = fieldAdjugateInverse(matrixRows(*this, n, n), true);
        StaticMatrix output_value;
        for (Int i = 0; i < n; ++i)
            for (Int j = 0; j < n; ++j)
                output_value(i, j) = (*rows)[i][j];
        return output_value;
    }
};

// 指定値でH行W列を初期化する。計算量 O(HW+1)。
template <Int H, Int W, class T> auto initMatrix(T value) {
    return StaticMatrix<H, W, T>(value);
}

// 配列をコピーしてH行W列の行列を作る。計算量 O(HW+1)。
template <std::size_t H, std::size_t W, class T>
auto initMatrix(const std::array<std::array<T, W>, H> &rows) {
    return StaticMatrix<H, W, T>(rows);
}

// 配列をコピーしてH行W列の行列を作る。計算量 O(HW+1)。
template <std::size_t H, std::size_t W, class T>
auto toMatrix(const std::array<std::array<T, W>, H> &rows) {
    return initMatrix(rows);
}

// H次の単位行列を返す。計算量 O(H²+1)。
template <Int H, Int W, class T> auto identity_matrix(Int n = H, T one = T(1), T zero = T(0)) {
    return StaticMatrix<H, W, T>::identity(n, one, zero);
}

// 行数を返す。計算量 O(1)。
template <Int H, Int W, class T> Int h(const StaticMatrix<H, W, T> &m) {
    return m.h();
}

// 列数を返す。計算量 O(1)。
template <Int H, Int W, class T> Int w(const StaticMatrix<H, W, T> &m) {
    return m.w();
}

// 要素のハッシュを O(1) として行列のハッシュを O(HW+1) で求める。
template <Int H, Int W, class T> auto hashValue(const StaticMatrix<H, W, T> &m) {
    return m.hash();
}

// H次正方行列を非負のk乗にする。計算量 O(H³ log(k+1)+H²+1)。
template <Int H, Int W, class T> auto pow(const StaticMatrix<H, W, T> &m, Int n) {
    return m.pow(n);
}

// 全要素の和を返す。計算量 O(HW+1)。
template <Int H, Int W, class T> T sum(const StaticMatrix<H, W, T> &m) {
    return m.sum();
}

// 左上h行w列の階数を返す。計算量 O(hw min(h,w)+h+w+1)。
template <Int H, Int W, class T>
Int rank(const StaticMatrix<H, W, T> &m, Int height = H, Int width = W) {
    return m.rank(height, width);
}

// 左上n次正方行列の行列式を返す。計算量 O(n³+n+1)。
template <Int H, Int W, class T> T determinant(const StaticMatrix<H, W, T> &m, Int n = H) {
    return m.determinant(n);
}

// 左上n次の対称な偶数次行列のhafnianを返す。計算量 O(n² 2^(n/2)+1)。
template <Int H, Int W, class T> T hafnian(const StaticMatrix<H, W, T> &m, Int n = H) {
    return m.hafnian(n);
}

// 左上n次正方行列の逆行列を返す。計算量 O(n³+HW+1)。
template <Int H, Int W, class T> auto inverse(const StaticMatrix<H, W, T> &m, Int n = H) {
    return m.inverse(n);
}

// 左上n次正方行列の余因子行列を返す。計算量 O(n³+HW+1)。
template <Int H, Int W, class T> auto adjugate(const StaticMatrix<H, W, T> &m, Int n = H) {
    return m.adjugate(n);
}

// 行列を文字列に変換する。要素の出力時間の合計と O(HW+H+1) が必要。
template <Int H, Int W, class T> auto to_string(const StaticMatrix<H, W, T> &m) {
    return m.str();
}

// 左上h行w列を用いたAx=bの特殊解と核の基底を返す。計算量 O((h+1)(w+1)²)。
template <Int H, Int W, class T, class B>
auto solveLinearSystem(const StaticMatrix<H, W, T> &m, const B &b, Int height = H, Int width = W) {
    return m.solveLinearSystem(b, height, width);
}
}

namespace std {
template <cplib::Int H, cplib::Int W, class T> struct hash<cplib::StaticMatrix<H, W, T>> {
    // 要素のハッシュを O(1) として行列のハッシュを O(HW+1) で求める。
    size_t operator()(const cplib::StaticMatrix<H, W, T> &m) const {
        return m.hash();
    }
};
}
