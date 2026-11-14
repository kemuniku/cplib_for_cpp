#pragma once
#include <cplib/matrix/field_matrix_ops.hpp>
#include <sstream>

namespace cplib {
// 動的版は元と同じ行ごとの配列。0行の場合も列数を保持する。
// 行数H、列数W、要素演算をO(1)として、要素ごとの算術・ビット演算と比較は O(HW+H)。
// 行列積 H×W と W×C は O(HWC+HC+H)。行・要素へのアクセスと行列の寸法取得は O(1)。
template <class T> class Matrix {
    std::vector<std::vector<T>> arr_;
    Int empty_width_ = 0;

public:
    using value_type = T;
    // 空行列を作る。計算量 O(1)。
    Matrix() = default;

    // 指定した値でH行W列を初期化する。計算量 O(HW+H)。
    Matrix(Int height, Int width, T value = T(0))
        : arr_(height, std::vector<T>(width, value)), empty_width_(width) {
        assert(height >= 0 && width >= 0);
    }

    // 行列をコピーして構築する。計算量 O(HW+H)。
    explicit Matrix(const std::vector<std::vector<T>> &rows) : arr_(rows) {
        for (const auto &row : rows)
            assert(row.size() == (rows.empty() ? 0 : rows[0].size()));
    }

    // N要素から行または列ベクトルを作る。計算量 O(N)。
    explicit Matrix(const std::vector<T> &values, bool vertical = false) {
        if (vertical) {
            empty_width_ = 1;
            for (T value : values)
                arr_.push_back({value});
        } else
            arr_.push_back(values);
    }

    Int h() const {
        return arr_.size();
    }

    Int w() const {
        return arr_.empty() ? empty_width_ : Int(arr_[0].size());
    }

    auto &operator[](Int row) {
        return arr_[row];
    }

    const auto &operator[](Int row) const {
        return arr_[row];
    }

    decltype(auto) operator()(Int row, Int col) {
        return arr_[row][col];
    }

    decltype(auto) operator()(Int row, Int col) const {
        return arr_[row][col];
    }

    bool operator==(const Matrix &b) const {
        return h() == b.h() && w() == b.w() && arr_ == b.arr_;
    }

    Matrix operator-() const {
        auto output_value = *this;
        for (Int i = 0; i < h(); ++i)
            for (Int j = 0; j < w(); ++j)
                output_value(i, j) = -(*this)(i, j);
        return output_value;
    }

    Matrix &operator*=(const Matrix &b) {
        assert(w() == b.h());
        Matrix output_value(h(), b.w(), T(0));
        for (Int i = 0; i < h(); ++i)
            for (Int j = 0; j < b.w(); ++j)
                for (Int k = 0; k < w(); ++k)
                    output_value(i, j) += (*this)(i, k) * b(k, j);
        *this = std::move(output_value);
        return *this;
    }

    Matrix &operator*=(T value) {
        for (Int i = 0; i < h(); ++i)
            for (Int j = 0; j < w(); ++j)
                (*this)(i, j) *= value;
        return *this;
    }

    friend Matrix operator*(Matrix a, const Matrix &b) {
        a *= b;
        return a;
    }

    friend Matrix operator*(Matrix a, T b) {
        a *= b;
        return a;
    }

    friend Matrix operator*(T a, Matrix b) {
        b *= a;
        return b;
    }

    Matrix &operator+=(const Matrix &b) {
        assert(h() == b.h() && w() == b.w());
        for (Int i = 0; i < h(); ++i)
            for (Int j = 0; j < w(); ++j) {
                T a = (*this)(i, j);
                T rhs = b(i, j);
                (*this)(i, j) = a + rhs;
            }
        return *this;
    }

    Matrix &operator+=(T b) {
        for (Int i = 0; i < h(); ++i)
            for (Int j = 0; j < w(); ++j) {
                T a = (*this)(i, j);
                (*this)(i, j) = a + b;
            }
        return *this;
    }

    friend Matrix operator+(Matrix a, const Matrix &b) {
        a += b;
        return a;
    }

    friend Matrix operator+(Matrix a, T b) {
        a += b;
        return a;
    }

    friend Matrix operator+(T a, Matrix matrix) {
        for (Int i = 0; i < matrix.h(); ++i)
            for (Int j = 0; j < matrix.w(); ++j) {
                T b = matrix(i, j);
                matrix(i, j) = a + b;
            }
        return matrix;
    }

    Matrix &operator-=(const Matrix &b) {
        assert(h() == b.h() && w() == b.w());
        for (Int i = 0; i < h(); ++i)
            for (Int j = 0; j < w(); ++j) {
                T a = (*this)(i, j);
                T rhs = b(i, j);
                (*this)(i, j) = a - rhs;
            }
        return *this;
    }

    Matrix &operator-=(T b) {
        for (Int i = 0; i < h(); ++i)
            for (Int j = 0; j < w(); ++j) {
                T a = (*this)(i, j);
                (*this)(i, j) = a - b;
            }
        return *this;
    }

    friend Matrix operator-(Matrix a, const Matrix &b) {
        a -= b;
        return a;
    }

    friend Matrix operator-(Matrix a, T b) {
        a -= b;
        return a;
    }

    friend Matrix operator-(T a, Matrix matrix) {
        for (Int i = 0; i < matrix.h(); ++i)
            for (Int j = 0; j < matrix.w(); ++j) {
                T b = matrix(i, j);
                matrix(i, j) = a - b;
            }
        return matrix;
    }

    Matrix &operator&=(const Matrix &b)
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

    Matrix &operator&=(T b)
        requires(std::integral<T>)
    {
        for (Int i = 0; i < h(); ++i)
            for (Int j = 0; j < w(); ++j) {
                T a = (*this)(i, j);
                (*this)(i, j) = a & b;
            }
        return *this;
    }

    friend Matrix operator&(Matrix a, const Matrix &b)
        requires(std::integral<T>)
    {
        a &= b;
        return a;
    }

    friend Matrix operator&(Matrix a, T b)
        requires(std::integral<T>)
    {
        a &= b;
        return a;
    }

    friend Matrix operator&(T a, Matrix matrix)
        requires(std::integral<T>)
    {
        for (Int i = 0; i < matrix.h(); ++i)
            for (Int j = 0; j < matrix.w(); ++j) {
                T b = matrix(i, j);
                matrix(i, j) = a & b;
            }
        return matrix;
    }

    Matrix &operator|=(const Matrix &b)
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

    Matrix &operator|=(T b)
        requires(std::integral<T>)
    {
        for (Int i = 0; i < h(); ++i)
            for (Int j = 0; j < w(); ++j) {
                T a = (*this)(i, j);
                (*this)(i, j) = a | b;
            }
        return *this;
    }

    friend Matrix operator|(Matrix a, const Matrix &b)
        requires(std::integral<T>)
    {
        a |= b;
        return a;
    }

    friend Matrix operator|(Matrix a, T b)
        requires(std::integral<T>)
    {
        a |= b;
        return a;
    }

    friend Matrix operator|(T a, Matrix matrix)
        requires(std::integral<T>)
    {
        for (Int i = 0; i < matrix.h(); ++i)
            for (Int j = 0; j < matrix.w(); ++j) {
                T b = matrix(i, j);
                matrix(i, j) = a | b;
            }
        return matrix;
    }

    Matrix &operator^=(const Matrix &b)
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

    Matrix &operator^=(T b)
        requires(std::integral<T>)
    {
        for (Int i = 0; i < h(); ++i)
            for (Int j = 0; j < w(); ++j) {
                T a = (*this)(i, j);
                (*this)(i, j) = a ^ b;
            }
        return *this;
    }

    friend Matrix operator^(Matrix a, const Matrix &b)
        requires(std::integral<T>)
    {
        a ^= b;
        return a;
    }

    friend Matrix operator^(Matrix a, T b)
        requires(std::integral<T>)
    {
        a ^= b;
        return a;
    }

    friend Matrix operator^(T a, Matrix matrix)
        requires(std::integral<T>)
    {
        for (Int i = 0; i < matrix.h(); ++i)
            for (Int j = 0; j < matrix.w(); ++j) {
                T b = matrix(i, j);
                matrix(i, j) = a ^ b;
            }
        return matrix;
    }

    Matrix &operator<<=(const Matrix &b)
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

    Matrix &operator<<=(T b)
        requires(std::integral<T>)
    {
        for (Int i = 0; i < h(); ++i)
            for (Int j = 0; j < w(); ++j) {
                T a = (*this)(i, j);
                (*this)(i, j) = T(std::make_unsigned_t<T>(a) << b);
            }
        return *this;
    }

    friend Matrix operator<<(Matrix a, const Matrix &b)
        requires(std::integral<T>)
    {
        a <<= b;
        return a;
    }

    friend Matrix operator<<(Matrix a, T b)
        requires(std::integral<T>)
    {
        a <<= b;
        return a;
    }

    friend Matrix operator<<(T a, Matrix matrix)
        requires(std::integral<T>)
    {
        for (Int i = 0; i < matrix.h(); ++i)
            for (Int j = 0; j < matrix.w(); ++j) {
                T b = matrix(i, j);
                matrix(i, j) = T(std::make_unsigned_t<T>(a) << b);
            }
        return matrix;
    }

    Matrix &operator>>=(const Matrix &b)
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

    Matrix &operator>>=(T b)
        requires(std::integral<T>)
    {
        for (Int i = 0; i < h(); ++i)
            for (Int j = 0; j < w(); ++j) {
                T a = (*this)(i, j);
                (*this)(i, j) = a >> b;
            }
        return *this;
    }

    friend Matrix operator>>(Matrix a, const Matrix &b)
        requires(std::integral<T>)
    {
        a >>= b;
        return a;
    }

    friend Matrix operator>>(Matrix a, T b)
        requires(std::integral<T>)
    {
        a >>= b;
        return a;
    }

    friend Matrix operator>>(T a, Matrix matrix)
        requires(std::integral<T>)
    {
        for (Int i = 0; i < matrix.h(); ++i)
            for (Int j = 0; j < matrix.w(); ++j) {
                T b = matrix(i, j);
                matrix(i, j) = a >> b;
            }
        return matrix;
    }

    Matrix &operator/=(const Matrix &b)
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

    Matrix &operator/=(T b)
        requires(std::integral<T>)
    {
        for (Int i = 0; i < h(); ++i)
            for (Int j = 0; j < w(); ++j) {
                T a = (*this)(i, j);
                (*this)(i, j) = a / b;
            }
        return *this;
    }

    friend Matrix operator/(Matrix a, const Matrix &b)
        requires(std::integral<T>)
    {
        a /= b;
        return a;
    }

    friend Matrix operator/(Matrix a, T b)
        requires(std::integral<T>)
    {
        a /= b;
        return a;
    }

    friend Matrix operator/(T a, Matrix matrix)
        requires(std::integral<T>)
    {
        for (Int i = 0; i < matrix.h(); ++i)
            for (Int j = 0; j < matrix.w(); ++j) {
                T b = matrix(i, j);
                matrix(i, j) = a / b;
            }
        return matrix;
    }

    Matrix &operator%=(const Matrix &b)
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

    Matrix &operator%=(T b)
        requires(std::integral<T>)
    {
        for (Int i = 0; i < h(); ++i)
            for (Int j = 0; j < w(); ++j) {
                T a = (*this)(i, j);
                (*this)(i, j) = a % b;
            }
        return *this;
    }

    friend Matrix operator%(Matrix a, const Matrix &b)
        requires(std::integral<T>)
    {
        a %= b;
        return a;
    }

    friend Matrix operator%(Matrix a, T b)
        requires(std::integral<T>)
    {
        a %= b;
        return a;
    }

    friend Matrix operator%(T a, Matrix matrix)
        requires(std::integral<T>)
    {
        for (Int i = 0; i < matrix.h(); ++i)
            for (Int j = 0; j < matrix.w(); ++j) {
                T b = matrix(i, j);
                matrix(i, j) = a % b;
            }
        return matrix;
    }

    // n次の単位行列を作る。計算量 O(n²+n)。
    static Matrix identity(Int n, T one = T(1), T zero = T(0)) {
        Matrix output_value(n, n, zero);
        for (Int i = 0; i < n; ++i)
            output_value(i, i) = one;
        return output_value;
    }

    // 正方行列を非負の指数kで累乗する。次元nとして O(n³ log(k+1)+n²+n)。
    Matrix pow(Int exponent) const {
        Matrix result = identity(h()), base = *this;
        while (exponent > 0) {
            if (exponent & 1)
                result *= base;
            base *= base;
            exponent >>= 1;
        }
        return result;
    }

    // 全要素の和を返す。計算量 O(HW+H)。
    T sum() const {
        T result = T(0);
        for (const auto &row : arr_)
            for (T value : row)
                result += value;
        return result;
    }

    // 行列を文字列に変換する。要素の出力時間の合計と O(HW+H) が必要。
    std::string str() const {
        std::ostringstream output_value;
        for (Int i = 0; i < h(); ++i) {
            if (i)
                output_value << '\n';
            for (Int j = 0; j < w(); ++j) {
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

    // 行列のハッシュを返す。要素のハッシュを O(1) として O(HW+H)。
    std::size_t hash() const {
        std::size_t seed = arr_.size();
        for (const auto &row : arr_) {
            seed ^= row.size() + 0x9e3779b9 + (seed << 6) + (seed >> 2);
            for (T value : row)
                seed ^= std::hash<T>{}(value) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
        }
        return seed;
    }

    // 階数を求める。O(h*w*min(h,w))。
    Int rank() const {
        return fieldRank(matrixRows(*this, h(), w()), w());
    }

    // 行列式を求める。空行列は1。O(n^3)。
    T determinant() const {
        assert(h() == w());
        return fieldDeterminant(matrixRows(*this, h(), w()));
    }

    // 対称な偶数次行列のhafnianを求める。O(n^2*2^(n/2))。
    T hafnian() const {
        assert(h() == w());
        return fieldHafnian(matrixRows(*this, h(), w()));
    }

    // Ax=bの特殊解と核の基底を返す。解なしはstd::nullopt。計算量 O((H+1)(W+1)²)。
    template <class B> auto solveLinearSystem(const B &b) const {
        return fieldSolve(matrixRows(*this, h(), w()), w(), b);
    }

    // 逆行列を返す。特異行列はstd::nullopt。O(n^3)。
    std::optional<Matrix> inverse() const {
        assert(h() == w());
        auto rows = fieldAdjugateInverse(matrixRows(*this, h(), w()), false);
        if (!rows)
            return std::nullopt;
        return Matrix(*rows);
    }

    // 特異行列を含む余因子行列を返す。O(n^3)。
    Matrix adjugate() const {
        assert(h() == w());
        return Matrix(*fieldAdjugateInverse(matrixRows(*this, h(), w()), true));
    }
};

// 行列のコピーから構築する。計算量 O(HW+H)。
template <class T> auto initMatrix(const std::vector<std::vector<T>> &rows) {
    return Matrix<T>(rows);
}

// 行列のコピーから構築する。計算量 O(HW+H)。
template <class T> auto toMatrix(const std::vector<std::vector<T>> &rows) {
    return Matrix<T>(rows);
}

// N要素の行または列ベクトルを構築する。計算量 O(N)。
template <class T> auto initMatrix(const std::vector<T> &values, bool vertical = false) {
    return Matrix<T>(values, vertical);
}

// 指定値でH行W列を構築する。計算量 O(HW+H)。
template <class T> auto initMatrix(Int h, Int w, T value) {
    return Matrix<T>(h, w, value);
}

// n次単位行列を返す。計算量 O(n²+n)。
template <class T> auto identity_matrix(Int n, T one = T(1), T zero = T(0)) {
    return Matrix<T>::identity(n, one, zero);
}

// 行数を返す。計算量 O(1)。
template <class T> Int h(const Matrix<T> &m) {
    return m.h();
}

// 列数を返す。計算量 O(1)。
template <class T> Int w(const Matrix<T> &m) {
    return m.w();
}

// 行列を文字列に変換する。要素の出力時間の合計と O(HW+H) が必要。
template <class T> auto to_string(const Matrix<T> &m) {
    return m.str();
}

// 要素のハッシュを O(1) として行列のハッシュを O(HW+H) で求める。
template <class T> auto hashValue(const Matrix<T> &m) {
    return m.hash();
}

// n次正方行列を非負のk乗にする。計算量 O(n³ log(k+1)+n²+n)。
template <class T> auto pow(const Matrix<T> &m, Int n) {
    return m.pow(n);
}

// 全要素の和を返す。計算量 O(HW+H)。
template <class T> T sum(const Matrix<T> &m) {
    return m.sum();
}

// 行列の階数を返す。計算量 O(HW min(H,W)+H+W)。
template <class T> Int rank(const Matrix<T> &m) {
    return m.rank();
}

// n次正方行列の行列式を返す。計算量 O(n³+n)。
template <class T> T determinant(const Matrix<T> &m) {
    return m.determinant();
}

// n次の対称な偶数次行列のhafnianを返す。計算量 O(n² 2^(n/2)+1)。
template <class T> T hafnian(const Matrix<T> &m) {
    return m.hafnian();
}

// Ax=bの特殊解と核の基底を返す。計算量 O((H+1)(W+1)²)。
template <class T, class B> auto solveLinearSystem(const Matrix<T> &m, const B &b) {
    return m.solveLinearSystem(b);
}

// n次正方行列の逆行列を返す。存在しない場合は空。計算量 O(n³+n)。
template <class T> auto inverse(const Matrix<T> &m) {
    return m.inverse();
}

// n次正方行列の余因子行列を返す。計算量 O(n³+n)。
template <class T> auto adjugate(const Matrix<T> &m) {
    return m.adjugate();
}
}

namespace std {
template <class T> struct hash<cplib::Matrix<T>> {
    // 要素のハッシュを O(1) として行列のハッシュを O(HW+H) で求める。
    size_t operator()(const cplib::Matrix<T> &m) const {
        return m.hash();
    }
};
}
