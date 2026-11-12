#pragma once
#include <cplib/matrix/bit_matrix_ops.hpp>
#include <bit>
#include <string_view>

namespace cplib {
// GF(2)上のH行W列の行列。次元と要素・行ビュー取得は O(1)、要素代入も O(1)。
// 64ビット語で保持し、等値比較は O(H ceil(W/64)+1)。
class MatrixMod2 {
    Int height_ = 0, width_ = 0, stride_ = 0;
    std::vector<UInt> words_;

    UInt &word(Int i, Int k) {
        return words_[i * stride_ + k];
    }

    UInt word(Int i, Int k) const {
        return words_[i * stride_ + k];
    }

    void swap_rows(Int i, Int j) {
        if (i != j)
            for (Int k = 0; k < stride_; ++k)
                std::swap(word(i, k), word(j, k));
    }

public:
    using value_type = bool;
    // 空行列を作る。計算量 O(1)。
    MatrixMod2() = default;

    // H行W列のGF(2)零行列を作る。行ごとの確保は行わない。計算量 O(H ceil(W/64)+1)。
    MatrixMod2(Int height, Int width)
        : height_(height), width_(width), stride_((width >> 6) + bool(width & 63)) {
        assert(height >= 0 && width >= 0);
        assert(!stride_ || height <= std::numeric_limits<Int>::max() / Int(sizeof(UInt)) / stride_);
        words_.resize(height * stride_);
    }

    // 整数配列の各要素を2で割った余りからH行W列の行列を作る。O(HW+H)。
    template <std::integral T>
    explicit MatrixMod2(const std::vector<std::vector<T>> &rows)
        : MatrixMod2(rows.size(), rows.empty() ? 0 : rows[0].size()) {
        for (Int i = 0; i < height_; ++i) {
            assert(Int(rows[i].size()) == width_);
            for (Int j = 0; j < width_; ++j)
                if (rows[i][j] & 1)
                    word(i, j >> 6) |= UInt(1) << (j & 63);
        }
    }

    Int h() const {
        return height_;
    }

    Int w() const {
        return width_;
    }

    bool operator()(Int i, Int j) const {
        assert(0 <= i && i < h() && 0 <= j && j < w());
        return (word(i, j >> 6) >> (j & 63)) & 1;
    }

    detail::MatrixBitReference operator()(Int i, Int j) {
        assert(0 <= i && i < h() && 0 <= j && j < w());
        return {&word(i, j >> 6), UInt(1) << (j & 63)};
    }

    struct Row {
        MatrixMod2 *matrix;
        Int row;

        auto operator[](Int col) {
            return (*matrix)(row, col);
        }
    };

    struct ConstRow {
        const MatrixMod2 *matrix;
        Int row;

        bool operator[](Int col) const {
            return (*matrix)(row, col);
        }
    };

    Row operator[](Int row) {
        return {this, row};
    }

    ConstRow operator[](Int row) const {
        return {this, row};
    }

    bool operator==(const MatrixMod2 &) const = default;

    // 行列を文字列に変換する。O(HW+H)。
    std::string str() const {
        std::string out;
        for (Int i = 0; i < h(); ++i) {
            if (i)
                out += '\n';
            for (Int j = 0; j < w(); ++j) {
                if (j)
                    out += ' ';
                out += (*this)(i, j) ? '1' : '0';
            }
        }
        return out;
    }

    // 指定行に文字列のビットを設定する。列数Wとして O(W+1)。
    void setRowBits(Int i, std::string_view values) {
        assert(0 <= i && i < h() && Int(values.size()) <= w());
        for (Int k = 0; k < stride_; ++k)
            word(i, k) = 0;
        for (Int j = 0; j < Int(values.size()); ++j)
            if (values[j] == '1')
                word(i, j >> 6) |= UInt(1) << (j & 63);
    }

    // 指定行の先頭width列を文字列で返す。O(width+1)。
    std::string rowBits(Int i, Int width) const {
        assert(0 <= i && i < h() && 0 <= width && width <= w());
        std::string out(width, '0');
        for (Int j = 0; j < width; ++j)
            out[j] = char('0' + ((word(i, j >> 6) >> (j & 63)) & 1));
        return out;
    }

    // 指定行を文字列で返す。列数Wとして O(W+1)。
    std::string rowBits(Int i) const {
        return rowBits(i, w());
    }

    // n次の単位行列を作る。O(n ceil(n/64)+n+1)。
    static MatrixMod2 identity(Int n) {
        MatrixMod2 out(n, n);
        for (Int i = 0; i < n; ++i)
            out(i, i) = true;
        return out;
    }

    // 転置行列を返す。O(HW+H+W)。
    MatrixMod2 transposed() const {
        MatrixMod2 out(w(), h());
        for (Int i = 0; i < h(); ++i)
            for (Int j = 0; j < w(); ++j)
                if ((*this)(i, j))
                    out(j, i) = true;
        return out;
    }

    // 小さい場合は転置とpopcount、大きい場合は8列ずつのFour Russians法。
    // H×WとW×Cの行列積を返す。O((H+1)(W+1)(C+1))。
    MatrixMod2 operator*(const MatrixMod2 &b) const {
        assert(w() == b.h());
        MatrixMod2 out(h(), b.w());
        if (h() < 40 || w() < 64) {
            auto bt = b.transposed();
            for (Int i = 0; i < h(); ++i)
                for (Int j = 0; j < b.w(); ++j) {
                    int parity = 0;
                    for (Int k = 0; k < stride_; ++k)
                        parity ^= std::popcount(word(i, k) & bt.word(j, k)) & 1;
                    if (parity)
                        out(i, j) = true;
                }
        } else {
            Int count = b.stride_;
            if (!count)
                return out;
            std::vector<UInt> table(256 * count);
            for (Int start = 0; start < w(); start += 8) {
                Int bits = std::min<Int>(8, w() - start);
                for (unsigned mask = 1; mask < (1u << bits); ++mask) {
                    Int previous = mask & (mask - 1), bit = std::countr_zero(mask);
                    for (Int k = 0; k < count; ++k)
                        table[mask * count + k] =
                            table[previous * count + k] ^ b.word(start + bit, k);
                }
                UInt mask = (UInt(1) << bits) - 1;
                for (Int i = 0; i < h(); ++i) {
                    Int index = (word(i, start >> 6) >> (start & 63)) & mask;
                    for (Int k = 0; k < count; ++k)
                        out.word(i, k) ^= table[index * count + k];
                }
            }
        }
        return out;
    }

    MatrixMod2 &operator*=(const MatrixMod2 &b) {
        *this = (*this) * b;
        return *this;
    }

    // n次正方行列を非負のk乗にする。O((n+1)³ log(k+2))。
    MatrixMod2 pow(Int exponent) const {
        assert(h() == w() && exponent >= 0);
        auto out = identity(h()), base = *this;
        for (; exponent; exponent >>= 1) {
            if (exponent & 1)
                out *= base;
            if (exponent > 1)
                base *= base;
        }
        return out;
    }

    // 階数をO(h*w+h*min(h,w)*ceil(w/64))で求める。元の行列は変更しない。空行列はO(1)。
    // 階数を返す。O((H+1)(W+1)(ceil(W/64)+1))。
    Int rank() const {
        if (!h() || !w())
            return 0;
        auto b = *this;
        Int result = 0;
        for (Int col = 0; col < w(); ++col) {
            Int pivot = result;
            while (pivot < h() && !b(pivot, col))
                ++pivot;
            if (pivot == h())
                continue;
            b.swap_rows(result, pivot);
            for (Int i = result + 1; i < h(); ++i)
                if (b(i, col))
                    for (Int k = 0; k < stride_; ++k)
                        b.word(i, k) ^= b.word(result, k);
            if (++result == h())
                break;
        }
        return result;
    }

    // n次正方行列の行列式を返す。O((n+1)²(ceil(n/64)+1))。
    bool determinant() const {
        assert(h() == w());
        return rank() == h();
    }

    // 右側の単位行列を語境界に置く動的版の掃き出し法。
    // 64bit単位の掃き出し法で逆行列を求める。O(n^2*ceil(n/64))。特異行列はstd::nullopt。
    // 元の行列は変更しない。作業領域はO(n*ceil(n/64))。
    // n次正方行列の逆行列を返す。特異行列は空。O((n+1)²(ceil(n/64)+1))。
    std::optional<MatrixMod2> inverse() const {
        assert(h() == w());
        Int n = h();
        if (!n)
            return MatrixMod2(0, 0);
        Int stride = 2 * stride_;
        assert(n <= std::numeric_limits<Int>::max() / Int(sizeof(UInt)) / stride);
        std::vector<UInt> data(n * stride);
        for (Int i = 0; i < n; ++i) {
            std::copy_n(words_.begin() + i * stride_, stride_, data.begin() + i * stride);
            data[i * stride + stride_ + (i >> 6)] = UInt(1) << (i & 63);
        }
        for (Int col = 0; col < n; ++col) {
            Int first = col >> 6, pivot = col;
            UInt mask = UInt(1) << (col & 63);
            while (pivot < n && !(data[pivot * stride + first] & mask))
                ++pivot;
            if (pivot == n)
                return std::nullopt;
            UInt *p = data.data() + col * stride;
            if (pivot != col)
                for (Int k = first; k < stride; ++k)
                    std::swap(p[k], data[pivot * stride + k]);
            for (Int i = 0; i < n; ++i) {
                if (i == col)
                    continue;
                UInt *row = data.data() + i * stride;
                if (row[first] & mask)
                    for (Int k = first; k < stride; ++k)
                        row[k] ^= p[k];
            }
        }
        MatrixMod2 out(n, n);
        for (Int i = 0; i < n; ++i)
            std::copy_n(data.begin() + i * stride + stride_, stride_,
                        out.words_.begin() + i * stride_);
        return out;
    }

    // ビット演算でAx=bの特殊解と核の基底を返す。O(h*min(h,w)*(w/64+1)+w^2)。
    // 元の行列は変更しない。解なしはstd::nullopt、基底の個数はw-rank。
    template <class B> auto solveLinearSystem(const B &b) const {
        assert(Int(b.size()) == h());
        auto rows = initBitLinearSystem(h(), w());
        Int stride = (w() >> 6) + 1;
        for (Int i = 0; i < h(); ++i) {
            for (Int k = 0; k < stride_; ++k)
                rows[i * stride + k] = word(i, k);
            if (b[i])
                rows[i * stride + (w() >> 6)] |= UInt(1) << (w() & 63);
        }
        return solveBitLinearSystem(rows, h(), w());
    }

    // GF(2)上の対称な偶数次行列のhafnianを求める。O(n^3)。
    bool hafnian() const {
        assert(h() == w());
        return fieldHafnian(matrixRows(*this, h(), w()));
    }

    // GF(2)上で特異行列も含めた余因子行列を求める。O(n^3)。
    MatrixMod2 adjugate() const {
        assert(h() == w());
        auto rows = *fieldAdjugateInverse(matrixRows(*this, h(), w()), true);
        return MatrixMod2(rows);
    }
};

// H行W列の零行列を作る。O(H ceil(W/64)+1)。
inline auto initMatrixMod2(Int h, Int w) {
    return MatrixMod2(h, w);
}

// 整数配列の偶奇からH行W列の行列を作る。O(HW+H)。
template <std::integral T> auto initMatrixMod2(const std::vector<std::vector<T>> &rows) {
    return MatrixMod2(rows);
}

// 整数配列の偶奇からH行W列の行列を作る。O(HW+H)。
template <std::integral T> auto toMatrixMod2(const std::vector<std::vector<T>> &rows) {
    return MatrixMod2(rows);
}

// 行数を返す。O(1)。
inline Int h(const MatrixMod2 &a) {
    return a.h();
}

// 列数を返す。O(1)。
inline Int w(const MatrixMod2 &a) {
    return a.w();
}

// H行W列の行列を文字列に変換する。O(HW+H)。
inline auto to_string(const MatrixMod2 &a) {
    return a.str();
}

// 指定行のビットを設定する。列数Wとして O(W+1)。
inline void setRowBits(MatrixMod2 &a, Int i, std::string_view values) {
    a.setRowBits(i, values);
}

// 指定行を文字列で返す。列数Wとして O(W+1)。
inline auto rowBits(const MatrixMod2 &a, Int i) {
    return a.rowBits(i);
}

// 指定行の先頭w列を文字列で返す。O(w+1)。
inline auto rowBits(const MatrixMod2 &a, Int i, Int w) {
    return a.rowBits(i, w);
}

// n次の単位行列を返す。O(n ceil(n/64)+n+1)。
inline auto identityMatrixMod2(Int n) {
    return MatrixMod2::identity(n);
}

// H行W列の転置行列を返す。O(HW+H+W)。
inline auto transposed(const MatrixMod2 &a) {
    return a.transposed();
}

// n次正方行列を非負のk乗にする。O((n+1)³ log(k+2))。
inline auto pow(const MatrixMod2 &a, Int n) {
    return a.pow(n);
}

// H行W列の階数を返す。O((H+1)(W+1)(ceil(W/64)+1))。
inline Int rank(const MatrixMod2 &a) {
    return a.rank();
}

// n次正方行列の行列式を返す。O((n+1)²(ceil(n/64)+1))。
inline bool determinant(const MatrixMod2 &a) {
    return a.determinant();
}

// n次正方行列の逆行列を返す。O((n+1)²(ceil(n/64)+1))。
inline auto inverse(const MatrixMod2 &a) {
    return a.inverse();
}

// n次対称偶数次行列のhafnianを返す。O(n³+1)。
inline bool hafnian(const MatrixMod2 &a) {
    return a.hafnian();
}

// n次正方行列の余因子行列を返す。O(n³+1)。
inline auto adjugate(const MatrixMod2 &a) {
    return a.adjugate();
}

// Ax=bの特殊解と核の基底を返す。行数H、列数Wとして O(H min(H,W)(W/64+1)+W²+H+W)。
template <class B> auto solveLinearSystem(const MatrixMod2 &a, const B &b) {
    return a.solveLinearSystem(b);
}
}
