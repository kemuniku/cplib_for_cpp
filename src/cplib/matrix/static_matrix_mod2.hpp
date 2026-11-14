#pragma once
#include <cplib/matrix/bit_matrix_ops.hpp>
#include <bit>
#include <memory>
#include <string_view>

namespace cplib {
template <Int H, Int W> class StaticMatrixMod2 {
    static_assert(H >= 0 && W >= 0);
    template <Int, Int> friend class StaticMatrixMod2;
    static constexpr Int Stride = (W + 63) / 64;
    std::array<std::array<UInt, Stride>, H> rows_{};

public:
    using value_type = bool;
    // 全要素が0の行列を作る。計算量 O(H ceil(W/64))。
    StaticMatrixMod2() = default;

    // 整数の2次元配列をGF(2)の行列に変換する。計算量 O(HW)。
    template <std::integral T>
    explicit StaticMatrixMod2(const std::array<std::array<T, W>, H> &rows) {
        for (Int i = 0; i < H; ++i)
            for (Int j = 0; j < W; ++j)
                if (rows[i][j] & 1)
                    rows_[i][j >> 6] |= UInt(1) << (j & 63);
    }

    // 行数を返す。計算量 O(1)。
    constexpr Int h() const {
        return H;
    }

    // 列数を返す。計算量 O(1)。
    constexpr Int w() const {
        return W;
    }

    // 指定要素を読み出す。計算量 O(1)。
    bool operator()(Int i, Int j) const {
        assert(0 <= i && i < H && 0 <= j && j < W);
        return (rows_[i][j >> 6] >> (j & 63)) & 1;
    }

    // 指定要素への書き込み可能な参照を返す。計算量 O(1)。
    detail::MatrixBitReference operator()(Int i, Int j) {
        assert(0 <= i && i < H && 0 <= j && j < W);
        return {&rows_[i][j >> 6], UInt(1) << (j & 63)};
    }

    struct Row {
        StaticMatrixMod2 *matrix;
        Int row;

        // 行の指定要素への書き込み可能な参照を返す。計算量 O(1)。
        auto operator[](Int col) {
            return (*matrix)(row, col);
        }
    };

    struct ConstRow {
        const StaticMatrixMod2 *matrix;
        Int row;

        // 行の指定要素を返す。計算量 O(1)。
        bool operator[](Int col) const {
            return (*matrix)(row, col);
        }
    };

    // 書き込み可能な行ビューを返す。計算量 O(1)。
    Row operator[](Int row) {
        return {this, row};
    }

    // 読み取り専用の行ビューを返す。計算量 O(1)。
    ConstRow operator[](Int row) const {
        return {this, row};
    }

    // 行列が等しいか判定する。計算量 O(H ceil(W/64))。
    bool operator==(const StaticMatrixMod2 &) const = default;

    // 行列を文字列に変換する。計算量 O(H(W+1))。
    std::string str() const {
        std::string output_value;
        for (Int i = 0; i < H; ++i) {
            if (i)
                output_value += '\n';
            for (Int j = 0; j < W; ++j) {
                if (j)
                    output_value += ' ';
                output_value += (*this)(i, j) ? '1' : '0';
            }
        }
        return output_value;
    }

    // 文字列の0・1で行を設定する。計算量 O(ceil(W/64)+values.size())。
    void setRowBits(Int i, std::string_view values) {
        assert(0 <= i && i < H && Int(values.size()) <= W);
        rows_[i].fill(0);
        for (Int j = 0; j < Int(values.size()); ++j)
            if (values[j] == '1')
                rows_[i][j >> 6] |= UInt(1) << (j & 63);
    }

    // 行の先頭 width 要素を0・1の文字列で返す。計算量 O(width)。
    std::string rowBits(Int i, Int width = W) const {
        assert(0 <= i && i < H && 0 <= width && width <= W);
        std::string output_value(width, '0');
        for (Int j = 0; j < width; ++j)
            output_value[j] = char('0' + ((rows_[i][j >> 6] >> (j & 63)) & 1));
        return output_value;
    }

    // 単位行列を返す。計算量 O(H ceil(W/64)+H)。
    static StaticMatrixMod2 identity()
        requires(H == W)
    {
        StaticMatrixMod2 output_value;
        for (Int i = 0; i < H; ++i)
            output_value(i, i) = true;
        return output_value;
    }

    // 転置行列を返す。計算量 O(HW+W ceil(H/64))。
    auto transposed() const {
        StaticMatrixMod2<W, H> output_value;
        for (Int i = 0; i < H; ++i)
            for (Int j = 0; j < W; ++j)
                if ((*this)(i, j))
                    output_value(j, i) = true;
        return output_value;
    }

    // GF(2)の行列積を返す。右行列の列数 C に対し、最悪 O((H+1)(W+1)(C+1))。
    template <Int C> auto operator*(const StaticMatrixMod2<W, C> &b) const {
        StaticMatrixMod2<H, C> output_value;
        if constexpr (H < 40 || W < 64) {
            auto bt = b.transposed();
            for (Int i = 0; i < H; ++i)
                for (Int j = 0; j < C; ++j) {
                    int parity = 0;
                    for (Int k = 0; k < Stride; ++k)
                        parity ^= std::popcount(rows_[i][k] & bt.rows_[j][k]) & 1;
                    if (parity)
                        output_value(i, j) = true;
                }
        } else {
            constexpr Int count = (C + 63) / 64;
            std::array<std::array<UInt, count>, 256> table{};
            for (Int start = 0; start < W; start += 8) {
                Int bits = std::min<Int>(8, W - start);
                for (unsigned mask = 1; mask < (1u << bits); ++mask) {
                    Int previous = mask & (mask - 1), bit = std::countr_zero(mask);
                    for (Int k = 0; k < count; ++k)
                        table[mask][k] = table[previous][k] ^ b.rows_[start + bit][k];
                }
                UInt mask = (UInt(1) << bits) - 1;
                for (Int i = 0; i < H; ++i) {
                    Int index = (rows_[i][start >> 6] >> (start & 63)) & mask;
                    for (Int k = 0; k < count; ++k)
                        output_value.rows_[i][k] ^= table[index][k];
                }
            }
        }
        return output_value;
    }

    // 右から行列を掛ける。計算量は行列積と同じ。
    StaticMatrixMod2 &operator*=(const StaticMatrixMod2 &b)
        requires(H == W)
    {
        *this = (*this) * b;
        return *this;
    }

    // 正方行列の非負整数乗を返す。行列積を O(log(exponent+1)) 回行う。
    StaticMatrixMod2 pow(Int exponent) const
        requires(H == W)
    {
        assert(exponent >= 0);
        auto output_value = identity(), base = *this;
        for (; exponent; exponent >>= 1) {
            if (exponent & 1)
                output_value *= base;
            if (exponent > 1)
                base *= base;
        }
        return output_value;
    }

    // 元の固定長版と同じくrankの作業行列はヒープに置く。
    // 階数を求める。O(HW + H min(H,W) ceil(W/64))。空行列は O(1)。
    Int rank() const {
        if constexpr (H == 0 || W == 0)
            return 0;
        auto storage = std::make_unique<StaticMatrixMod2>(*this);
        auto &b = *storage;
        Int result = 0;
        for (Int col = 0; col < W; ++col) {
            Int pivot = result;
            while (pivot < H && !b(pivot, col))
                ++pivot;
            if (pivot == H)
                continue;
            std::swap(b.rows_[result], b.rows_[pivot]);
            for (Int i = result + 1; i < H; ++i)
                if (b(i, col))
                    for (Int k = 0; k < Stride; ++k)
                        b.rows_[i][k] ^= b.rows_[result][k];
            if (++result == H)
                break;
        }
        return result;
    }

    // 正方行列の行列式を返す。計算量は rank と同じ。
    bool determinant() const
        requires(H == W)
    {
        return rank() == H;
    }

    // 左右の固定長行列を別々に掃き出す。動的版の拡大配列とは別実装。
    // 逆行列を返す。正方行列の一辺を N として O(N² ceil(N/64))。
    std::optional<StaticMatrixMod2> inverse() const
        requires(H == W)
    {
        auto left = *this, right = identity();
        for (Int col = 0; col < H; ++col) {
            Int pivot = col;
            while (pivot < H && !left(pivot, col))
                ++pivot;
            if (pivot == H)
                return std::nullopt;
            std::swap(left.rows_[col], left.rows_[pivot]);
            std::swap(right.rows_[col], right.rows_[pivot]);
            for (Int i = 0; i < H; ++i)
                if (i != col && left(i, col)) {
                    for (Int k = 0; k < Stride; ++k)
                        left.rows_[i][k] ^= left.rows_[col][k];
                    for (Int k = 0; k < Stride; ++k)
                        right.rows_[i][k] ^= right.rows_[col][k];
                }
        }
        return right;
    }

    // 左上h行w列でAx=bをビット演算で解く。O(h*min(h,w)*(w div 64+1)+w^2)。
    // height/widthの省略時はH/W。元の行列は変更せず、解なしはstd::nulloptを返す。
    template <class B> auto solveLinearSystem(const B &b, Int height = H, Int width = W) const {
        assert(0 <= height && height <= H && 0 <= width && width <= W && Int(b.size()) == height);
        auto rows = initBitLinearSystem(height, width);
        Int stride = (width >> 6) + 1, full = width >> 6, tail = width & 63;
        for (Int i = 0; i < height; ++i) {
            for (Int k = 0; k < full; ++k)
                rows[i * stride + k] = rows_[i][k];
            if (tail)
                rows[i * stride + full] = rows_[i][full] & ((UInt(1) << tail) - 1);
            if (b[i])
                rows[i * stride + full] |= UInt(1) << tail;
        }
        return solveBitLinearSystem(rows, height, width);
    }

    // GF(2)上の対称な偶数次行列のhafnianを求める。O(n^3)。
    bool hafnian() const {
        assert(H == W);
        return fieldHafnian(matrixRows(*this, H, W));
    }

    // GF(2)上で特異行列も含めた余因子行列を求める。O(n^3)。
    StaticMatrixMod2 adjugate() const {
        assert(H == W);
        auto rows = *fieldAdjugateInverse(matrixRows(*this, H, W), true);
        StaticMatrixMod2 output_value;
        for (Int i = 0; i < H; ++i)
            for (Int j = 0; j < W; ++j)
                output_value(i, j) = bool(rows[i][j]);
        return output_value;
    }
};

// 零行列を作る。計算量は StaticMatrixMod2 の既定コンストラクタと同じ。
template <Int H, Int W> auto initStaticMatrixMod2() {
    return StaticMatrixMod2<H, W>{};
}

// 二次元配列をGF(2)の固定長行列に変換する。計算量 O(HW)。
template <std::size_t H, std::size_t W, std::integral T>
auto initStaticMatrixMod2(const std::array<std::array<T, W>, H> &rows) {
    return StaticMatrixMod2<H, W>(rows);
}

// 二次元配列をGF(2)の固定長行列に変換する。計算量 O(HW)。
template <std::size_t H, std::size_t W, std::integral T>
auto toStaticMatrixMod2(const std::array<std::array<T, W>, H> &rows) {
    return initStaticMatrixMod2(rows);
}

// 単位行列を作る。計算量は StaticMatrixMod2::identity() と同じ。
template <Int N> auto identityStaticMatrixMod2() {
    return StaticMatrixMod2<N, N>::identity();
}

// 行数を返す。計算量 O(1)。
template <Int H, Int W> Int h(const StaticMatrixMod2<H, W> &a) {
    return a.h();
}

// 列数を返す。計算量 O(1)。
template <Int H, Int W> Int w(const StaticMatrixMod2<H, W> &a) {
    return a.w();
}

// 行列の文字列を返す。計算量は a.str と同じ。
template <Int H, Int W> auto to_string(const StaticMatrixMod2<H, W> &a) {
    return a.str();
}

// 文字列で行を設定する。計算量は a.setRowBits と同じ。
template <Int H, Int W> void setRowBits(StaticMatrixMod2<H, W> &a, Int i, std::string_view values) {
    a.setRowBits(i, values);
}

// 行列の行のビット文字列を返す。計算量は a.rowBits と同じ。
template <Int H, Int W> auto rowBits(const StaticMatrixMod2<H, W> &a, Int i, Int width = W) {
    return a.rowBits(i, width);
}

// 行列の転置行列を返す。計算量は a.transposed と同じ。
template <Int H, Int W> auto transposed(const StaticMatrixMod2<H, W> &a) {
    return a.transposed();
}

// 行列の非負整数乗を返す。計算量は a.pow と同じ。
template <Int N> auto pow(const StaticMatrixMod2<N, N> &a, Int exponent) {
    return a.pow(exponent);
}

// 行列の階数を返す。計算量は a.rank と同じ。
template <Int H, Int W> Int rank(const StaticMatrixMod2<H, W> &a) {
    return a.rank();
}

// 行列の行列式を返す。計算量は a.determinant と同じ。
template <Int N> bool determinant(const StaticMatrixMod2<N, N> &a) {
    return a.determinant();
}

// 行列の逆行列を返す。計算量は a.inverse と同じ。
template <Int N> auto inverse(const StaticMatrixMod2<N, N> &a) {
    return a.inverse();
}

// 行列のhafnianを返す。計算量は a.hafnian と同じ。
template <Int H, Int W> bool hafnian(const StaticMatrixMod2<H, W> &a) {
    return a.hafnian();
}

// 行列の余因子行列を返す。計算量は a.adjugate と同じ。
template <Int H, Int W> auto adjugate(const StaticMatrixMod2<H, W> &a) {
    return a.adjugate();
}

// 左上の係数行列で連立方程式を解く。計算量は a.solveLinearSystem と同じ。
template <Int H, Int W, class B>
auto solveLinearSystem(const StaticMatrixMod2<H, W> &a, const B &b, Int height = H, Int width = W) {
    return a.solveLinearSystem(b, height, width);
}
}
