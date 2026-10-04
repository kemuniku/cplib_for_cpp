#pragma once
#include <cplib/matrix/field_matrix_ops.hpp>

namespace cplib {
// 右辺を含むh行(w+1)列の作業領域を確保する。O(h*(w div 64+1))。
inline std::vector<UInt> initBitLinearSystem(Int height, Int width) {
    assert(height >= 0 && width >= 0);
    Int stride = (width >> 6) + 1;
    assert(height <= std::numeric_limits<Int>::max() / Int(sizeof(UInt)) / stride);
    return std::vector<UInt>(height * stride);
}

// 拡大行列を64bit単位で掃き出す。O(h min(h,w)(w/64+1)+w²)。
// 各行はL = width / 64 + 1ワードで、右辺はwidth列目。rowsを変更する。
// 解なしならstd::nulloptを返し、解があれば特殊解と核の基底を返す。
inline std::optional<LinearSystemSolution<bool>> solveBitLinearSystem(std::vector<UInt> &rows,
                                                                      Int height, Int width) {
    assert(height >= 0 && width >= 0);
    Int stride = (width >> 6) + 1;
    assert(height <= std::numeric_limits<Int>::max() / stride &&
           Int(rows.size()) == height * stride);
    std::vector<Int> pivots;
    for (Int col = 0; col < width; ++col) {
        Int rank = pivots.size();
        if (rank == height)
            break;
        Int first = col >> 6;
        UInt mask = UInt(1) << (col & 63);
        Int pivot = rank;
        while (pivot < height && !(rows[pivot * stride + first] & mask))
            ++pivot;
        if (pivot == height)
            continue;
        UInt *p = rows.data() + rank * stride;
        if (pivot != rank)
            for (Int k = first; k < stride; ++k)
                std::swap(p[k], rows[pivot * stride + k]);
        for (Int i = 0; i < height; ++i) {
            if (i == rank)
                continue;
            UInt *row = rows.data() + i * stride;
            if (row[first] & mask)
                for (Int k = first; k < stride; ++k)
                    row[k] ^= p[k];
        }
        pivots.push_back(col);
    }
    for (Int i = pivots.size(); i < height; ++i)
        if (rows[i * stride + (width >> 6)] & (UInt(1) << (width & 63)))
            return std::nullopt;
    LinearSystemSolution<bool> solution;
    solution.particular.resize(width);
    std::vector<bool> is_pivot(width);
    for (Int i = 0; i < Int(pivots.size()); ++i) {
        Int col = pivots[i];
        is_pivot[col] = true;
        solution.particular[col] = (rows[i * stride + (width >> 6)] >> (width & 63)) & 1;
    }
    for (Int free = 0; free < width; ++free) {
        if (is_pivot[free])
            continue;
        std::vector<bool> vector(width);
        vector[free] = true;
        for (Int i = 0; i < Int(pivots.size()); ++i)
            vector[pivots[i]] = (rows[i * stride + (free >> 6)] >> (free & 63)) & 1;
        solution.basis.push_back(std::move(vector));
    }
    return solution;
}

namespace detail {
struct MatrixBitReference {
    UInt *word;
    UInt mask;

    operator bool() const {
        return (*word & mask) != 0;
    }

    MatrixBitReference &operator=(bool value) {
        if (value)
            *word |= mask;
        else
            *word &= ~mask;
        return *this;
    }

    template <std::integral T> MatrixBitReference &operator=(T value) {
        return *this = bool(value & 1);
    }

    MatrixBitReference &operator=(const MatrixBitReference &other) {
        return *this = bool(other);
    }
};
}
}
