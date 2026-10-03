#pragma once
#include <cplib/common.hpp>

namespace cplib {
class Cumsum2D {
    std::vector<std::vector<Int>> B;

public:
    // 二次元累積和を構築する。O(HW)。
    explicit Cumsum2D(std::span<const std::vector<Int>> x) {
        Int h = x.size(), w = h ? x[0].size() : 0;
        B.assign(h + 1, std::vector<Int>(w + 1));
        for (Int i = 1; i <= h; ++i)
            for (Int j = 1; j <= w; ++j)
                B[i][j] = B[i - 1][j] + B[i][j - 1] - B[i - 1][j - 1] + x[i - 1][j - 1];
    }

    // 半開矩形の総和を返す。O(1)。
    Int query(Int il, Int ir, Int jl, Int jr) const {
        return B[ir][jr] - B[ir][jl] - B[il][jr] + B[il][jl];
    }
};

inline Cumsum2D toCumSum2D(std::span<const std::vector<Int>> x) {
    return Cumsum2D(x);
}

inline Int query(const Cumsum2D &s, Int il, Int ir, Int jl, Int jr) {
    return s.query(il, ir, jl, jr);
}
}
