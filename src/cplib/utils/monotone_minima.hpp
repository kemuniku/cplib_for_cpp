#pragma once
#include <cplib/common.hpp>

namespace cplib {
// 単調な左端最小列を探索する。O(height+width log(height+1))回の比較。
// better(row, oldCol, newCol)は新しい列が真に小さいときtrue。幅0なら各行に-1。
template <class Better> std::vector<Int> monotoneMinima(Int height, Int width, Better better) {
    assert(height >= 0 && width >= 0);
    std::vector<Int> answer(height, -1);
    if (!height || !width)
        return answer;
    auto solve = [&](auto &&self, Int top, Int bottom, Int left, Int right) -> void {
        if (top >= bottom)
            return;
        Int row = (top + bottom) / 2, best = left;
        for (Int col = left + 1; col <= right; ++col)
            if (better(row, best, col))
                best = col;
        answer[row] = best;
        self(self, top, row, left, best);
        self(self, row + 1, bottom, best, right);
    };
    solve(solve, 0, height, 0, width - 1);
    return answer;
}
}
