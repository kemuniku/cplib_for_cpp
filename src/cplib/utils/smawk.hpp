#pragma once
#include <cplib/common.hpp>

namespace cplib {
// 全単調行列の左端最小列をO(height+width)回の比較で求める。幅0は-1。
// better(row, oldCol, newCol)は新しい列が真に小さいときtrue。同値なら左端を返す。
template <class Better> std::vector<Int> smawk(Int height, Int width, Better better) {
    assert(height >= 0 && width >= 0);
    std::vector<Int> answer(height, width ? 0 : -1);
    if (!height || !width)
        return answer;
    Int capacity = width;
    for (Int count = height; count > 0; count /= 2)
        capacity += std::min(count, width);
    std::vector<Int> columns(capacity);
    std::iota(columns.begin(), columns.begin() + width, 0);
    auto solve = [&](auto &&self, Int first, Int step, Int count, Int offset, Int w) -> void {
        Int reduced = offset + w, size = 0;
        for (Int p = offset; p < reduced; ++p) {
            Int col = columns[p];
            while (size > 0 && better(first + (size - 1) * step, columns[reduced + size - 1], col))
                --size;
            if (size < count)
                columns[reduced + size++] = col;
        }
        if (count > 1)
            self(self, first + step, step * 2, count / 2, reduced, size);
        Int left = 0;
        for (Int i = 0; i < count; i += 2) {
            Int row = first + i * step,
                right = i + 1 < count ? answer[row + step] : columns[reduced + size - 1],
                best = columns[reduced + left];
            while (columns[reduced + left] < right) {
                Int col = columns[reduced + (++left)];
                if (better(row, best, col))
                    best = col;
            }
            answer[row] = best;
        }
    };
    solve(solve, 0, 1, height, 0, width);
    return answer;
}
}
