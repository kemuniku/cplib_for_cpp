#pragma once
#include <cplib/collections/private/convex_hull_trick_impl.hpp>

namespace cplib {
class ConvexHullTrickMonotoneSlope {
    CHTMonotoneHull hull_;

public:
    // 傾きが単調な最小値CHTを初期化します。increasing=falseは広義単調減少、trueは広義単調増加です。O(1)。
    explicit ConvexHullTrickMonotoneSlope(bool increasing = false) : hull_(increasing) {
    }

    // ax+bを追加します。傾きは指定した向きに単調である必要があります。償却O(1)。
    void add_line(Int a, Int b) {
        hull_.chtAddLine(a, b);
    }

    // 任意座標での最小値を二分探索。O(log N)。
    // 任意の整数座標xでの最小値を返します。空の場合はassert。O(log N)。
    Int get_min(Int x) const {
        assert(!hull_.lines.empty());
        std::size_t l = 0, r = hull_.lines.size() - 1;
        while (l < r) {
            auto m = l + (r - l) / 2;
            if (chtValue(hull_.lines[m], x) >= chtValue(hull_.lines[m + 1], x))
                l = m + 1;
            else
                r = m;
        }
        return chtAnswer(chtValue(hull_.lines[l], x));
    }
};

inline ConvexHullTrickMonotoneSlope initConvexHullTrickMonotoneSlope(bool slopeIncreasing = false) {
    return ConvexHullTrickMonotoneSlope(slopeIncreasing);
}
}
