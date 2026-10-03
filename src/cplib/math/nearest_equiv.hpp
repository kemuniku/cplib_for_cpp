#pragma once
#include <cplib/common.hpp>

namespace cplib {
// y≡x (mod m), y>=l を満たす最小値を返す。O(1)。
inline Int nearest_equiv(Int xx, Int ll, Int mm) {
    assert(mm != 0);
    __int128 x = xx, l = ll, m = mm;
    if (m < 0)
        m = -m;
    return static_cast<Int>(x < l ? x + (l - x + m - 1) / m * m : x - (x - l) / m * m);
}
}
