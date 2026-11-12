#pragma once
#include <cplib/common.hpp>

namespace cplib {
// 各判定が成立する最小操作数を返す。初期成立は0、不成立はM+1。判定は単調であること。
// 各巡回の先頭でresetしapply(0),apply(1),...を適用する。Q=0なら呼び出さない。
// コールバックを除く時間 O((M+Q) log(M+2))、領域 O(M+Q)。終了時の状態は復元しない。
template <class Reset, class Apply, class Check>
std::vector<Int> parallelBinarySearch(Int updateCount, Int queryCount,
                                     Reset reset, Apply apply, Check check) {
    assert(0 <= updateCount && updateCount < std::numeric_limits<Int>::max());
    assert(queryCount >= 0);
    std::vector<Int> upper(queryCount, updateCount + 1), lower(queryCount);
    if (queryCount == 0)
        return upper;
    std::vector<Int> heads(updateCount + 1), next(queryCount);
    while (true) {
        std::fill(heads.begin(), heads.end(), -1);
        Int last = -1;
        for (Int q = 0; q < queryCount; ++q) {
            if (lower[q] == upper[q])
                continue;
            Int mid = lower[q] + (upper[q] - lower[q]) / 2;
            next[q] = heads[mid];
            heads[mid] = q;
            last = std::max(last, mid);
        }
        if (last < 0)
            break;
        reset();
        for (Int count = 0; count <= last; ++count) {
            if (count > 0)
                apply(count - 1);
            for (Int q = heads[count]; q >= 0; q = next[q]) {
                if (check(q))
                    upper[q] = count;
                else
                    lower[q] = count + 1;
            }
        }
    }
    return upper;
}
}
