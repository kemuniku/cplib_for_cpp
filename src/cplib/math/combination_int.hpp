#pragma once
#include <cplib/utils/constants.hpp>
namespace cplib {
// nCrをlimitで打ち切って返す。O(min(r,n-r))、既定limitでは約32回以下。
inline Int ncr_int(Int n, Int r, Int limit=INF64) {
    if (limit<=0) return limit;
    if (n<0 || r<0 || n<r) return 0;
    Int rr=std::min(r,n-r); __int128 res=1;
    for (Int i=1;i<=rr;++i) { res=res*(n-rr+i)/i; if (res>=limit) return limit; }
    return static_cast<Int>(res);
}
}
