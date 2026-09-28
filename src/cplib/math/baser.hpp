#pragma once
#include <cplib/common.hpp>
namespace cplib {
// 下位桁から基数変換する。size省略時の0は空列。O(桁数)。
inline std::vector<Int> baser(Int num, Int base, Int size=-1) {
    std::vector<Int> result;
    if (size==-1) { while (num>0) { result.push_back(num%base); num/=base; } }
    else for (Int i=0;i<size;++i) { result.push_back(num%base); num/=base; }
    return result;
}
}
