#pragma once
#include <cplib/common.hpp>
namespace cplib {
// 128ビット中間値で加算の剰余を求める。O(1)。
inline Int add(Int a, Int b, Int m) { return (static_cast<__int128>(a) + b) % m; }
// 128ビット中間値で乗算の剰余を求める。O(1)。
inline Int mul(Int a, Int b, Int m) { return static_cast<__int128>(a) * b % m; }
}
