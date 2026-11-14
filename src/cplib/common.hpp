#pragma once
#include <algorithm>
#include <array>
#include <bit>
#include <cassert>
#include <cstdint>
#include <limits>
#include <numeric>
#include <span>
#include <string>
#include <utility>
#include <vector>

namespace cplib {
using Int = std::int64_t;
using UInt = std::uint64_t;

// size()を持つ列の長さを返す。計算量はsize()と同じ。
// 括弧付きの宣言で、テンプレート側のlenマクロとも併用できる。
template <class T>
    requires requires(const T &value) { value.size(); }
constexpr Int (len)(const T &value) {
    return static_cast<Int>(value.size());
}
}
