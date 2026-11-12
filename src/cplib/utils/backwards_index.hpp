#pragma once
#include <cplib/common.hpp>

namespace cplib {
// Nimの^kに対応する末尾からの添字。^1は最後の要素。
struct BackwardsIndex {
    Int value;

    // 末尾からk番目を表す添字を作る。計算量 O(1)。
    explicit constexpr BackwardsIndex(Int k) : value(k) {
    }
};

// 末尾からk番目を表す添字を返す。計算量 O(1)。
constexpr BackwardsIndex from_end(Int k) {
    return BackwardsIndex(k);
}

// 列長から末尾基準の添字を先頭基準に変換する。計算量 O(1)。
constexpr Int resolve_index(Int length, BackwardsIndex index) {
    return length - index.value;
}

// 先頭基準の添字をそのまま返す。計算量 O(1)。
constexpr Int resolve_index(Int, Int index) {
    return index;
}

template <class Left = Int, class Right = Int> struct ClosedSlice {
    Left a;
    Right b;
};
template <class Left, class Right> ClosedSlice(Left, Right) -> ClosedSlice<Left, Right>;

// 両端を含む区間の端点を保持する。計算量 O(1)。
template <class Left, class Right> constexpr auto closed_slice(Left a, Right b) {
    return ClosedSlice<Left, Right>{a, b};
}
}

// 整数添字版へ転送するオーバーロードをクラス内に追加する。追加計算量len()+O(1)。
#define CPLIB_BACKWARDS_INDEX_OVERLOADS                                                            \
    decltype(auto) operator[](::cplib::BackwardsIndex index) {                                     \
        return (*this)[len() - index.value];                                                       \
    }                                                                                              \
    decltype(auto) operator[](::cplib::BackwardsIndex index) const {                               \
        return (*this)[len() - index.value];                                                       \
    }
