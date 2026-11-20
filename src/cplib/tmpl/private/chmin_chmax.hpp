#pragma once

namespace cplib {
namespace detail {
namespace template_support {

// yが小さければxを更新し、更新したか返す。比較を1回、必要な場合に代入を1回行う。
template <class T> inline bool chmin(T &x, const T &y) {
    if (x > y) {
        x = y;
        return true;
    }
    return false;
}

// yが大きければxを更新し、更新したか返す。比較を1回、必要な場合に代入を1回行う。
template <class T> inline bool chmax(T &x, const T &y) {
    if (x < y) {
        x = y;
        return true;
    }
    return false;
}

} // namespace template_support
} // namespace detail
} // namespace cplib
