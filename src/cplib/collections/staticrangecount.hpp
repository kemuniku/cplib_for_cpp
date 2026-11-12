#pragma once
#include <cplib/common.hpp>
#include <cplib/utils/backwards_index.hpp>
#include <unordered_map>

namespace cplib {
template <class T> class StaticRangeCount {
    Int length_;
    std::unordered_map<T, std::vector<Int>> t;

public:
    // 配列から出現位置を構築する。要素数nとして期待O(n)、最悪O(n²)。
    explicit StaticRangeCount(std::span<const T> v) : length_(v.size()) {
        for (Int i = 0; i < Int(v.size()); ++i)
            t[v[i]].push_back(i);
    }

    // 要素数を返す。O(1)。
    Int size() const {
        return length_;
    }

    // [l,r)内のxの出現回数を返す。期待O(log(n+1))、ハッシュ衝突時は最悪O(n)。
    Int count(Int l, Int r, const T &x) const {
        auto it = t.find(x);
        if (it == t.end())
            return 0;
        const auto &a = it->second;
        return std::lower_bound(a.begin(), a.end(), r) - std::lower_bound(a.begin(), a.end(), l);
    }

    // 両端を含む区間内のxの出現回数を返す。期待O(log(n+1))、ハッシュ衝突時は最悪O(n)。
    Int count(ClosedSlice<Int, Int> range, const T &x) const {
        auto it = t.find(x);
        if (it == t.end())
            return 0;
        const auto &a = it->second;
        return std::upper_bound(a.begin(), a.end(), range.b) -
               std::lower_bound(a.begin(), a.end(), range.a);
    }
};

// 配列から構築する。要素数nとして期待O(n)、最悪O(n²)。
template <class T> auto initStaticRangeCount(std::span<const T> v) {
    return StaticRangeCount<T>(v);
}

// 配列から構築する。要素数nとして期待O(n)、最悪O(n²)。
template <class T> auto initStaticRangeCount(const std::vector<T> &v) {
    return StaticRangeCount<T>(v);
}

// [l,r)内のxの出現回数を返す。期待O(log(n+1))、最悪O(n)。
template <class T> Int count(const StaticRangeCount<T> &s, Int l, Int r, const T &x) {
    return s.count(l, r, x);
}

// 両端を含む区間内のxの出現回数を返す。期待O(log(n+1))、最悪O(n)。
template <class T>
Int count(const StaticRangeCount<T> &s, ClosedSlice<Int, Int> range, const T &x) {
    return s.count(range, x);
}
}
