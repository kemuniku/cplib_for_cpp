#pragma once
#include <cplib/common.hpp>
#include <cplib/utils/backwards_index.hpp>
#include <unordered_map>

namespace cplib {
template <class T> class StaticRangeCount {
    Int length_;
    std::unordered_map<T, std::vector<Int>> t;

public:
    explicit StaticRangeCount(std::span<const T> v) : length_(v.size()) {
        for (Int i = 0; i < Int(v.size()); ++i)
            t[v[i]].push_back(i);
    }

    Int size() const {
        return length_;
    }

    // 構築は期待O(N)、区間内の出現回数はO(log N)。
    Int count(Int l, Int r, const T &x) const {
        auto it = t.find(x);
        if (it == t.end())
            return 0;
        const auto &a = it->second;
        return std::lower_bound(a.begin(), a.end(), r) - std::lower_bound(a.begin(), a.end(), l);
    }

    Int count(ClosedSlice<Int, Int> range, const T &x) const {
        auto it = t.find(x);
        if (it == t.end())
            return 0;
        const auto &a = it->second;
        return std::upper_bound(a.begin(), a.end(), range.b) -
               std::lower_bound(a.begin(), a.end(), range.a);
    }
};

template <class T> auto initStaticRangeCount(std::span<const T> v) {
    return StaticRangeCount<T>(v);
}

template <class T> auto initStaticRangeCount(const std::vector<T> &v) {
    return StaticRangeCount<T>(v);
}

template <class T> Int count(const StaticRangeCount<T> &s, Int l, Int r, const T &x) {
    return s.count(l, r, x);
}

template <class T>
Int count(const StaticRangeCount<T> &s, ClosedSlice<Int, Int> range, const T &x) {
    return s.count(range, x);
}
}
