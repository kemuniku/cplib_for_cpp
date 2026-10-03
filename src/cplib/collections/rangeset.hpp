#pragma once
#include <cplib/collections/avlset.hpp>
#include <tuple>

namespace cplib {
// 区間の値を AVL 集合で管理。更新 O((消す区間数+1) log N)、取得 O(log N)。
template <class T> struct RangeSet {
    using Segment = std::tuple<Int, Int, T>;
    AvlSortedMultiSet<Segment> st;
    T default_value;

    explicit RangeSet(T value) : default_value(value) {
        st.incl({std::numeric_limits<Int>::min(), std::numeric_limits<Int>::max(), value});
    }

    void update(Int l, Int r, const T &value) {
        for (;;) {
            auto x = st.ge({l, std::numeric_limits<Int>::min(), default_value});
            if (!x)
                break;
            auto [a, b, c] = *x;
            if (b <= r)
                st.excl(*x);
            else if (a < r) {
                st.excl(*x);
                if (c == value)
                    r = b;
                else
                    st.incl({r, b, c});
                break;
            } else {
                if (a == r && c == value) {
                    st.excl(*x);
                    r = b;
                }
                break;
            }
        }
        auto x = st.le({l, std::numeric_limits<Int>::max(), default_value});
        if (x) {
            auto [a, b, c] = *x;
            if (r < b) {
                if (c != value) {
                    st.excl(*x);
                    st.incl({a, l, c});
                    st.incl({r, b, c});
                } else
                    return;
            } else if (l < b) {
                st.excl(*x);
                if (c != value)
                    st.incl({a, l, c});
                else
                    l = a;
            } else if (l == b && c == value) {
                st.excl(*x);
                l = a;
            }
        }
        st.incl({l, r, value});
    }

    Segment get_segment(Int x) const {
        return st.le({x, std::numeric_limits<Int>::max(), default_value}).value();
    }
};

template <class T> auto initRangeSet(T value) {
    return RangeSet<T>(value);
}
}
