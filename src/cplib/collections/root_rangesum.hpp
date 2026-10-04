#pragma once
#include <cplib/utils/backwards_index.hpp>
#include <cmath>
#include <sstream>

namespace cplib {
// 平方分割の区間和。構築 O(N)、点更新 O(1)、和と境界検索 O(B+N/B)。
template <class T> class RootRangeSum {
    Int blocksize_;
    std::vector<T> arr_, blocks_;
    T e_;

public:
    explicit RootRangeSum(std::vector<T> v, Int bsize = 0, T e = T(0))
        : blocksize_(bsize > 0 ? bsize : std::max(Int(std::sqrt(v.size())), Int(1))),
          arr_(std::move(v)), blocks_((arr_.size() + blocksize_ - 1) / blocksize_, e), e_(e) {
        for (Int i = 0; i < len(); ++i)
            blocks_[i / blocksize_] = blocks_[i / blocksize_] + arr_[i];
    }

    Int len() const {
        return arr_.size();
    }

    Int size() const {
        return len();
    }

    void update(Int i, const T &value) {
        assert(0 <= i && i < len());
        blocks_[i / blocksize_] = blocks_[i / blocksize_] + value - arr_[i];
        arr_[i] = value;
    }

    // 半開区間[l,r)の演算結果を返します。
    T get(Int l, Int r) const {
        assert(0 <= l && l <= r && r <= len());
        T result = e_;
        Int bl = l / blocksize_, br = r / blocksize_;
        if (bl == br) {
            for (Int i = l; i < r; ++i)
                result = result + arr_[i];
            return result;
        }
        for (Int i = l; i < (bl + 1) * blocksize_; ++i)
            result = result + arr_[i];
        for (Int b = bl + 1; b < br; ++b)
            result = result + blocks_[b];
        for (Int i = br * blocksize_; i < r; ++i)
            result = result + arr_[i];
        return result;
    }

    template <class L, class R> T get(ClosedSlice<L, R> s) const {
        return get(resolve_index(len(), s.a), resolve_index(len(), s.b) + 1);
    }

    template <class L, class R> T operator[](ClosedSlice<L, R> s) const {
        return get(s);
    }

    T operator[](Int i) const {
        assert(0 <= i && i < len());
        return arr_[i];
    }

    struct Reference {
        RootRangeSum *owner;
        Int i;

        operator T() const {
            return owner->arr_[i];
        }

        Reference &operator=(const T &v) {
            owner->update(i, v);
            return *this;
        }

        Reference &operator=(const Reference &v) {
            return *this = T(v);
        }
    };

    Reference operator[](Int i) {
        assert(0 <= i && i < len());
        return {this, i};
    }

    CPLIB_BACKWARDS_INDEX_OVERLOADS
    template <class F> Int max_right(Int l, F f) const {
        assert(0 <= l && l <= len() && f(e_));
        if (l == len())
            return len();
        T sm = e_;
        Int bl = l / blocksize_;
        for (Int i = l; i < std::min((bl + 1) * blocksize_, len()); ++i) {
            if (!f(sm + arr_[i]))
                return i;
            sm = sm + arr_[i];
        }
        for (Int b = bl + 1; b < Int(blocks_.size()); ++b) {
            if (!f(sm + blocks_[b])) {
                for (Int i = b * blocksize_; i < (b + 1) * blocksize_; ++i) {
                    if (i >= len())
                        return i;
                    if (!f(sm + arr_[i]))
                        return i;
                    sm = sm + arr_[i];
                }
            } else
                sm = sm + blocks_[b];
        }
        return len();
    }

    template <class F> Int min_left(Int r, F f) const {
        assert(0 <= r && r <= len() && f(e_));
        if (!r)
            return 0;
        T sm = e_;
        Int br = (r - 1) / blocksize_;
        for (Int i = r; i-- > br * blocksize_;) {
            if (!f(sm + arr_[i]))
                return i + 1;
            sm = sm + arr_[i];
        }
        for (Int b = br; b-- > 0;) {
            if (!f(sm + blocks_[b])) {
                for (Int i = std::min((b + 1) * blocksize_, len()); i-- > b * blocksize_;) {
                    if (!f(sm + arr_[i]))
                        return i + 1;
                    sm = sm + arr_[i];
                }
            } else
                sm = sm + blocks_[b];
        }
        return 0;
    }

    std::string to_string() const {
        std::ostringstream out;
        out << "@[";
        for (Int i = 0; i < len(); ++i) {
            if (i)
                out << ", ";
            out << arr_[i];
        }
        return out.str() + "]";
    }
};

template <class T> auto initrangesum(const std::vector<T> &v, Int bsize = 0, T e = T(0)) {
    return RootRangeSum<T>(v, bsize, e);
}
}
