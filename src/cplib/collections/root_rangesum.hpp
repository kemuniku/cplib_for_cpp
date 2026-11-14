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
    // 配列をブロック長bsizeで平方分割する。省略時は約√n。要素数nとしてO(n)。
    explicit RootRangeSum(std::vector<T> v, Int bsize = 0, T e = T(0))
        : blocksize_(bsize > 0 ? bsize : std::max(Int(std::sqrt(v.size())), Int(1))),
          arr_(std::move(v)), blocks_((arr_.size() + blocksize_ - 1) / blocksize_, e), e_(e) {
        for (Int i = 0; i < size(); ++i)
            blocks_[i / blocksize_] = blocks_[i / blocksize_] + arr_[i];
    }

    // 要素数を返す。O(1)。
    Int size() const {
        return arr_.size();
    }

    // i番目の値を設定する。O(1)回の要素演算。
    void update(Int i, const T &value) {
        assert(0 <= i && i < size());
        blocks_[i / blocksize_] = blocks_[i / blocksize_] + value - arr_[i];
        arr_[i] = value;
    }

    // 要素数n、ブロック長BとしてO(B+n/B)回の要素演算。
    // 半開区間[l,r)の演算結果を返します。
    T get(Int l, Int r) const {
        assert(0 <= l && l <= r && r <= size());
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

    // 両端を含む区間の和を返す。要素数n、ブロック長BとしてO(B+n/B)。
    template <class L, class R> T get(ClosedSlice<L, R> s) const {
        return get(resolve_index(size(), s.a), resolve_index(size(), s.b) + 1);
    }

    // 両端を含む区間の和を返す。要素数n、ブロック長BとしてO(B+n/B)。
    template <class L, class R> T operator[](ClosedSlice<L, R> s) const {
        return get(s);
    }

    // i番目の値を返す。O(1)。
    T operator[](Int i) const {
        assert(0 <= i && i < size());
        return arr_[i];
    }

    struct Reference {
        RootRangeSum *owner;
        Int i;

        // 対象の値を取得する。O(1)。
        operator T() const {
            return owner->arr_[i];
        }

        // 対象の値を設定する。O(1)回の要素演算。
        Reference &operator=(const T &v) {
            owner->update(i, v);
            return *this;
        }

        // 他の位置の値を読み取って代入する。O(1)回の要素演算。
        Reference &operator=(const Reference &v) {
            return *this = T(v);
        }
    };

    // i番目を読み書きする参照を返す。O(1)。
    Reference operator[](Int i) {
        assert(0 <= i && i < size());
        return {this, i};
    }

    CPLIB_BACKWARDS_INDEX_OVERLOADS
    // f(和[l,r))が真となる最大のrを返す。fは単位元で真かつ単調。要素数n、ブロック長BとしてO(B+n/B)回の演算と述語呼出し。
    template <class F> Int max_right(Int l, F f) const {
        assert(0 <= l && l <= size() && f(e_));
        if (l == size())
            return size();
        T sm = e_;
        Int bl = l / blocksize_;
        for (Int i = l; i < std::min((bl + 1) * blocksize_, size()); ++i) {
            if (!f(sm + arr_[i]))
                return i;
            sm = sm + arr_[i];
        }
        for (Int b = bl + 1; b < Int(blocks_.size()); ++b) {
            if (!f(sm + blocks_[b])) {
                for (Int i = b * blocksize_; i < (b + 1) * blocksize_; ++i) {
                    if (i >= size())
                        return i;
                    if (!f(sm + arr_[i]))
                        return i;
                    sm = sm + arr_[i];
                }
            } else
                sm = sm + blocks_[b];
        }
        return size();
    }

    // f(和[l,r))が真となる最小のlを返す。fは単位元で真かつ単調。要素数n、ブロック長BとしてO(B+n/B)回の演算と述語呼出し。
    template <class F> Int min_left(Int r, F f) const {
        assert(0 <= r && r <= size() && f(e_));
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
                for (Int i = std::min((b + 1) * blocksize_, size()); i-- > b * blocksize_;) {
                    if (!f(sm + arr_[i]))
                        return i + 1;
                    sm = sm + arr_[i];
                }
            } else
                sm = sm + blocks_[b];
        }
        return 0;
    }

    // 全要素を文字列化する。要素数nとしてO(n+出力長)と各要素の文字列化の計算量。
    std::string to_string() const {
        std::ostringstream out;
        out << "@[";
        for (Int i = 0; i < size(); ++i) {
            if (i)
                out << ", ";
            out << arr_[i];
        }
        return out.str() + "]";
    }
};

// 配列から平方分割区間和を構築する。要素数nとしてO(n)。
template <class T> auto initrangesum(const std::vector<T> &v, Int bsize = 0, T e = T(0)) {
    return RootRangeSum<T>(v, bsize, e);
}
}
