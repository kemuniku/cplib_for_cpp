#pragma once
#include <cplib/collections/private/swag_base.hpp>

namespace cplib {
template <class T> class QSWAG : public detail::SwagBase<T> {
    using Base = detail::SwagBase<T>;

public:
    template <class Op> QSWAG(Op op, T e) : Base(op, e) {
    }

    // 末尾に追加する。償却O(1)回のモノイド演算。
    void push(const T &x) {
        this->pushbottom(x);
    }

    // 先頭から削除する。償却O(1)回のモノイド演算。
    T pop() {
        if (this->top.empty()) {
            if (this->bottom.empty())
                throw std::out_of_range("QSWAG empty");
            while (!this->bottom.empty())
                this->pushtop(this->popbottom());
        }
        return this->poptop();
    }

    std::string to_string() const {
        auto t = this->top;
        std::reverse(t.begin(), t.end());
        return detail::swag_string(t) + detail::swag_string(this->bottom);
    }
};

template <class Op, class T> QSWAG<T> initSWAG(Op op, T e) {
    return QSWAG<T>(op, e);
}

template <class T> void push(QSWAG<T> &s, const T &x) {
    s.push(x);
}

template <class T> T pop(QSWAG<T> &s) {
    return s.pop();
}

template <class T> T fold(const QSWAG<T> &s) {
    return s.fold();
}

template <class T> Int len(const QSWAG<T> &s) {
    return s.len();
}

template <class T> std::string to_string(const QSWAG<T> &s) {
    return s.to_string();
}

// 単調な述語を満たす各左端の最大右端を求める。O(n)回のモノイド演算。
template <class Range, class Op, class T, class Predicate>
std::vector<Int> get_maxrights(const Range &v, Op op, T e, Predicate f) {
    assert(f(e));
    auto s = initSWAG(op, e);
    Int r = 0;
    std::vector<Int> result;
    for (Int l = 0; l < Int(v.size()); ++l) {
        if (l > r)
            r = l;
        while (r < Int(v.size()) && f(s.fold()))
            s.push(v[r++]);
        result.push_back(f(s.fold()) ? Int(v.size()) : r - 1);
        if (s.len())
            s.pop();
    }
    return result;
}
}
