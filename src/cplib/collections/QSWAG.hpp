#pragma once
#include <cplib/collections/private/swag_base.hpp>

namespace cplib {
template <class T> class QSWAG : public detail::SwagBase<T> {
    using Base = detail::SwagBase<T>;

public:
    // モノイド演算と単位元を指定して空の列を構築する。O(1)。
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

    // 全要素を文字列化する。要素数をnとしてO(n)と各要素の文字列化の計算量。
    std::string to_string() const {
        auto t = this->top;
        std::reverse(t.begin(), t.end());
        return detail::swag_string(t) + detail::swag_string(this->bottom);
    }
};

// モノイド演算と単位元を指定して空の列を生成する。O(1)。
template <class Op, class T> QSWAG<T> initQSWAG(Op op, T e) {
    return QSWAG<T>(op, e);
}

// 末尾に追加する。償却O(1)回のモノイド演算。
template <class T> void push(QSWAG<T> &s, const T &x) {
    s.push(x);
}

// 先頭から取り出す。償却O(1)回のモノイド演算。
template <class T> T pop(QSWAG<T> &s) {
    return s.pop();
}

// 全要素の総積を返す。モノイド演算を1回行う。
template <class T> T fold(const QSWAG<T> &s) {
    return s.fold();
}

// 全要素を文字列化する。要素数をnとしてO(n)と各要素の文字列化の計算量。
template <class T> std::string to_string(const QSWAG<T> &s) {
    return s.to_string();
}

// 単調な述語を満たす各左端の最大右端を求める。O(n)回のモノイド演算。
template <class Range, class Op, class T, class Predicate>
std::vector<Int> get_maxrights(const Range &v, Op op, T e, Predicate f) {
    assert(f(e));
    auto s = initQSWAG(op, e);
    Int r = 0;
    std::vector<Int> result;
    for (Int l = 0; l < Int(v.size()); ++l) {
        if (l > r)
            r = l;
        while (r < Int(v.size()) && f(s.fold()))
            s.push(v[r++]);
        result.push_back(f(s.fold()) ? Int(v.size()) : r - 1);
        if (s.size())
            s.pop();
    }
    return result;
}
}
