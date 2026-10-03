#pragma once
#include <cplib/collections/private/swag_base.hpp>

namespace cplib {
template <class T> class SWAG : public detail::SwagBase<T> {
    using Base = detail::SwagBase<T>;

public:
    template <class Op> SWAG(Op op, T e) : Base(op, e) {
    }

    void addFirst(const T &x) {
        this->pushtop(x);
    }

    void addLast(const T &x) {
        this->pushbottom(x);
    }

    // 半分ずつ再分配して先頭を削除する。償却O(1)回のモノイド演算。
    T popFirst() {
        if (this->top.empty()) {
            if (this->bottom.empty())
                throw std::out_of_range("SWAG empty");
            std::vector<T> stack;
            std::size_t half = this->bottom.size() / 2;
            for (std::size_t i = 0; i < half; ++i)
                stack.push_back(this->popbottom());
            while (!this->bottom.empty())
                this->pushtop(this->popbottom());
            while (!stack.empty()) {
                this->pushbottom(stack.back());
                stack.pop_back();
            }
        }
        return this->poptop();
    }

    // 半分ずつ再分配して末尾を削除する。償却O(1)回のモノイド演算。
    T popLast() {
        if (this->bottom.empty()) {
            if (this->top.empty())
                throw std::out_of_range("SWAG empty");
            std::vector<T> a, b;
            std::size_t half = this->top.size() / 2;
            for (std::size_t i = 0; i < half; ++i)
                a.push_back(this->poptop());
            while (!this->top.empty())
                b.push_back(this->poptop());
            std::reverse(b.begin(), b.end());
            while (!a.empty()) {
                this->pushtop(a.back());
                a.pop_back();
            }
            while (!b.empty()) {
                this->pushbottom(b.back());
                b.pop_back();
            }
        }
        return this->popbottom();
    }

    std::string to_string() const {
        auto v = this->top;
        std::reverse(v.begin(), v.end());
        v.insert(v.end(), this->bottom.begin(), this->bottom.end());
        return "swag" + detail::swag_string(v);
    }
};

template <class Op, class T> SWAG<T> initSWAG(Op op, T e) {
    return SWAG<T>(op, e);
}

template <class T> void addFirst(SWAG<T> &s, const T &x) {
    s.addFirst(x);
}

template <class T> void addLast(SWAG<T> &s, const T &x) {
    s.addLast(x);
}

template <class T> T popFirst(SWAG<T> &s) {
    return s.popFirst();
}

template <class T> T popLast(SWAG<T> &s) {
    return s.popLast();
}

template <class T> T fold(const SWAG<T> &s) {
    return s.fold();
}

template <class T> Int len(const SWAG<T> &s) {
    return s.len();
}

template <class T> std::string to_string(const SWAG<T> &s) {
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
            s.addLast(v[r++]);
        result.push_back(f(s.fold()) ? Int(v.size()) : r - 1);
        if (s.len())
            s.popFirst();
    }
    return result;
}
}
