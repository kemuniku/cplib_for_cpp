#pragma once
#include <cplib/collections/QSWAG.hpp>

namespace cplib {
template <class T> class SWAG : public detail::SwagBase<T> {
    using Base = detail::SwagBase<T>;

public:
    // モノイド演算と単位元を指定して空の列を構築する。O(1)。
    template <class Op> SWAG(Op op, T e) : Base(op, e) {
    }

    // 先頭に追加する。償却O(1)回のモノイド演算。
    void addFirst(const T &x) {
        this->pushtop(x);
    }

    // 末尾に追加する。償却O(1)回のモノイド演算。
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

    // 全要素を文字列化する。要素数をnとしてO(n)と各要素の文字列化の計算量。
    std::string to_string() const {
        auto v = this->top;
        std::reverse(v.begin(), v.end());
        v.insert(v.end(), this->bottom.begin(), this->bottom.end());
        return "swag" + detail::swag_string(v);
    }
};

// モノイド演算と単位元を指定して空の列を生成する。O(1)。
template <class Op, class T> SWAG<T> initSWAG(Op op, T e) {
    return SWAG<T>(op, e);
}

// 演算と単位元から空のSWAGを作る。計算量 O(1)。
template <class Op, class T> auto newSwagWith(Op op, T e) {
    return initSWAG(op, e);
}

// 先頭に追加する。償却O(1)回のモノイド演算。
template <class T> void addFirst(SWAG<T> &s, const T &x) {
    s.addFirst(x);
}

// 末尾に追加する。償却O(1)回のモノイド演算。
template <class T> void addLast(SWAG<T> &s, const T &x) {
    s.addLast(x);
}

// 先頭から取り出す。償却O(1)回のモノイド演算。
template <class T> T popFirst(SWAG<T> &s) {
    return s.popFirst();
}

// 末尾から取り出す。償却O(1)回のモノイド演算。
template <class T> T popLast(SWAG<T> &s) {
    return s.popLast();
}

// 全要素の総積を返す。モノイド演算を1回行う。
template <class T> T fold(const SWAG<T> &s) {
    return s.fold();
}

// 全要素を文字列化する。要素数をnとしてO(n)と各要素の文字列化の計算量。
template <class T> std::string to_string(const SWAG<T> &s) {
    return s.to_string();
}

}
