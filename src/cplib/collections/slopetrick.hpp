#pragma once
#include <cplib/utils/constants.hpp>
#include <functional>

namespace cplib {
// 折れ点を二つのヒープで管理。追加 O(log N)、平行移動 O(1)、値の評価 O(N)。
class SlopeTrick {
    std::vector<Int> left_, right_;
    Int minimum_ = 0, left_add_ = 0, right_add_ = 0;

    static void push(std::vector<Int> &h, Int x) {
        h.push_back(x);
        std::push_heap(h.begin(), h.end(), std::greater<Int>());
    }

    static Int pushpop(std::vector<Int> &h, Int x) {
        if (x <= h.front())
            return x;
        std::pop_heap(h.begin(), h.end(), std::greater<Int>());
        Int result = h.back();
        h.back() = x;
        std::push_heap(h.begin(), h.end(), std::greater<Int>());
        return result;
    }

public:
    // 定数関数f(x)=aで初期化する。
    explicit SlopeTrick(Int a = 0) : minimum_(a) {
        clearL();
        clearR();
    }

    // f(x)をmin_{y>=x} f(y)に置き換える。
    void clearL() {
        left_.clear();
        push(left_, INF64 + left_add_);
    }

    // f(x)をmin_{y<=x} f(y)に置き換える。
    void clearR() {
        right_.clear();
        push(right_, INF64 - right_add_);
    }

    Int min() const {
        return minimum_;
    }

    void add_all(Int a) {
        minimum_ += a;
    }

    // f(x)にmax(x-a,0)を加算します。
    void add_x_minus_a(Int a) {
        minimum_ += std::max(-left_.front() + left_add_ - a, Int(0));
        Int x = -pushpop(left_, -a + left_add_) + left_add_;
        push(right_, x - right_add_);
    }

    // f(x)にmax(a-x,0)を加算します。
    void add_a_minus_x(Int a) {
        minimum_ += std::max(a - right_.front() - right_add_, Int(0));
        Int x = pushpop(right_, a - right_add_) + right_add_;
        push(left_, -x + left_add_);
    }

    void add_abs(Int a) {
        add_x_minus_a(a);
        add_a_minus_x(a);
    }

    Int min_index() const {
        return -left_.front() + left_add_;
    }

    void shift(Int a) {
        shift(a, a);
    }

    void shift(Int a, Int b) {
        left_add_ += a;
        right_add_ += b;
    }

    Int get_value(Int x) const {
        Int result = minimum_;
        for (Int y : left_)
            result += std::max(Int(0), -y + left_add_ - x);
        for (Int y : right_)
            result += std::max(Int(0), x - y - right_add_);
        return result;
    }
};

inline SlopeTrick initSlopeTrick(Int a) {
    return SlopeTrick(a);
}
}
