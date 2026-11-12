#pragma once
#include <cplib/common.hpp>
#include <functional>
#include <iterator>

namespace cplib {
// 添字と要素を述語に渡し、全要素が条件を満たすか判定する。最大N回の述語呼出し。
template <class Range, class Pred> bool allItidx(const Range &s, Pred pred) {
    Int i = 0;
    for (const auto &x : s)
        if (!pred(i++, x))
            return false;
    return true;
}

// 添字と要素を述語に渡し、最初に条件を満たす添字を返す。なければ-1。最大N回の述語呼出し。
template <class Range, class Pred> Int findItIdx(const Range &s, Pred pred) {
    Int i = 0;
    for (const auto &x : s) {
        if (pred(i, x))
            return i;
        ++i;
    }
    return -1;
}

// 添字と要素を述語に渡し、条件を満たすものが存在するか判定する。最大N回の述語呼出し。
template <class Range, class Pred> bool anyItIdx(const Range &s, Pred pred) {
    return findItIdx(s, pred) != -1;
}

// 添字と要素に関数を適用した列を作る。N回の関数呼出しと O(N) の追加操作。
template <class Range, class Op> auto mapItIdx(const Range &s, Op op) {
    using T = std::decay_t<std::invoke_result_t<Op, Int, decltype(*std::begin(s))>>;
    std::vector<T> result;
    result.reserve(std::size(s));
    Int i = 0;
    for (const auto &x : s)
        result.push_back(op(i++, x));
    return result;
}

// 各添字に関数を適用して長さnの列を作る。n回の関数呼出しと O(n) の追加操作。
template <class Init> auto newSeqWithIdx(Int n, Init init) {
    assert(n >= 0);
    using T = std::decay_t<std::invoke_result_t<Init, Int>>;
    std::vector<T> result;
    result.reserve(n);
    for (Int i = 0; i < n; ++i)
        result.push_back(init(i));
    return result;
}
}
