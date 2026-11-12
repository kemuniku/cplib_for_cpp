#pragma once
#include <cplib/collections/private/avl_sequence.hpp>

namespace cplib {
template <class T> using AvlMonoid = detail::AvlSequence<T, detail::NoAvlAction, true, false>;
// 結合演算opと単位元eを持つAVL列を作る。時間・領域 O(N)。
template <class T, class Op> auto initAvlMonoid(const std::vector<T> &values, Op op, T e) {
    return AvlMonoid<T>(values, op, e, {}, {}, {});
}
// 配列を集約付きAVL列に変換する。時間・領域 O(N)。
template <class T, class Op> auto toAvlMonoid(const std::vector<T> &values, Op op, T e) {
    return initAvlMonoid(values, op, e);
}
// 演算を指定して集約付きAVL列を作る。時間・領域 O(N)。
template <class T, class Op> auto newAvlMonoidWith(const std::vector<T> &values, Op op, T e) {
    return initAvlMonoid(values, op, e);
}
}
