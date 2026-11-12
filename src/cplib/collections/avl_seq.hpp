#pragma once
#include <cplib/collections/private/avl_sequence.hpp>

namespace cplib {
template <class T> using AvlSeq = detail::AvlSequence<T, detail::NoAvlAction, false, false>;
// 配列からAVL列を構築する。時間・領域 O(N)。
template <class T> auto initAvlSeq(const std::vector<T> &values) {
    return AvlSeq<T>(values, {}, T{}, {}, {}, {});
}
// 配列をAVL列に変換する。時間・領域 O(N)。
template <class T> auto toAvlSeq(const std::vector<T> &values) {
    return initAvlSeq(values);
}
}
