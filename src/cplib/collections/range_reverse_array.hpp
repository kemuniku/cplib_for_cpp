#pragma once
#include <cplib/collections/private/range_reverse_treap.hpp>

namespace cplib {
template <class T>
using RangeReverseArray = detail::RangeReverseTreap<T, detail::NoTreapAction, false, false>;

// 初期列から構築する。要素数をnとしてO(n+1)回の演算・作用・合成。
template <class T> auto initRangeReverseArray(std::span<const T> values) {
    return RangeReverseArray<T>(values, {}, T{}, {}, {}, detail::NoTreapAction{});
}

// 初期列から構築する。要素数をnとしてO(n+1)回の演算・作用・合成。
template <class T> auto initRangeReverseArray(const std::vector<T> &values) {
    return RangeReverseArray<T>(values, {}, T{}, {}, {}, detail::NoTreapAction{});
}

// 初期列から構築する。要素数をnとしてO(n+1)回の演算・作用・合成。
template <class V> auto toRangeReverseArray(const V &values) {
    return initRangeReverseArray(values);
}
}
