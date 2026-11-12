#pragma once
#include <cplib/collections/private/range_reverse_treap.hpp>

namespace cplib {
template <class T>
using RangeReverseArrayMonoid = detail::RangeReverseTreap<T, detail::NoTreapAction, true, false>;

// 初期列から構築する。要素数をnとしてO(n+1)回の演算・作用・合成。
template <class T, class V, class Op>
auto initRangeReverseArrayMonoid(const V &values, Op op, T e) {
    return RangeReverseArrayMonoid<T>(values, op, e, {}, {}, detail::NoTreapAction{});
}

// 初期列から構築する。要素数をnとしてO(n+1)回の演算・作用・合成。
template <class V, class Op, class T> auto toRangeReverseArrayMonoid(const V &values, Op op, T e) {
    return initRangeReverseArrayMonoid(values, op, e);
}

// 初期列から構築する。要素数をnとしてO(n+1)回の演算・作用・合成。
template <class V, class Op, class T>
auto newRangeReverseArrayMonoidWith(const V &values, Op op, T e) {
    return initRangeReverseArrayMonoid(values, op, e);
}
}
