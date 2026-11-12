#pragma once
#include <cplib/collections/private/dynamic_retroactive_impl.hpp>

namespace cplib {
template <class K, class T, class S = void>
using DynamicRetroactivePriorityQueue = detail::RetroactiveAVL<K, T, S, false>;

// 空の動的時刻付き遡及優先度付きキューを構築する。O(1)。
template <class K, class T> auto initDynamicRetroactivePriorityQueue(SortOrder order = Ascending) {
    return DynamicRetroactivePriorityQueue<K, T>(order);
}

// 残存要素のモノイド集約に対応する空のキューを構築する。演算・値のコピーがO(1)ならO(1)。
template <class K, class T, class S, class Op, class Lift>
auto initDynamicRetroactivePriorityQueue(Op op, S e, Lift lift, SortOrder order = Ascending) {
    return DynamicRetroactivePriorityQueue<K, T, S>(op, e, lift, order);
}
}
