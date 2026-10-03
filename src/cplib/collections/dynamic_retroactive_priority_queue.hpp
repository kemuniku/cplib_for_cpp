#pragma once
#include <cplib/collections/private/dynamic_retroactive_impl.hpp>

namespace cplib {
template <class K, class T, class S = void>
using DynamicRetroactivePriorityQueue = detail::RetroactiveAVL<K, T, S, false>;

template <class K, class T> auto initDynamicRetroactivePriorityQueue(SortOrder order = Ascending) {
    return DynamicRetroactivePriorityQueue<K, T>(order);
}

template <class K, class T, class S, class Op, class Lift>
auto initDynamicRetroactivePriorityQueue(Op op, S e, Lift lift, SortOrder order = Ascending) {
    return DynamicRetroactivePriorityQueue<K, T, S>(op, e, lift, order);
}
}
