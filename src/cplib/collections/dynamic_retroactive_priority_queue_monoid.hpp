#pragma once
#include <cplib/collections/dynamic_retroactive_priority_queue.hpp>

namespace cplib {
template <class K, class T, class S>
using DynamicRetroactivePriorityQueueMonoid = DynamicRetroactivePriorityQueue<K, T, S>;

template <class K, class T, class S, class Op, class Lift>
auto initDynamicRetroactivePriorityQueueMonoid(Op op, S e, Lift lift, SortOrder order = Ascending) {
    return DynamicRetroactivePriorityQueueMonoid<K, T, S>(op, e, lift, order);
}
}
