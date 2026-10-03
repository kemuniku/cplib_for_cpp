#pragma once
#include <cplib/collections/private/dynamic_retroactive_impl.hpp>

namespace cplib {
template <class T>
using IndexedRetroactivePriorityQueue = detail::RetroactiveAVL<Int, T, void, true>;

template <class T> auto initIndexedRetroactivePriorityQueue(SortOrder order = Ascending) {
    return IndexedRetroactivePriorityQueue<T>(order);
}
}
