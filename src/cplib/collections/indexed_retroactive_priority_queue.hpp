#pragma once
#include <cplib/collections/private/dynamic_retroactive_impl.hpp>

namespace cplib {
template <class T>
using IndexedRetroactivePriorityQueue = detail::RetroactiveAVL<Int, T, void, true>;

// 空の添字指定型の遡及優先度付きキューを構築する。O(1)。
template <class T> auto initIndexedRetroactivePriorityQueue(SortOrder order = Ascending) {
    return IndexedRetroactivePriorityQueue<T>(order);
}
}
