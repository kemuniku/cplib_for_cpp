#pragma once
#include <cplib/collections/private/range_reverse_treap.hpp>
namespace cplib {
template<class T> using RangeReverseArray=detail::RangeReverseTreap<T,detail::NoTreapAction,false,false>;
template<class T> auto initRangeReverseArray(std::span<const T> values){return RangeReverseArray<T>(values,{},T{}, {},{},detail::NoTreapAction{});}
template<class T> auto initRangeReverseArray(const std::vector<T>& values){return RangeReverseArray<T>(values,{},T{}, {},{},detail::NoTreapAction{});}
template<class V> auto toRangeReverseArray(const V& values){return initRangeReverseArray(values);}
}
