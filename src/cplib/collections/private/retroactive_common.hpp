#pragma once
#include <cplib/common.hpp>
#include <cplib/collections/private/sort_order.hpp>
#include <optional>
#include <queue>
#include <sstream>
#include <type_traits>
namespace cplib {
enum QueueDebugKind{qdkNone,qdkPush,qdkPop};
template<class K,class T> struct QueueValue{K time;T value;};
template<class K,class T> struct QueueDebugEntry{K time{};QueueDebugKind kind=qdkNone;std::optional<T> value;std::optional<QueueValue<K,T>> popped;};
template<class K,class T> struct QueueDelta{std::vector<QueueValue<K,T>> added,removed;};
template<class K,class T> void replayQueueDebug(std::vector<QueueDebugEntry<K,T>>& entries,SortOrder order){auto cmp=[&](Int a,Int b){auto& x=*entries[a].value;auto& y=*entries[b].value;if(x<y)return order==Descending;if(y<x)return order==Ascending;return a>b;};std::priority_queue<Int,std::vector<Int>,decltype(cmp)> heap(cmp);for(Int i=0;i<Int(entries.size());++i){auto& e=entries[i];e.popped.reset();if(e.kind==qdkPush)heap.push(i);else if(e.kind==qdkPop&&!heap.empty()){Int j=heap.top();heap.pop();e.popped=QueueValue<K,T>{entries[j].time,*entries[j].value};}}}
template<class R> std::string formatQueueDebug(const R& entries){std::ostringstream out;bool first=true;for(auto& e:entries){if(!first)out<<'\n';first=false;out<<e.time<<": ";if(e.kind==qdkNone)out<<"noop";else if(e.kind==qdkPush)out<<"push("<<*e.value<<')';else if(!e.popped)out<<"pop -> empty";else out<<"pop -> "<<e.popped->value<<" (push at "<<e.popped->time<<')';}return out.str();}
namespace detail {
template<class T> T queue_add(T a,T b){if constexpr(std::is_integral_v<T>&&std::is_signed_v<T>){using U=std::make_unsigned_t<T>;return std::bit_cast<T>(U(U(a)+U(b)));}else return a+b;}
template<class T> T queue_sub(T a,T b){if constexpr(std::is_integral_v<T>&&std::is_signed_v<T>){using U=std::make_unsigned_t<T>;return std::bit_cast<T>(U(U(a)-U(b)));}else return a-b;}
}
}
