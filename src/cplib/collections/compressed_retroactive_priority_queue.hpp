#pragma once
#include <cplib/collections/retroactive_priority_queue.hpp>
namespace cplib {
template<class K,class T> class CompressedRetroactivePriorityQueue {
 std::vector<K> coords;RetroactivePriorityQueue<T> queue;
 static std::vector<K> coordinates(std::span<const K> times){std::vector<K> c(times.begin(),times.end());std::sort(c.begin(),c.end());c.erase(std::unique(c.begin(),c.end()),c.end());return c;}
 Int find(const K& t)const{auto it=std::lower_bound(coords.begin(),coords.end(),t);return it==coords.end()||!(*it==t)?-1:it-coords.begin();}
 Int index(const K& t)const{Int i=find(t);assert(i>=0);return i;}
 QueueDelta<K,T> convert(const QueueDelta<Int,T>& d)const{QueueDelta<K,T> out;for(auto e:d.added)out.added.push_back({coords[e.time],e.value});for(auto e:d.removed)out.removed.push_back({coords[e.time],e.value});return out;}
 std::vector<QueueDebugEntry<K,T>> convert(const std::vector<QueueDebugEntry<Int,T>>& entries)const{std::vector<QueueDebugEntry<K,T>> out;for(auto& e:entries){QueueDebugEntry<K,T> c;c.time=coords[e.time];c.kind=e.kind;c.value=e.value;if(e.popped)c.popped=QueueValue<K,T>{coords[e.popped->time],e.popped->value};out.push_back(c);}return out;}
public:
 CompressedRetroactivePriorityQueue(std::span<const K> times,SortOrder order=Ascending):coords(coordinates(times)),queue(coords.size(),order){}
 auto setPush(K t,T v){return convert(queue.setPush(index(t),v));}auto setPop(K t){return convert(queue.setPop(index(t)));}auto erase(K t){return convert(queue.erase(index(t)));}
 Int len()const{return queue.len();}T sum()const{return queue.sum();}T poppedSum()const{return queue.poppedSum();}auto peek()const{return queue.peek();}bool isRemaining(K t)const{Int i=find(t);return i>=0&&queue.isRemaining(i);}
 auto debugOperations()const{return convert(queue.debugOperations());}auto debugTimeline()const{return convert(queue.debugTimeline());}std::string debugDump()const{return formatQueueDebug(debugTimeline());}
};
template<class K,class T> auto initCompressedRetroactivePriorityQueue(std::span<const K> times,SortOrder order=Ascending){return CompressedRetroactivePriorityQueue<K,T>(times,order);}
}
