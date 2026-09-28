#pragma once
#include <cplib/str/static_string.hpp>
#include <cplib/collections/waveletmatrix.hpp>
#include <cplib/utils/itertools.hpp>
#include <unordered_map>
namespace cplib {
namespace detail {struct StaticSearchEntry{virtual ~StaticSearchEntry()=default;};inline std::unordered_map<const void*,std::shared_ptr<StaticSearchEntry>> staticSearchCache;inline const void* staticSearchLastKey=nullptr;inline std::shared_ptr<StaticSearchEntry> staticSearchLast;}
template<class T> class StaticStringSearch {
 struct Cache{Int first=0,last=0,l=0,r=0,length=0,count=0;};
 struct Index{WaveletMatrix wm;std::vector<Cache> cache;Int shift=0;explicit Index(const std::vector<Int>& positions):wm(positions){Int slots=positions.size();while(slots>4096){slots=(slots+1)>>1;++shift;}cache.resize(slots);}};
 std::shared_ptr<Index> index;
 std::pair<Int,Int> suffixRange(const StaticString<T>& pattern)const{Int rank=base->RSA[pattern.l];return {base->RMQ.minLeft(rank,std::int32_t(pattern.len())),base->RMQ.maxRight(rank,std::int32_t(pattern.len()))+1};}
public:
 std::shared_ptr<StaticStringBase<T>> base;
 explicit StaticStringSearch(std::shared_ptr<StaticStringBase<T>> b):base(std::move(b)){std::vector<Int> positions(base->SA.begin(),base->SA.end());index=std::make_shared<Index>(positions);}
 Int count(const StaticString<T>& s,const StaticString<T>& t)const{assert(s.base==base&&t.base==base);Int m=t.len();if(!m)return s.len()+1;if(m>s.len())return 0;Int rank=base->RSA[t.l];auto& c=index->cache[rank>>index->shift];if(c.length==m&&c.l==s.l&&c.r==s.r&&c.first<=rank&&rank<c.last)return c.count;auto [first,last]=suffixRange(t);Int out=index->wm.range_freq(first,last,s.l,Int(s.r)-m+1);c={first,last,s.l,s.r,m,out};return out;}
 bool contains(const StaticString<T>& s,const StaticString<T>& t)const{return count(s,t)>0;}
 // Copy the small search handle into the coroutine so temporary handles are safe.
 static Generator<Int> positions(StaticStringSearch search,StaticString<T> s,StaticString<T> t){assert(s.base==search.base&&t.base==search.base);Int m=t.len();if(!m){for(Int p=0;p<=s.len();++p)co_yield p;}else if(m<=s.len()){auto [first,last]=search.suffixRange(t);auto& wm=search.index->wm;Int begin=wm.range_lowerbound(first,last,s.l),end=wm.range_lowerbound(first,last,Int(s.r)-m+1);for(Int k=begin;k<end;++k)co_yield wm.kth_smallest(first,last,k)-s.l;}}
 auto findAll(StaticString<T> s,StaticString<T> t)const{return positions(*this,s,t);}
 struct View{StaticStringSearch search;StaticString<T> target;Int count(const StaticString<T>& t)const{return search.count(target,t);}bool contains(const StaticString<T>& t)const{return count(t)>0;}auto findAll(StaticString<T> t)const{return search.findAll(target,t);}};
 View operator[](StaticString<T> target)const{assert(target.base==base);return {*this,target};}
};
template<class T> struct StaticStringSearchCacheEntry:detail::StaticSearchEntry{StaticStringSearch<T> search;explicit StaticStringSearchCacheEntry(std::shared_ptr<StaticStringBase<T>> b):search(std::move(b)){}};
inline void clearStaticStringSearchCache(){detail::staticSearchLastKey=nullptr;detail::staticSearchLast.reset();detail::staticSearchCache.clear();}
template<class T> void clearStaticStringSearchCache(const std::shared_ptr<StaticStringBase<T>>& b){if(detail::staticSearchLastKey==b.get()){detail::staticSearchLastKey=nullptr;detail::staticSearchLast.reset();}detail::staticSearchCache.erase(b.get());}
template<class T> StaticStringSearch<T> initStaticStringSearch(const std::shared_ptr<StaticStringBase<T>>& b){auto key=b.get();if(detail::staticSearchLastKey==key&&detail::staticSearchLast)return std::static_pointer_cast<StaticStringSearchCacheEntry<T>>(detail::staticSearchLast)->search;auto it=detail::staticSearchCache.find(key);std::shared_ptr<detail::StaticSearchEntry> entry;if(it==detail::staticSearchCache.end()){entry=std::make_shared<StaticStringSearchCacheEntry<T>>(b);detail::staticSearchCache[key]=entry;}else entry=it->second;detail::staticSearchLastKey=key;detail::staticSearchLast=entry;return std::static_pointer_cast<StaticStringSearchCacheEntry<T>>(entry)->search;}
template<class T> Int count(const StaticString<T>& s,const StaticString<T>& t){assert(s.base==t.base);if(!t.len())return s.len()+1;if(t.len()>s.len())return 0;return initStaticStringSearch(s.base).count(s,t);}
template<class T> bool contains(const StaticString<T>& s,const StaticString<T>& t){return count(s,t)>0;}
template<class T> auto findAll(StaticString<T> s,StaticString<T> t){assert(s.base==t.base);return initStaticStringSearch(s.base).findAll(s,t);}
}
