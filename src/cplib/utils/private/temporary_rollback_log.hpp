#pragma once
#include <cplib/utils/private/auto_rollback.hpp>
#include <unordered_map>
namespace cplib {
struct TemporaryCheckpoint{Int position,parentStart;};
struct TemporaryIndexCache{std::vector<Int> stamps;const void* owner=nullptr;const void* base=nullptr;std::size_t size=0;Int length=0;Int indexCapacity()const{return stamps.size();}};
class TemporaryRollbackLog {
 struct Location{void* address;std::size_t size;bool operator==(const Location&)const=default;};
 struct Hash{std::size_t operator()(Location x)const{return std::hash<void*>{}(x.address)^(x.size*0x9e3779b9);}};
 struct Entry{Location location;Int previous;std::size_t offset;Int* stamp;};
 std::vector<Entry> entries;std::unordered_map<Location,Int,Hash> latest;Int scopeStart=0;std::vector<unsigned char> saved;std::vector<TemporaryIndexCache*> caches;
 template<class T> void save(T* p,Int previous,Int* stamp=nullptr){static_assert(std::is_trivially_copyable_v<T>);auto offset=saved.size();saved.resize(offset+sizeof(T));std::memcpy(saved.data()+offset,p,sizeof(T));entries.push_back({{p,sizeof(T)},previous,offset,stamp});}
public:
 TemporaryRollbackLog()=default;TemporaryRollbackLog(const TemporaryRollbackLog&)=delete;TemporaryRollbackLog& operator=(const TemporaryRollbackLog&)=delete;~TemporaryRollbackLog(){restore(0);}
 Int len()const{return entries.size();}
 template<class T> void remember(T* p){Location key{p,sizeof(T)};auto it=latest.find(key);Int prev=it==latest.end()?-1:it->second;if(prev>=scopeStart)return;Int i=len();save(p,prev);latest[key]=i;}
 template<class T> void rememberIndexed(T* p,T* base,Int length,TemporaryIndexCache& cache){if(!cache.owner){if(length>Int(cache.stamps.size()))cache.stamps.resize(length);cache.owner=this;cache.base=base;cache.size=sizeof(T);cache.length=length;caches.push_back(&cache);}auto offset=reinterpret_cast<std::uintptr_t>(p)-reinterpret_cast<std::uintptr_t>(base);if(cache.owner!=this||cache.base!=base||cache.size!=sizeof(T)||cache.length!=length||offset/sizeof(T)>=UInt(length)||offset%sizeof(T)){remember(p);return;}auto& stamp=cache.stamps[offset/sizeof(T)];Int prev=stamp-1;if(prev>=scopeStart)return;Int i=len();save(p,prev,&stamp);stamp=i+1;}
 template<class T> T& write(T& x){remember(&x);return x;}
 template<class T,class V> void set(T& x,V&& v){write(x)=std::forward<V>(v);}
 void restore(Int position){assert(0<=position&&position<=len());while(len()>position){auto e=entries.back();entries.pop_back();std::memcpy(e.location.address,saved.data()+e.offset,e.location.size);saved.resize(e.offset);if(e.stamp)*e.stamp=e.previous+1;else if(e.previous<0)latest.erase(e.location);else latest[e.location]=e.previous;}if(entries.empty()){for(auto p:caches)p->owner=nullptr;caches.clear();}}
 TemporaryCheckpoint beginTemporary(){TemporaryCheckpoint c{len(),scopeStart};scopeStart=len();return c;}
 void endTemporary(TemporaryCheckpoint c){restore(c.position);scopeStart=c.parentStart;}
};
}
