#pragma once
#include <cplib/collections/unionfind.hpp>
namespace cplib {
class PartialPersistentUnionFind {
    std::vector<std::int32_t> par_or_siz;
    std::vector<Int> time;
    std::vector<std::vector<Int>> size_time,size_value;
    Int last=-1;
public:
    // 部分永続集合を構築する。O(n)。
    explicit PartialPersistentUnionFind(Int n):par_or_siz(n,-1),time(n,-1),size_time(n,{-1}),size_value(n,{1}){}
    Int root(Int x,Int t)const{while(time[x]!=-1 && time[x]<=t)x=par_or_siz[x];return x;}
    Int root(Int x)const{return root(x,last);}
    // 非減少の時刻で結合する。O(log n)。
    bool unite(Int u,Int v,Int t){assert(last<=t);last=t;u=root(u,t);v=root(v,t);if(u==v)return false;if(par_or_siz[u]>par_or_siz[v])std::swap(u,v);par_or_siz[u]+=par_or_siz[v];par_or_siz[v]=u;size_time[u].push_back(t);size_value[u].push_back(-par_or_siz[u]);time[v]=t;return true;}
    Int unite(Int u,Int v){unite(u,v,last+1);return last;}
    bool issame(Int u,Int v,Int t)const{return root(u,t)==root(v,t);}
    bool issame(Int u,Int v)const{return issame(u,v,last);}
    Int size(Int x,Int t)const{assert(t>=-1);x=root(x,t);return size_value[x][std::upper_bound(size_time[x].begin(),size_time[x].end(),t)-size_time[x].begin()-1];}
    Int size(Int x)const{return size(x,last);}
    // 初めて連結になった時刻を返す。非連結は-2、自身は-1。O(log n)。
    Int when_unite(Int u,Int v)const{std::vector<Int> tu{u},tv{v};while(par_or_siz[u]>=0){u=par_or_siz[u];tu.push_back(u);}while(par_or_siz[v]>=0){v=par_or_siz[v];tv.push_back(v);}if(u!=v)return -2;while(!tu.empty()&&!tv.empty()&&tu.back()==tv.back()){tu.pop_back();tv.pop_back();}Int r=-1;if(!tu.empty())r=std::max(r,time[tu.back()]);if(!tv.empty())r=std::max(r,time[tv.back()]);return r;}
};
inline PartialPersistentUnionFind initPartialPersistentUnionFind(Int n){return PartialPersistentUnionFind(n);}
inline Int size(const PartialPersistentUnionFind& u,Int x,Int t){return u.size(x,t);}
inline Int size(const PartialPersistentUnionFind& u,Int x){return u.size(x);}
inline Int when_unite(const PartialPersistentUnionFind& u,Int x,Int y){return u.when_unite(x,y);}
}
