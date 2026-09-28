#pragma once
#include <cplib/collections/unionfind.hpp>
namespace cplib {
template<class T=Int> class WeightedUnionFind {
    std::vector<std::int32_t> par_or_siz;
    std::vector<T> potential_diff;
public:
    Int count;
    // ポテンシャル付き集合を構築する。O(n)。
    explicit WeightedUnionFind(Int n):par_or_siz(n,-1),potential_diff(n,T(0)),count(n){}
    Int root(Int x){if(par_or_siz[x]<0)return x;Int p=par_or_siz[x],r=root(p);potential_diff[x]+=potential_diff[p];par_or_siz[x]=r;return r;}
    T potential(Int x){root(x);return potential_diff[x];}
    bool issame(Int x,Int y){return root(x)==root(y);}
    T diff(Int x,Int y){Int rx=root(x),ry=root(y);assert(rx==ry);(void)rx;(void)ry;return potential_diff[y]-potential_diff[x];}
    // potential[y]-potential[x]=wを満たすよう結合する。矛盾ならfalse。償却O(α(n))。
    bool unite(Int x,Int y,T w){w=w+potential(x)-potential(y);x=root(x);y=root(y);if(x==y)return w==T(0);if(par_or_siz[x]>par_or_siz[y]){std::swap(x,y);w=-w;}par_or_siz[x]+=par_or_siz[y];par_or_siz[y]=x;--count;potential_diff[y]=w;return true;}
    Int siz(Int x){return -par_or_siz[root(x)];}
};
template<class T=Int> auto initWeightedUnionFind(Int n){return WeightedUnionFind<T>(n);}
template<class T> T potential(WeightedUnionFind<T>& u,Int x){return u.potential(x);}
template<class T> T diff(WeightedUnionFind<T>& u,Int x,Int y){return u.diff(x,y);}
}
