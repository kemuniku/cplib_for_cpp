#pragma once
#include <cplib/graph/graph.hpp>
#include <cplib/collections/unionfind.hpp>
#include <cplib/utils/backwards_index.hpp>
namespace cplib {
struct MergeTree {
    std::vector<Int> ein,eout,et,ret;UnWeightedUnDirectedGraph tree;UnionFind uf;std::vector<Int> now;Int N,alr_query=0;std::vector<std::pair<Int,Int>> v;
    // 先読みした併合列から葉のEuler順を構築する。O((N+Q)α(N))。
    MergeTree(Int n,std::span<const std::pair<Int,Int>> queries):tree(n+queries.size()+1),uf(n),now(n),N(n),v(queries.begin(),queries.end()){
        auto buildUf=initUnionFind(n);std::vector<Int> buildNow(n);std::iota(buildNow.begin(),buildNow.end(),0);std::iota(now.begin(),now.end(),0);for(Int i=0;i<Int(v.size());++i){auto [a,b]=v[i];Int x=buildNow[buildUf.root(a)],y=buildNow[buildUf.root(b)];tree.add_edge(x,n+i);if(x!=y)tree.add_edge(y,n+i);buildUf.unite(a,b);buildNow[buildUf.root(a)]=n+i;}std::vector<bool> done(n);for(Int i=0;i<n;++i)if(!done[buildUf.root(i)]){done[buildUf.root(i)]=true;tree.add_edge(n+v.size(),buildNow[buildUf.root(i)]);}
        ein.assign(tree.len,-1);eout.assign(tree.len,-1);auto dfs=[&](auto&& self,Int x,Int p)->void{ein[x]=et.size();if(x<n)et.push_back(x);for(Int y:tree[x])if(y!=p)self(self,y,x);eout[x]=et.size();};dfs(dfs,n+v.size(),-1);ret.assign(n,-1);for(Int i=0;i<n;++i)ret[et[i]]=i;
    }
    void unite(Int u,Int w){assert(alr_query<Int(v.size()));assert((v[alr_query].first==u||v[alr_query].first==w)&&(v[alr_query].second==w||v[alr_query].second==u));uf.unite(u,w);now[uf.root(u)]=N+alr_query++;}
    Int get_id(Int x){return now[uf.root(x)];}
    auto get_range(Int x){Int id=get_id(x);return closed_slice(ein[id],eout[id]-1);}
    template<class T> std::vector<T> make_seq(const std::vector<T>& values)const{assert(Int(values.size())==N);std::vector<T> out(N);for(Int i=0;i<N;++i)out[ret[i]]=values[i];return out;}
    template<class T> std::vector<T> restore_seq(const std::vector<T>& values)const{assert(Int(values.size())==N);std::vector<T> out(N);for(Int i=0;i<N;++i)out[et[i]]=values[i];return out;}
    Int index(Int x)const{return ret[x];}
};
inline auto initMergeTree(Int n,std::span<const std::pair<Int,Int>> v){return MergeTree(n,v);}
inline Int get_id(MergeTree& t,Int x){return t.get_id(x);}inline auto get_range(MergeTree& t,Int x){return t.get_range(x);}inline Int index(const MergeTree& t,Int x){return t.index(x);}
template<class T> auto make_seq(const MergeTree& t,const std::vector<T>& v){return t.make_seq(v);}template<class T> auto restore_seq(const MergeTree& t,const std::vector<T>& v){return t.restore_seq(v);}
}
