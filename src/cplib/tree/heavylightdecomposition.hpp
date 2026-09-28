#pragma once
#include <cplib/graph/graph.hpp>
#include <unordered_set>
namespace cplib {
struct HeavyLightDecomposition {
    Int N=0;
    std::vector<Int> P,PP,PD,D,I,rangeL,rangeR;
    // 親より子が後の順序と、元の列挙順を保つ子配列から残りの配列を構築。O(N)。
    void build(const std::vector<std::vector<Int>>& child,Int root){
        std::vector<Int> size(N,1),heavy(N,-1);PP.resize(N);std::iota(PP.begin(),PP.end(),0);
        for(Int i=N-1;i>0;--i){Int v=I[i],p=P[v];size[p]+=size[v];if(heavy[p]==-1||size[heavy[p]]<size[v])heavy[p]=v;}
        for(Int v:I)if(heavy[v]!=-1)PP[heavy[v]]=v;
        PD.assign(N,N);PD[root]=0;D.assign(N,0);
        for(Int v:I)if(v!=root){PP[v]=PP[PP[v]];PD[v]=std::min(PD[PP[v]],PD[P[v]]+1);D[v]=D[P[v]]+1;}
        rangeL.assign(N,0);rangeR.resize(N);
        for(Int p:I){rangeR[p]=rangeL[p]+size[p];Int ir=rangeR[p];for(Int v:child[p])if(v!=heavy[p]){ir-=size[v];rangeL[v]=ir;}if(heavy[p]!=-1)rangeL[heavy[p]]=rangeL[p]+1;}
        for(Int v=0;v<N;++v)I[rangeL[v]]=v;
    }
    Int numVertices() const{return N;}
    Int depth(Int v) const{return D[v];}
    Int toSeq(Int v) const{return rangeL[v];}
    template<class T> std::vector<T> toSeq(const std::vector<T>& values) const{assert(Int(values.size())==N);std::vector<T> result(N);for(Int i=0;i<N;++i)result[toSeq(i)]=values[i];return result;}
    Int toVtx(Int i) const{return I[i];}
    Int toSeq2In(Int v) const{return rangeL[v]*2-D[v];}
    Int toSeq2Out(Int v) const{return rangeR[v]*2-D[v]-1;}
    Int parentOf(Int v) const{return P[v];}
    Int heavyRootOf(Int v) const{return PP[v];}
    Int heavyChildOf(Int v) const{if(toSeq(v)==N-1)return -1;Int c=toVtx(toSeq(v)+1);return PP[v]==PP[c]?c:-1;}
    // Heavy pathを辿ってLCAを取得。O(log N)。
    Int lca(Int u,Int v) const{if(PD[u]<PD[v])std::swap(u,v);while(PD[u]>PD[v])u=P[PP[u]];while(PP[u]!=PP[v]){u=P[PP[u]];v=P[PP[v]];}return D[u]>D[v]?v:u;}
    Int dist(Int u,Int v) const{return D[u]+D[v]-2*D[lca(u,v)];}
    // 祖先から子孫への半開区間。逆向き区間は反転配列の添字。O(log N)。
    std::vector<std::pair<Int,Int>> path(Int r,Int c,bool include_root,bool reverse_path) const{
        Int k=PD[c]-PD[r]+1;if(k<=0)return {};std::vector<std::pair<Int,Int>> result(k);
        for(Int i=0;i<k-1;++i){result[i]={rangeL[PP[c]],rangeL[c]+1};c=P[PP[c]];}
        if(PP[r]!=PP[c]||D[r]>D[c])return {};
        result.back()={rangeL[r]+Int(!include_root),rangeL[c]+1};if(result.back().first==result.back().second)result.pop_back();
        if(reverse_path)for(auto& [l,u]:result){Int old=l;l=N-u;u=N-old;}else std::reverse(result.begin(),result.end());return result;
    }
    std::vector<std::tuple<Int,Int,bool>> pathWithDirection(Int u,Int v) const{Int a=lca(u,v);std::vector<std::tuple<Int,Int,bool>> result;for(auto [l,r]:path(a,u,true,true))result.emplace_back(l,r,true);for(auto [l,r]:path(a,v,false,false))result.emplace_back(l,r,false);return result;}
    std::vector<std::pair<Int,Int>> path(Int u,Int v) const{std::vector<std::pair<Int,Int>> result;for(auto [l,r,up]:pathWithDirection(u,v))result.emplace_back(up?N-r:l,up?N-l:r);return result;}
    std::pair<Int,Int> subtree(Int v) const{return {rangeL[v],rangeR[v]};}
    auto subtreeV(Int v) const{return std::span<const Int>(I).subspan(rangeL[v],rangeR[v]-rangeL[v]);}
    Int median(Int x,Int y,Int z) const{return lca(x,y)^lca(y,z)^lca(x,z);}
    Int la(Int u,Int v,Int d) const{if(d<0)return -1;Int g=lca(u,v),distance=D[u]-2*D[g]+D[v];if(distance<d)return -1;Int p=u;if(D[u]-D[g]<d){p=v;d=distance-d;}while(D[p]-D[PP[p]]<d){d-=D[p]-D[PP[p]]+1;p=P[PP[p]];}return I[rangeL[p]-d];}
    struct ChildrenRange {
        const HeavyLightDecomposition* h;Int first,last;
        struct Iterator{const HeavyLightDecomposition* h;Int index;Int operator*() const{return h->I[index];}Iterator& operator++(){Int v=h->I[index];index+=h->rangeR[v]-h->rangeL[v];return *this;}bool operator!=(const Iterator& b) const{return index!=b.index;}};
        Iterator begin() const{return {h,first};}Iterator end() const{return {h,last};}
    };
    // 子の列挙は追加領域O(1)、時間O(子の数+1)。
    auto children(Int v) const{return ChildrenRange{this,rangeL[v]+1,rangeR[v]};}
};
// 親配列からO(N)で構築する。親の子は元実装と同じ逆番号順に走査。
inline auto initHldFromParent(std::span<const Int> parent,Int root){
    HeavyLightDecomposition h;h.N=parent.size();assert(0<=root&&root<h.N);h.P.assign(parent.begin(),parent.end());h.P[root]=-1;std::vector<std::vector<Int>> child(h.N);
    for(Int v=h.N-1;v>=0;--v)if(v!=root){assert(0<=h.P[v]&&h.P[v]<h.N);child[h.P[v]].push_back(v);}
    h.I.reserve(h.N);h.I.push_back(root);for(std::size_t i=0;i<h.I.size();++i)for(Int v:child[h.I[i]])h.I.push_back(v);assert(Int(h.I.size())==h.N);h.build(child,root);return h;
}
template<UnDirectedGraph G> auto initHld(const G& g,Int root){
    HeavyLightDecomposition h;h.N=g.len;assert(0<=root&&root<h.N);h.P.assign(h.N,-1);h.I.reserve(h.N);h.I.push_back(root);std::vector<std::vector<Int>> child(h.N);
    for(std::size_t i=0;i<h.I.size();++i){Int p=h.I[i];for(auto [v,c]:g.to_and_cost(p))if(v!=h.P[p]){h.I.push_back(v);h.P[v]=p;child[p].push_back(v);}}
    h.build(child,root);return h;
}
namespace detail {
struct HldPairHash{std::size_t operator()(std::pair<Int,Int> p) const{return std::hash<Int>{}(p.first)^(std::hash<Int>{}(p.second)+0x9e3779b97f4a7c15ULL+(UInt(p.first)<<6)+(UInt(p.first)>>2));}};
template<class Adj> auto hld_undirected(const Adj& adj){
    Int n;if constexpr(requires{adj.len;})n=adj.len;else n=adj.size();auto g=initUnWeightedUnDirectedStaticGraph(n);std::unordered_set<std::pair<Int,Int>,HldPairHash> seen;
    for(Int v=0;v<n;++v){auto add=[&](Int u){if(seen.emplace(v,u).second){g.add_edge(v,u);seen.emplace(u,v);}};if constexpr(GraphTypes<Adj>){for(auto [u,c]:adj.to_and_cost(v))add(u);}else for(Int u:adj[v])add(u);}g.build();return g;
}
}
template<DirectedGraph G> auto initHld(const G& g,Int root){return initHld(detail::hld_undirected(g),root);}
inline auto initHld(const std::vector<std::vector<Int>>& adj,Int root){return initHld(detail::hld_undirected(adj),root);}
template<UnDirectedGraph G> auto initHldFromForest(const G& g){Int n=g.len;std::vector<Int> parent(n+1,-1),stack;for(Int root=0;root<n;++root){if(parent[root]!=-1)continue;parent[root]=n;stack.push_back(root);while(!stack.empty()){Int v=stack.back();stack.pop_back();for(auto [u,c]:g.to_and_cost(v))if(parent[u]==-1){parent[u]=v;stack.push_back(u);}}}return initHldFromParent(parent,n);}
template<DirectedGraph G> auto initHldFromForest(const G& g){return initHldFromForest(detail::hld_undirected(g));}
inline auto initHldFromForest(const std::vector<std::vector<Int>>& adj){return initHldFromForest(detail::hld_undirected(adj));}
namespace detail {
inline auto hld_aux_vertices(const HeavyLightDecomposition& h,std::span<const Int> input){std::vector<Int> v(input.begin(),input.end());auto cmp=[&](Int a,Int b){return h.toSeq(a)<h.toSeq(b);};std::sort(v.begin(),v.end(),cmp);Int n=v.size();for(Int i=0;i<n-1;++i)v.push_back(h.lca(v[i],v[i+1]));std::sort(v.begin(),v.end(),cmp);v.erase(std::unique(v.begin(),v.end()),v.end());return v;}
}
// 補助木を期待O(K log K+K log N)で構築。元の頂点番号をラベルに保持。
inline auto initAuxiliaryTree(const HeavyLightDecomposition& h,std::span<const Int> input){auto v=detail::hld_aux_vertices(h,input);auto result=initUnWeightedUnDirectedTableGraph(v);std::vector<Int> stack;for(Int u:v){while(!stack.empty()&&h.toSeq2Out(stack.back())<h.toSeq2In(u))stack.pop_back();if(!stack.empty())result.add_edge(stack.back(),u);stack.push_back(u);}return result;}
template<class S=Int> auto initAuxiliaryWeightedTree(const HeavyLightDecomposition& h,std::span<const Int> input){auto v=detail::hld_aux_vertices(h,input);auto result=initWeightedUnDirectedTableGraph<S>(v);std::vector<Int> stack;for(Int u:v){while(!stack.empty()&&h.toSeq2Out(stack.back())<h.toSeq2In(u))stack.pop_back();if(!stack.empty())result.add_edge(stack.back(),u,S(h.depth(u)-h.depth(stack.back())));stack.push_back(u);}return result;}
inline Int numVertices(const HeavyLightDecomposition& h){return h.numVertices();}
inline Int depth(const HeavyLightDecomposition& h,Int v){return h.depth(v);}
inline Int toSeq(const HeavyLightDecomposition& h,Int v){return h.toSeq(v);}
template<class T> auto toSeq(const HeavyLightDecomposition& h,const std::vector<T>& v){return h.toSeq(v);}
inline Int toVtx(const HeavyLightDecomposition& h,Int v){return h.toVtx(v);}
inline Int toSeq2In(const HeavyLightDecomposition& h,Int v){return h.toSeq2In(v);}
inline Int toSeq2Out(const HeavyLightDecomposition& h,Int v){return h.toSeq2Out(v);}
inline Int parentOf(const HeavyLightDecomposition& h,Int v){return h.parentOf(v);}
inline Int heavyRootOf(const HeavyLightDecomposition& h,Int v){return h.heavyRootOf(v);}
inline Int heavyChildOf(const HeavyLightDecomposition& h,Int v){return h.heavyChildOf(v);}
inline Int lca(const HeavyLightDecomposition& h,Int u,Int v){return h.lca(u,v);}
inline Int dist(const HeavyLightDecomposition& h,Int u,Int v){return h.dist(u,v);}
inline auto path(const HeavyLightDecomposition& h,Int u,Int v){return h.path(u,v);}
inline auto path(const HeavyLightDecomposition& h,Int u,Int v,bool inc,bool rev){return h.path(u,v,inc,rev);}
inline auto pathWithDirection(const HeavyLightDecomposition& h,Int u,Int v){return h.pathWithDirection(u,v);}
inline auto subtree(const HeavyLightDecomposition& h,Int v){return h.subtree(v);}
inline auto subtreeV(const HeavyLightDecomposition& h,Int v){return h.subtreeV(v);}
inline auto children(const HeavyLightDecomposition& h,Int v){return h.children(v);}
inline Int median(const HeavyLightDecomposition& h,Int x,Int y,Int z){return h.median(x,y,z);}
inline Int la(const HeavyLightDecomposition& h,Int u,Int v,Int d){return h.la(u,v,d);}
}
