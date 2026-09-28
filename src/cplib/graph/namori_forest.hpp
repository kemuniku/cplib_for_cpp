#pragma once
#include <cplib/tree/heavylightdecomposition.hpp>
#include <cplib/utils/constants.hpp>
namespace cplib {
class NamoriForest {
    HeavyLightDecomposition tree;std::vector<Int> rootNo,roots,comp,cyclesize;Int superRoot;
public:
    // 各連結成分が木または単一閉路グラフである森を分解する。O(V+E)。
    explicit NamoriForest(const UnWeightedUnDirectedGraph& graph):superRoot(graph.len){Int n=graph.len;std::vector<Int> stack,sizes(n);rootNo.assign(n,-1);roots.assign(n,-1);comp.assign(n,-1);auto forest=initUnWeightedUnDirectedGraph(n+1);for(Int i=0;i<n;++i){sizes[i]=graph.edges[i].size();if(sizes[i]==1)stack.push_back(i);}while(!stack.empty()){Int i=stack.back();stack.pop_back();for(Int j:graph[i])if(sizes[j]!=1&&--sizes[j]==1)stack.push_back(j);}std::vector<bool> onCycle(n),usedCycle(n);for(Int i=0;i<n;++i)onCycle[i]=sizes[i]>1;for(Int i=0;i<n;++i)for(Int j:graph[i])if(i<j&&!(onCycle[i]&&onCycle[j]))forest.add_edge(i,j);
        for(Int s=0;s<n;++s){if(!onCycle[s]||usedCycle[s])continue;std::vector<Int> cyc;Int now=s,prev=-1;while(true){cyc.push_back(now);usedCycle[now]=true;Int next=-1;for(Int to:graph[now])if(onCycle[to]&&to!=prev){next=to;break;}if(next==-1||next==s)break;prev=now;now=next;}Int cid=cyclesize.size();cyclesize.push_back(cyc.size());for(Int idx=0;idx<Int(cyc.size());++idx){Int r=cyc[idx];rootNo[r]=idx;forest.add_edge(r,superRoot);auto dfs=[&](auto&& self,Int x,Int p)->void{roots[x]=r;comp[x]=cid;for(Int y:graph[x])if(y!=p&&!onCycle[y])self(self,y,x);};dfs(dfs,r,-1);}}
        for(Int s=0;s<n;++s)if(comp[s]==-1){Int cid=cyclesize.size();cyclesize.push_back(0);roots[s]=s;comp[s]=cid;forest.add_edge(s,superRoot);stack={s};while(!stack.empty()){Int x=stack.back();stack.pop_back();for(Int y:graph[x])if(comp[y]==-1){roots[y]=s;comp[y]=cid;stack.push_back(y);}}}
        tree=initHld(forest,superRoot);
    }
    // 別成分なら(INF64,INF64)、経路が1本なら第2要素INF64。O(log V)。
    std::pair<Int,Int> dist(Int u,Int v)const{if(comp[u]!=comp[v])return {INF64,INF64};Int lca=tree.lca(u,v);if(lca!=superRoot)return {tree.dist(u,v),INF64};Int x=roots[u],y=roots[v],cid=comp[u],a=std::abs(rootNo[x]-rootNo[y]),b=cyclesize[cid]-a,base=tree.depth(u)+tree.depth(v)-2;return {std::min(a,b)+base,std::max(a,b)+base};}
    bool incycle(Int x)const{return roots[x]==x&&rootNo[x]!=-1;}Int root(Int x)const{return roots[x];}bool same_tree(Int x,Int y)const{return roots[x]==roots[y];}bool same_component(Int x,Int y)const{return comp[x]==comp[y];}Int component(Int x)const{return comp[x];}
};
inline auto initNamoriForest(const UnWeightedUnDirectedGraph& g){return NamoriForest(g);}
inline auto dist(const NamoriForest& g,Int u,Int v){return g.dist(u,v);}inline bool incycle(const NamoriForest& g,Int x){return g.incycle(x);}inline Int root(const NamoriForest& g,Int x){return g.root(x);}inline bool same_tree(const NamoriForest& g,Int x,Int y){return g.same_tree(x,y);}inline bool same_component(const NamoriForest& g,Int x,Int y){return g.same_component(x,y);}inline Int component(const NamoriForest& g,Int x){return g.component(x);}
}
