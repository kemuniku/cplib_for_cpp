#pragma once
#include <cplib/tree/heavylightdecomposition.hpp>
#include <cplib/utils/constants.hpp>
#include <stdexcept>
namespace cplib {
class NamoriGraph {
    HeavyLightDecomposition tree;std::vector<Int> rootNo,roots;Int cyclesize=0;
public:
    // 単一の閉路に木が付いた連結グラフを葉刈りして分解する。O(V+E)。
    explicit NamoriGraph(const UnWeightedUnDirectedGraph& graph){Int n=graph.len;std::vector<Int> stack,sizes(n);rootNo.assign(n,-1);roots.resize(n);auto forest=initUnWeightedUnDirectedGraph(n+1);for(Int i=0;i<n;++i){sizes[i]=graph.edges[i].size();if(sizes[i]==1)stack.push_back(i);}while(!stack.empty()){Int i=stack.back();stack.pop_back();for(Int j:graph[i])if(sizes[j]!=1){forest.add_edge(i,j);if(--sizes[j]==1)stack.push_back(j);}}
        for(Int i=0;i<n;++i)if(sizes[i]!=1){Int now=i,bef=-1,tmp=1;rootNo[i]=0;while(true){forest.add_edge(now,n);Int next=-1;auto dfs=[&](auto&& self,Int x,Int p)->void{roots[x]=now;for(Int y:graph[x])if(p!=y){if(sizes[y]==1)self(self,y,x);else if(bef!=y)next=y;}};dfs(dfs,now,-1);if(next==-1||next==i)break;bef=now;now=next;rootNo[now]=tmp++;}tree=initHld(forest,n);cyclesize=tmp;return;}
        throw std::invalid_argument("NamoriGraphには閉路を持つ連結グラフが必要です");
    }
    // 単純経路の短い距離・長い距離。経路が1本なら第2要素はINF64。O(log V)。
    std::pair<Int,Int> dist(Int u,Int v)const{Int a=tree.lca(u,v);if(a!=Int(rootNo.size()))return {tree.dist(u,v),INF64};Int x=roots[u],y=roots[v],d=std::abs(rootNo[x]-rootNo[y]),other=cyclesize-d,base=tree.depth(u)+tree.depth(v)-2;return {std::min(d,other)+base,std::max(d,other)+base};}
    bool incycle(Int x)const{return roots[x]==x;}Int root(Int x)const{return roots[x];}bool same_tree(Int x,Int y)const{return roots[x]==roots[y];}
};
inline auto initNamoriGraph(const UnWeightedUnDirectedGraph& g){return NamoriGraph(g);}
inline auto dist(const NamoriGraph& g,Int u,Int v){return g.dist(u,v);}inline bool incycle(const NamoriGraph& g,Int x){return g.incycle(x);}inline Int root(const NamoriGraph& g,Int x){return g.root(x);}inline bool same_tree(const NamoriGraph& g,Int x,Int y){return g.same_tree(x,y);}
}
