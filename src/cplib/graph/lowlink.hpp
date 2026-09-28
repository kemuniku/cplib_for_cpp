#pragma once
#include <cplib/graph/graph.hpp>
namespace cplib {
struct LowLink {
    std::vector<Int> ord,low,parent,preorder,postorder,articulation;
    std::vector<bool> is_articulation;
    std::vector<std::pair<Int,Int>> bridges;
};
// 非再帰DFS。多重辺・自己ループ対応。O(V+E)時間・領域。
template<UnDirectedGraph G> LowLink initLowLink(const G& g){
    if constexpr(StaticGraphTypes<G>)g.static_graph_initialized_check();
    Int n=g.len;LowLink out;out.ord.assign(n,-1);out.low.resize(n);out.parent.assign(n,-1);out.is_articulation.resize(n);std::vector<Int> next(n),children(n),stack;std::vector<bool> skippedParent(n);
    for(Int root=0;root<n;++root){if(out.ord[root]!=-1)continue;out.ord[root]=out.preorder.size();out.low[root]=out.ord[root];out.preorder.push_back(root);stack.push_back(root);
        while(!stack.empty()){Int v=stack.back();auto adj=g.adjacency(v);if(next[v]<Int(adj.size())){Int to=adj[next[v]++].dst;if(to==out.parent[v]&&!skippedParent[v]){skippedParent[v]=true;continue;}if(out.ord[to]==-1){out.parent[to]=v;++children[v];out.ord[to]=out.preorder.size();out.low[to]=out.ord[to];out.preorder.push_back(to);stack.push_back(to);}else out.low[v]=std::min(out.low[v],out.ord[to]);}
            else{stack.pop_back();out.postorder.push_back(v);Int p=out.parent[v];if(p==-1)out.is_articulation[v]=children[v]>1;else{out.low[p]=std::min(out.low[p],out.low[v]);if(out.low[v]>out.ord[p])out.bridges.emplace_back(p,v);if(out.parent[p]!=-1&&out.low[v]>=out.ord[p])out.is_articulation[p]=true;}}
        }
    }for(Int v=0;v<n;++v)if(out.is_articulation[v])out.articulation.push_back(v);return out;
}
}
