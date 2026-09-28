#pragma once
#include <cplib/graph/graph.hpp>
#include <random>
namespace cplib {
inline constexpr UInt TREE_HASH_MOD=(UInt(1)<<61)-1;
struct TreeHashResult{std::vector<UInt> subtree,all_roots;};
namespace detail {
struct TreeHashState{UInt hash;Int height;};
inline constexpr TreeHashState treeHashIdentity{1,-1};
inline std::mt19937_64 treeHashRandom{std::random_device{}()};
inline std::vector<UInt> treeHashDepth;
// 元実装と同じ、128bit整数を使わない61bit剰余乗算。
inline UInt treeHashMul(UInt a,UInt b){constexpr UInt mask31=(UInt(1)<<31)-1,mask30=(UInt(1)<<30)-1;UInt au=a>>31,al=a&mask31,bu=b>>31,bl=b&mask31,mid=al*bu+au*bl,value=au*bu*2+(mid>>30)+((mid&mask30)<<31)+al*bl,result=(value>>61)+(value&TREE_HASH_MOD);if(result>=TREE_HASH_MOD)result-=TREE_HASH_MOD;return result;}
inline auto treeHashMerge(TreeHashState a,TreeHashState b){return TreeHashState{treeHashMul(a.hash,b.hash),std::max(a.height,b.height)};}
inline auto treeHashVertex(TreeHashState a){TreeHashState result{a.hash+treeHashDepth[a.height+1],a.height+1};if(result.hash>=TREE_HASH_MOD)result.hash-=TREE_HASH_MOD;return result;}
template<bool Reroot,GraphTypes G> TreeHashResult treeHashImpl(const G& g,Int root){
    TreeHashResult result;Int n=g.len;if(n==0)return result;assert(0<=root&&root<n);
    while(Int(treeHashDepth.size())<n)treeHashDepth.push_back(std::uniform_int_distribution<UInt>(0,TREE_HASH_MOD-1)(treeHashRandom));
    std::vector<Int> parent(n,-2),order{root};parent[root]=-1;std::vector<std::vector<Int>> children(n);
    for(std::size_t i=0;i<order.size();++i){Int u=order[i];for(auto [v,c]:g.to_and_cost(u)){if(v==parent[u])continue;assert(parent[v]==-2);parent[v]=u;children[u].push_back(v);order.push_back(v);}}
    assert(Int(order.size())==n);std::vector<TreeHashState> down(n);result.subtree.resize(n);
    for(Int i=n-1;i>=0;--i){Int u=order[i];auto value=treeHashIdentity;for(Int v:children[u])value=treeHashMerge(value,down[v]);down[u]=treeHashVertex(value);result.subtree[u]=down[u].hash;}
    if constexpr(Reroot){std::vector<TreeHashState> up(n),prefix;up[root]=treeHashIdentity;result.all_roots.resize(n);
        for(Int u:order){Int count=children[u].size();prefix.resize(count+1);prefix[0]=up[u];for(Int i=0;i<count;++i)prefix[i+1]=treeHashMerge(prefix[i],down[children[u][i]]);result.all_roots[u]=treeHashVertex(prefix[count]).hash;auto suffix=treeHashIdentity;for(Int i=count-1;i>=0;--i){Int v=children[u][i];up[v]=treeHashVertex(treeHashMerge(prefix[i],suffix));suffix=treeHashMerge(down[v],suffix);}}
    }return result;
}
}
// 高さごとの乱数を共用。構築・全頂点のハッシュ取得とも時間・空間O(N)。
template<GraphTypes G> auto subtree_hash(const G& g,Int root=0){return detail::treeHashImpl<false>(g,root).subtree;}
template<UnDirectedGraph G> auto tree_hash(const G& g,Int root=0){return detail::treeHashImpl<true>(g,root);}
template<UnDirectedGraph G> auto all_roots_hash(const G& g,Int root=0){return tree_hash(g,root).all_roots;}
}
