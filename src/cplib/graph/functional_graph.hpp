#pragma once
#include <cplib/tree/heavylightdecomposition.hpp>
namespace cplib {
// 入次数削除でサイクルを抽出し、仮想根に繋いだHLDと逆像の索引をO(N)で構築。
class FunctionalGraph {
 std::vector<std::vector<Int>> depth_tin_;
 std::vector<std::vector<std::vector<Int>>> cycle_depth_;
 std::vector<Int> component_size_;
public:
 HeavyLightDecomposition tree;
 std::vector<Int> cycle_number,cycle_idx,roots;
 std::vector<std::vector<Int>> cycle;
 explicit FunctionalGraph(std::span<const Int> next){
  Int n=next.size();std::vector<Int> removed,indegree(n),parent(n+1,-1);cycle_number.assign(n,-1);cycle_idx.assign(n,-1);roots.resize(n);
  for(Int v:next){assert(0<=v&&v<n);++indegree[v];}
  for(Int v=0;v<n;++v)if(!indegree[v])removed.push_back(v);
  for(std::size_t i=0;i<removed.size();++i){Int v=removed[i],u=next[v];parent[v]=u;if(!--indegree[u])removed.push_back(u);}
  for(Int v=0;v<n;++v)if(indegree[v]&&cycle_number[v]<0){Int cid=cycle.size(),now=v;cycle.emplace_back();while(cycle_number[now]<0){cycle_number[now]=cid;cycle_idx[now]=cycle.back().size();roots[now]=now;cycle.back().push_back(now);parent[now]=n;now=next[now];}}
  for(auto it=removed.rbegin();it!=removed.rend();++it)roots[*it]=roots[next[*it]];
  tree=initHldFromParent(parent,n);depth_tin_.resize(n);
  for(Int tin=0;tin<tree.N;++tin){Int v=tree.toVtx(tin);if(v<n)depth_tin_[tree.depth(v)-1].push_back(tin);}
  cycle_depth_.resize(cycle.size());component_size_.resize(cycle.size());
  for(std::size_t cid=0;cid<cycle.size();++cid)cycle_depth_[cid].resize(cycle[cid].size());
  for(Int d=0;d<n;++d)for(Int tin:depth_tin_[d]){Int v=tree.toVtx(tin),r=roots[v],cid=cycle_number[r],size=cycle[cid].size(),residue=(cycle_idx[r]-d%size+size)%size;cycle_depth_[cid][residue].push_back(d);++component_size_[cid];}
 }
 bool incycle(Int x)const{return cycle_number[x]!=-1;}
 Int depth(Int x)const{return tree.depth(x)-1;}
 Int root(Int x)const{return roots[x];}
 Int cyclesize(Int x)const{return cycle[cycle_number[roots[x]]].size();}
 // HLDによる木部分はO(log N)、サイクル上はO(1)。
 Int movekth(Int x,Int count)const{assert(count>=0);Int d=depth(x);if(d>=count)return tree.la(x,cycle_number.size(),count);Int r=roots[x],cid=cycle_number[r],size=cycle[cid].size();return cycle[cid][(cycle_idx[r]+(count-d)%size)%size];}
 Int canmove_size(Int x)const{return depth(x)+cyclesize(x);}
 Int reachable_to_size(Int x)const{if(incycle(x))return component_size_[cycle_number[x]];auto [l,r]=tree.subtree(x);return r-l;}
 Int dist(Int u,Int v)const{if(cycle_number[roots[u]]!=cycle_number[roots[v]])return -1;Int a=tree.lca(u,v);if(a==v)return depth(u)-depth(v);if(a==Int(cycle_number.size())&&incycle(v)){Int x=cycle_idx[roots[u]],y=cycle_idx[v];return (y-x+cyclesize(v))%cyclesize(v)+depth(u);}return -1;}
 std::vector<Int> get_cycle(Int x)const{return cycle[cycle_number[roots[x]]];}
 // 各サイクルを1頂点に縮約した無向森、元頂点からの対応、根。O(N)。
 auto compressed_forest()const{
  Int n=cycle_number.size(),count=0;std::vector<Int> compressed(n,-1),new_roots(cycle.size());
  for(std::size_t cid=0;cid<cycle.size();++cid){new_roots[cid]=count;for(Int v:cycle[cid])compressed[v]=count;++count;}
  for(Int v=0;v<n;++v)if(compressed[v]<0)compressed[v]=count++;
  auto forest=initUnWeightedUnDirectedGraph(count);for(Int v=0;v<n;++v)if(!incycle(v))forest.add_edge(compressed[v],compressed[tree.P[v]]);
  return std::tuple{std::move(forest),std::move(compressed),std::move(new_roots)};
 }
 Int next(Int x)const{if(!incycle(x))return tree.P[x];Int cid=cycle_number[x];return cycle[cid][(cycle_idx[x]+1)%cycle[cid].size()];}
 std::vector<Int> walk(Int x,Int k)const{assert(0<=x&&x<Int(cycle_number.size())&&0<=k&&k<std::numeric_limits<Int>::max());std::vector<Int> out(k+1);for(Int i=0;i<=k;++i){out[i]=x;if(i<k)x=next(x);}return out;}
 // k回の遷移でxへ到達する始点の数を二分探索で求める。O(log N)。
 Int count_kth(Int x,Int k)const{
  assert(0<=x&&x<Int(cycle_number.size())&&k>=0);
  if(!incycle(x)){Int d=depth(x);if(k>=Int(depth_tin_.size())-d)return 0;auto [l,r]=tree.subtree(x);const auto& tins=depth_tin_[d+k];return std::lower_bound(tins.begin(),tins.end(),r)-std::lower_bound(tins.begin(),tins.end(),l);}
  Int cid=cycle_number[x],size=cycle[cid].size(),residue=(cycle_idx[x]-k%size+size)%size;const auto& values=cycle_depth_[cid][residue];return std::upper_bound(values.begin(),values.end(),k)-values.begin();
 }
};
using Functional_Graph=FunctionalGraph;
inline auto initFunctionalGraph(std::span<const Int> next){return FunctionalGraph(next);}
inline auto initFunctionalGraph(const UnWeightedDirectedGraph& graph){std::vector<Int> next(graph.len);for(Int v=0;v<graph.len;++v)for(auto [u,c]:graph.to_and_cost(v))next[v]=u;return FunctionalGraph(next);}
inline bool incycle(const FunctionalGraph& f,Int x){return f.incycle(x);}
inline Int movekth(const FunctionalGraph& f,Int x,Int k){return f.movekth(x,k);}
inline Int cyclesize(const FunctionalGraph& f,Int x){return f.cyclesize(x);}
inline Int canmove_size(const FunctionalGraph& f,Int x){return f.canmove_size(x);}
inline Int reachable_to_size(const FunctionalGraph& f,Int x){return f.reachable_to_size(x);}
inline Int depth(const FunctionalGraph& f,Int x){return f.depth(x);}
inline Int dist(const FunctionalGraph& f,Int u,Int v){return f.dist(u,v);}
inline auto get_cycle(const FunctionalGraph& f,Int x){return f.get_cycle(x);}
inline Int root(const FunctionalGraph& f,Int x){return f.root(x);}
inline auto compressed_forest(const FunctionalGraph& f){return f.compressed_forest();}
inline auto walk(const FunctionalGraph& f,Int x,Int k){return f.walk(x,k);}
inline Int count_kth(const FunctionalGraph& f,Int x,Int k){return f.count_kth(x,k);}
}
