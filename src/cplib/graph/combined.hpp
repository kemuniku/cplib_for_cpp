#pragma once
#include <cplib/graph/private/combined_hld.hpp>
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
  Int n=next.size();std::vector<Int> removed,stack,indegree(n);auto old_tree=initUnWeightedUnDirectedStaticGraph(n+1);cycle_number.assign(n,-1);cycle_idx.assign(n,-1);roots.resize(n);
  for(Int v:next){assert(0<=v&&v<n);++indegree[v];}
  for(Int v=0;v<n;++v)if(!indegree[v])stack.push_back(v);
  while(!stack.empty()){Int v=stack.back(),u=next[v];stack.pop_back();removed.push_back(v);old_tree.add_edge(v,u);if(!--indegree[u])stack.push_back(u);}
  for(Int v=0;v<n;++v)if(indegree[v]&&cycle_number[v]<0){Int cid=cycle.size(),now=v;cycle.emplace_back();while(cycle_number[now]<0){cycle_number[now]=cid;cycle_idx[now]=cycle.back().size();roots[now]=now;cycle.back().push_back(now);old_tree.add_edge(now,n);now=next[now];}}
  for(auto it=removed.rbegin();it!=removed.rend();++it)roots[*it]=roots[next[*it]];
  old_tree.build();tree=initHld(old_tree,n);depth_tin_.resize(n);
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
inline Int count_kth(const FunctionalGraph& f,Int x,Int k){return f.count_kth(x,k);}
}

#include <cplib/graph/private/combined_segtree.hpp>
namespace cplib {
// HLDの逆順とサイクル順を別々のセグメント木に保持し、非可換な積の順序を保つ。
class FunctionalGraph_with_op:public FunctionalGraph {
 std::function<Int(Int,Int)> op_;Int e_;
 SegmentTree<Int> st_hld_,st_cycle_;
 std::vector<Int> offsets_;
 Int index(Int v)const{return tree.N-1-tree.toSeq(v);}
 Int value(Int v)const{return st_hld_[index(v)];}
public:
 template<class Op> FunctionalGraph_with_op(FunctionalGraph graph,const std::vector<Int>& values,Op op,Int e):FunctionalGraph(std::move(graph)),op_(op),e_(e),st_hld_(0,op,e),st_cycle_(0,op,e){
  assert(values.size()==cycle_number.size());std::vector<Int> hld(tree.N,e),cyc;for(Int v=0;v<Int(values.size());++v)hld[index(v)]=values[v];st_hld_=SegmentTree<Int>(hld,op,e);
  for(const auto& vertices:cycle){offsets_.push_back(cyc.size());for(Int v:vertices)cyc.push_back(values[v]);}st_cycle_=SegmentTree<Int>(cyc,op,e);
 }
 void set(Int x,const Int& value){st_hld_.update(index(x),value);if(incycle(x))st_cycle_.update(offsets_[cycle_number[x]]+cycle_idx[x],value);}
 Int get(Int x)const{return value(x);}
 struct Reference{FunctionalGraph_with_op* owner;Int index;operator Int()const{return owner->get(index);}Reference& operator=(const Int& v){owner->set(index,v);return *this;}Reference& operator=(const Reference& v){return *this=Int(v);}};
 Reference operator[](Int x){return {this,x};}Int operator[](Int x)const{return value(x);}
 // 可換・冪等な積について全始点の到達可能頂点の積をO(N)で求める。
 std::vector<Int> prod_reachable_idempotent_all()const{
  Int n=cycle_number.size();std::vector<Int> result(n,e_);
  for(const auto& vertices:cycle){Int product=e_;for(Int v:vertices)product=op_(product,value(v));for(Int v:vertices)result[v]=product;}
  for(Int v:tree.I)if(v<n&&!incycle(v))result[v]=op_(value(v),result[tree.P[v]]);
  return result;
 }
 // k回の移動で訪れる値の積。木のHLD分解と周回積の二乗法でO(log² N+log k)。
 Int prod(Int start,Int k)const{
  assert(k>=0);
  Int result=e_;Int r=roots[start],last=r;bool ends_in_tree=depth(start)>k;if(ends_in_tree)last=movekth(start,k);
  for(auto [l,u]:tree.path(last,start,ends_in_tree,true))result=op_(result,st_hld_.get(l,u));
  if(ends_in_tree)return result;
  Int cid=cycle_number[r],size=cyclesize(start),offset=offsets_[cid],begin=cycle_idx[r];UInt count=UInt(k-depth(start))+1,full=count/size,remainder=count%size;
  Int product=op_(st_cycle_.get(offset+begin,offset+size),st_cycle_.get(offset,offset+begin));
  for(;full;full>>=1){if(full&1)result=op_(result,product);if(full>1)product=op_(product,product);}
  Int first=std::min<Int>(remainder,size-begin);result=op_(result,st_cycle_.get(offset+begin,offset+begin+first));if(UInt(first)<remainder)result=op_(result,st_cycle_.get(offset,offset+remainder-first));return result;
 }
 // 単調な判定の最後に成功した移動回数を求める。O(log² N+log L)。
 template<class Predicate> Int move_while(Predicate predicate,Int x,Int limit)const{
  assert(0<=limit&&limit<std::numeric_limits<Int>::max());Int count=limit+1,used=1;Int accumulated=value(x);if(!predicate(accumulated))return 0;if(used==count)return limit;
  Int r=roots[x];bool first_segment=true;
  for(auto [l,u]:tree.path(r,x,true,true)){Int begin=l+Int(first_segment);first_segment=false;Int end=begin+std::min(u-begin,count-used);if(begin<end){Int bound=st_hld_.max_right(begin,[&](Int v){return predicate(op_(accumulated,v));});if(bound<end)return used+bound-begin-1;accumulated=op_(accumulated,st_hld_.get(begin,end));used+=end-begin;}if(used==count)return limit;}
  Int cid=cycle_number[r],size=cyclesize(x),offset=offsets_[cid],begin=(cycle_idx[r]+1)%size;
  auto range=[&](Int l,Int u){return st_cycle_.get(offset+l,offset+u);};Int one_cycle=op_(range(begin,size),range(0,begin));Int max_cycles=(count-used)/size;
  if(max_cycles){std::vector<Int> powers{one_cycle};for(Int block=1;block<=max_cycles/2;block*=2)powers.push_back(op_(powers.back(),powers.back()));Int accepted=0;for(Int i=Int(powers.size())-1;i>=0;--i){Int cycles=Int(1)<<i;if(cycles<=max_cycles-accepted){Int candidate=op_(accumulated,powers[i]);if(predicate(candidate)){accumulated=candidate;accepted+=cycles;}}}used+=accepted*size;if(used==count)return limit;}
  auto consume=[&](Int l,Int u){if(l==u)return true;Int bound=st_cycle_.max_right(offset+l,[&](Int v){return predicate(op_(accumulated,v));});if(bound<offset+u){used+=bound-offset-l;return false;}accumulated=op_(accumulated,range(l,u));used+=u-l;return true;};
  Int rest=std::min(count-used,size),first=std::min(rest,size-begin);if(!consume(begin,begin+first))return used-1;rest-=first;if(rest&&!consume(0,rest))return used-1;return used-1;
 }
};
template<class Source,class Op> auto initFunctionalGraph_with_op(const Source& source,const std::vector<Int>& values,Op op,Int e){if constexpr(std::is_same_v<Source,FunctionalGraph>)return FunctionalGraph_with_op(source,values,op,e);else return FunctionalGraph_with_op(initFunctionalGraph(source),values,op,e);}
inline void set(FunctionalGraph_with_op& f,Int x,Int value){f.set(x,value);}
inline auto prod_reachable_idempotent_all(const FunctionalGraph_with_op& f){return f.prod_reachable_idempotent_all();}
inline Int prod(const FunctionalGraph_with_op& f,Int x,Int k){return f.prod(x,k);}
template<class F> Int move_while(const FunctionalGraph_with_op& f,F predicate,Int x,Int limit){return f.move_while(predicate,x,limit);}
}
