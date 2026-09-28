#pragma once
#include <cplib/tree/private/link_cut_tree_base.hpp>
namespace cplib::detail {
template<class S,class F> struct LazySubtreeLinkCutTreeState:LinkCutTreeOperations<LazySubtreeLinkCutTreeState<S,F>,S> {
 struct Node {Int left=0,right=0,parent=0;bool rev=false;S value,prod,rprod,virtualValue,all;F lazy,cancel,pathLazy;bool hasPathLazy=false;};
 std::vector<Node> nodes;std::vector<Int> stack;std::function<S(S,S)> merge;S defaultValue;std::function<S(F,S)> mapping;std::function<F(F,F)> composition;F id;std::function<S(S)> inverse;std::function<F(F)> inverseAction;
 LazySubtreeLinkCutTreeState(std::span<const S> values,std::function<S(S,S)> op,S e,std::function<S(F,S)> map,std::function<F(F,F)> comp,F identity,std::function<S(S)> inv,std::function<F(F)> invAction):merge(std::move(op)),defaultValue(e),mapping(std::move(map)),composition(std::move(comp)),id(identity),inverse(std::move(inv)),inverseAction(std::move(invAction)){assert(inverse&&inverseAction);nodes.reserve(values.size()+1);nodes.push_back({0,0,0,false,e,e,e,e,e,id,id,id,false});for(const auto& x:values)nodes.push_back({0,0,0,false,x,x,x,e,x,id,id,id,false});}
 // 補助木とvirtual child全体へ作用し、累積タグを更新する。
 void applyAll(Int v,F f){if(!v)return;auto& x=nodes[v];x.value=mapping(f,x.value);x.prod=mapping(f,x.prod);x.rprod=mapping(f,x.rprod);x.virtualValue=mapping(f,x.virtualValue);x.all=mapping(f,x.all);x.lazy=composition(f,x.lazy);}
 // パスだけへの作用は、全体集約に差分として反映する。
 void applyPathNode(Int v,F f){if(!v)return;auto& x=nodes[v];S old=x.prod;x.value=mapping(f,x.value);x.prod=mapping(f,old);x.rprod=mapping(f,x.rprod);x.all=merge(x.all,merge(inverse(old),x.prod));if(x.hasPathLazy)x.pathLazy=composition(f,x.pathLazy);else{x.pathLazy=f;x.hasPathLazy=true;}}
 void push(Int v){if(!v)return;Int p=nodes[v].parent;if(p){bool pending=true;if constexpr(requires(F a,F b){a==b;})pending=!(nodes[p].lazy==nodes[v].cancel);if(pending){applyAll(v,composition(nodes[p].lazy,inverseAction(nodes[v].cancel)));nodes[v].cancel=nodes[p].lazy;}}if(nodes[v].rev){for(Int child:{nodes[v].left,nodes[v].right})this->toggle(child);nodes[v].rev=false;}if(nodes[v].hasPathLazy){applyPathNode(nodes[v].left,nodes[v].pathLazy);applyPathNode(nodes[v].right,nodes[v].pathLazy);nodes[v].pathLazy=id;nodes[v].hasPathLazy=false;}}
 // 旧親の未受領作用を反映してから付け替え、新親の過去の作用を打ち消す。
 void setParent(Int v,Int parent){if(!v)return;push(v);nodes[v].parent=parent;nodes[v].cancel=nodes[parent].lazy;}
 void pull(Int v){Int l=nodes[v].left,r=nodes[v].right;push(l);push(r);nodes[v].prod=merge(merge(nodes[l].prod,nodes[v].value),nodes[r].prod);nodes[v].rprod=merge(merge(nodes[r].rprod,nodes[v].value),nodes[l].rprod);nodes[v].all=merge(merge(nodes[l].all,nodes[r].all),merge(nodes[v].value,nodes[v].virtualValue));}
 void addVirtual(Int v,Int child){if(!child)return;push(child);nodes[v].virtualValue=merge(nodes[v].virtualValue,nodes[child].all);}
 void removeVirtual(Int v,Int child){if(!child)return;push(child);nodes[v].virtualValue=merge(nodes[v].virtualValue,inverse(nodes[child].all));}
 void pathApply(Int u,Int v,F f){this->check(u);this->check(v);this->makeRoot(u);this->accessNode(v+1);applyPathNode(v+1,f);}
 void componentApply(Int v,F f){this->check(v);this->accessNode(v+1);applyAll(v+1,f);}
 void subtreeApply(Int v,Int parent,F f){this->cut(v,parent);componentApply(v,f);this->link(v,parent);}
};
}
namespace cplib {
// 集約・作用は可換群、mappingは群作用であること。パス・部分木・成分の更新は償却O(log N)。
template<class S,class F> class LazySubtreeLinkCutTree:public detail::LinkCutTreeHandle<detail::LazySubtreeLinkCutTreeState<S,F>,S> {
 using State=detail::LazySubtreeLinkCutTreeState<S,F>;using Base=detail::LinkCutTreeHandle<State,S>;
public:
 LazySubtreeLinkCutTree(std::span<const S> values,std::function<S(S,S)> merge,S defaultValue,std::function<S(F,S)> mapping,std::function<F(F,F)> composition,F id,std::function<S(S)> inverse,std::function<F(F)> inverseAction):Base(std::make_shared<State>(values,std::move(merge),std::move(defaultValue),std::move(mapping),std::move(composition),std::move(id),std::move(inverse),std::move(inverseAction))){}
 LazySubtreeLinkCutTree(Int n,std::function<S(S,S)> merge,S defaultValue,std::function<S(F,S)> mapping,std::function<F(F,F)> composition,F id,std::function<S(S)> inverse,std::function<F(F)> inverseAction):LazySubtreeLinkCutTree(std::vector<S>(checked(n),defaultValue),std::move(merge),defaultValue,std::move(mapping),std::move(composition),std::move(id),std::move(inverse),std::move(inverseAction)){}
 void pathApply(Int u,Int v,F f)const{this->impl_->pathApply(u,v,std::move(f));}void componentApply(Int v,F f)const{this->impl_->componentApply(v,std::move(f));}void subtreeApply(Int v,Int parent,F f)const{this->impl_->subtreeApply(v,parent,std::move(f));}
private:static Int checked(Int n){assert(n>=0);return n;}
};
template<class V,class Op,class S,class Map,class Comp,class F,class Inv,class InvAction> auto initLazySubtreeLinkCutTree(V&& values,Op merge,S defaultValue,Map mapping,Comp composition,F id,Inv inverse,InvAction inverseAction){return LazySubtreeLinkCutTree<S,F>(std::forward<V>(values),std::move(merge),std::move(defaultValue),std::move(mapping),std::move(composition),std::move(id),std::move(inverse),std::move(inverseAction));}
template<class... A> auto newLazySubtreeLinkCutTreeWith(A&&... args){return initLazySubtreeLinkCutTree(std::forward<A>(args)...);}
}
