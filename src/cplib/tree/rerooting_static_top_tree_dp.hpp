#pragma once
#include <cplib/tree/static_top_tree.hpp>
#include <functional>
namespace cplib {
// 順向き・逆向きクラスタを保持し、構造を変えず根を変えたDPを求める。更新・問い合わせO(log N)。
template<class Forward,class Backward> class RerootingStaticTopTreeDP {
 using Op=std::function<Forward(Forward,Forward)>;using Reverse=std::function<Backward(Backward,Backward)>;using Rake=std::function<Backward(Backward,Forward)>;
 struct State {std::vector<Forward> forward;std::vector<Backward> backward;Op compress,rake;Reverse compressReverse;Rake rakeAtRoot,rakeAtEnd;};std::shared_ptr<State> state_;
 void recalculate(Int node){auto& s=*state_;const auto& x=tree->nodes[node];if(x.kind==sttCompress){s.forward[node]=s.compress(s.forward[x.left],s.forward[x.right]);s.backward[node]=s.compressReverse(s.backward[x.right],s.backward[x.left]);}else if(x.kind==sttRake){s.forward[node]=s.rake(s.forward[x.left],s.forward[x.right]);s.backward[node]=s.rakeAtEnd(s.backward[x.left],s.forward[x.right]);}}
public:
 StaticTopTree tree;
 RerootingStaticTopTreeDP(StaticTopTree t,std::span<const Forward> forward,std::span<const Backward> backward,Op compress,Op rake,Reverse reverse,Rake atRoot,Rake atEnd):state_(std::make_shared<State>(State{std::vector<Forward>(t->nodes.size()),std::vector<Backward>(t->nodes.size()),std::move(compress),std::move(rake),std::move(reverse),std::move(atRoot),std::move(atEnd)})),tree(std::move(t)){assert(Int(forward.size())==tree->numVertices&&Int(backward.size())==tree->numVertices);std::copy(forward.begin(),forward.end(),state_->forward.begin());std::copy(backward.begin(),backward.end(),state_->backward.begin());for(Int node=tree->numVertices;node<Int(tree->nodes.size());++node)recalculate(node);}
 void set(Int v,Forward forward,Backward backward){assert(v>=0&&v<tree->numVertices);state_->forward[v]=std::move(forward);state_->backward[v]=std::move(backward);for(Int node=tree->nodes[v].parent;node!=-1;node=tree->nodes[node].parent)recalculate(node);}
 Forward getAll()const{return state_->forward[tree->root];}
 Backward prod(Int v)const{assert(v>=0&&v<tree->numVertices);const auto& s=*state_;std::vector<Int> path;for(Int node=v;tree->nodes[node].parent!=-1;node=tree->nodes[node].parent)path.push_back(node);Backward upper{};Forward lower{};bool hasUpper=false,hasLower=false;for(auto it=path.rbegin();it!=path.rend();++it){Int child=*it;const auto& x=tree->nodes[tree->nodes[child].parent];if(x.kind==sttCompress){if(child==x.left){lower=hasLower?s.compress(s.forward[x.right],lower):s.forward[x.right];hasLower=true;}else{upper=hasUpper?s.compressReverse(s.backward[x.left],upper):s.backward[x.left];hasUpper=true;}}else if(x.kind==sttRake){assert(hasUpper);if(child==x.left)upper=s.rakeAtRoot(upper,s.forward[x.right]);else{Forward rest=hasLower?s.compress(s.forward[x.left],lower):s.forward[x.left];upper=s.rakeAtRoot(upper,rest);hasLower=false;}}else assert(false);}
  Backward result=s.backward[v];if(hasUpper)result=s.compressReverse(result,upper);if(hasLower)result=s.rakeAtRoot(result,lower);return result;
 }
};
template<class F,class B,class C,class R,class CR,class RR,class RE> auto initRerootingStaticTopTreeDP(StaticTopTree tree,const std::vector<F>& forward,const std::vector<B>& backward,C compress,R rake,CR reverse,RR atRoot,RE atEnd){return RerootingStaticTopTreeDP<F,B>(std::move(tree),forward,backward,std::move(compress),std::move(rake),std::move(reverse),std::move(atRoot),std::move(atEnd));}
template<class F,class B> auto getAll(const RerootingStaticTopTreeDP<F,B>& dp){return dp.getAll();}
template<class F,class B> auto prod(const RerootingStaticTopTreeDP<F,B>& dp,Int v){return dp.prod(v);}
template<class F,class B> void set(RerootingStaticTopTreeDP<F,B>& dp,Int v,F forward,B backward){dp.set(v,std::move(forward),std::move(backward));}
}
