#pragma once
#include <cplib/tree/static_top_tree.hpp>
#include <functional>
namespace cplib {
// 枝の順序・合法な結合順序に依存しない木DP。値の更新O(log N)、全体取得O(1)。
template<class Forward> class StaticTopTreeDP {
 using Operation=std::function<Forward(Forward,Forward)>;
 struct State {std::vector<Forward> values;Operation compress,rake;};std::shared_ptr<State> state_;
 void recalculate(Int node){const auto& x=tree->nodes[node];if(x.kind==sttCompress)state_->values[node]=state_->compress(state_->values[x.left],state_->values[x.right]);else if(x.kind==sttRake)state_->values[node]=state_->rake(state_->values[x.left],state_->values[x.right]);}
public:
 StaticTopTree tree;
 StaticTopTreeDP(StaticTopTree t,std::span<const Forward> values,Operation compress,Operation rake):state_(std::make_shared<State>(State{std::vector<Forward>(t->nodes.size()),std::move(compress),std::move(rake)})),tree(std::move(t)){assert(Int(values.size())==tree->numVertices);std::copy(values.begin(),values.end(),state_->values.begin());for(Int node=values.size();node<Int(tree->nodes.size());++node)recalculate(node);}
 void set(Int v,Forward value){assert(v>=0&&v<tree->numVertices);state_->values[v]=std::move(value);for(Int node=tree->nodes[v].parent;node!=-1;node=tree->nodes[node].parent)recalculate(node);}
 Forward getAll()const{return state_->values[tree->root];}
};
template<class T,class C,class R> auto initStaticTopTreeDP(StaticTopTree tree,const std::vector<T>& values,C compress,R rake){return StaticTopTreeDP<T>(std::move(tree),values,std::move(compress),std::move(rake));}
template<class T> void set(StaticTopTreeDP<T>& dp,Int v,T value){dp.set(v,std::move(value));}
template<class T> T getAll(const StaticTopTreeDP<T>& dp){return dp.getAll();}
}
