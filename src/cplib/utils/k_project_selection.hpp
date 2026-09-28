#pragma once
#include <cplib/utils/project_selection.hpp>
namespace cplib {
template<class Cost> struct KProjectSelectionResult{bool feasible=false;Cost min_cost{};std::vector<Int> assignment;};
template<class Cost=Int> class KProjectSelection {
    std::vector<Int> sizes,starts;ProjectSelection<Cost> binary{0};
    void checkIndex(Int i)const{if(i<0||i>=Int(sizes.size()))throw std::invalid_argument("variable index out of range");}
    Int ge(Int i,Int lower)const{checkIndex(i);if(lower<0||lower>sizes[i])throw std::invalid_argument("lower threshold out of range");if(lower==0)return 1;if(lower==sizes[i])return 0;return starts[i]+lower-1;}
    Int gt(Int i,Int upper)const{checkIndex(i);if(upper<-1||upper>=sizes[i])throw std::invalid_argument("upper threshold out of range");return ge(i,upper+1);}
public:
    explicit KProjectSelection(std::span<const Int> ranges):sizes(ranges.begin(),ranges.end()),starts(ranges.size()){
        Int count=2;for(Int i=0;i<Int(sizes.size());++i){Int k=sizes[i];if(k<=0)throw std::invalid_argument("domain size must be positive");if(count>std::numeric_limits<Int>::max()-(k-1))throw std::overflow_error("threshold count overflow");starts[i]=count;count+=k-1;}binary=initProjectSelection<Cost>(count);binary.force(0,false);binary.force(1,true);for(Int i=0;i<Int(sizes.size());++i)for(Int t=2;t<sizes[i];++t)binary.imply(starts[i]+t-1,starts[i]+t-2);
    }
    void add_unary_cost(Int i,std::span<const Cost> costs){checkIndex(i);if(Int(costs.size())!=sizes[i])throw std::invalid_argument("unary table length mismatch");std::vector<Cost> diffs(costs.size()-1);for(Int t=1;t<Int(costs.size());++t)diffs[t-1]=detail::psSub(costs[t],costs[t-1]);binary.add_unary_cost(0,costs[0],costs[0]);for(Int t=1;t<Int(costs.size());++t)binary.add_unary_cost(ge(i,t),Cost(0),diffs[t-1]);}
    void add_cost(Int i,Int value,Cost w){checkIndex(i);if(value<0||value>=sizes[i])throw std::invalid_argument("value out of range");std::vector<Cost> costs(sizes[i]);costs[value]=w;add_unary_cost(i,costs);}
    void add_gain(Int i,Int value,Cost w){add_cost(i,value,detail::psSub(Cost(0),w));}
    // Monge費用表を閾値変数へ展開。O(sizes[i]*sizes[j])。
    void add_pair_cost(Int i,Int j,const std::vector<std::vector<Cost>>& costs){checkIndex(i);checkIndex(j);Int ki=sizes[i],kj=sizes[j];if(Int(costs.size())!=ki)throw std::invalid_argument("pair table row mismatch");for(const auto& row:costs)if(Int(row.size())!=kj)throw std::invalid_argument("pair table column mismatch");if(i==j){std::vector<Cost> diagonal(ki);for(Int a=0;a<ki;++a)diagonal[a]=costs[a][a];add_unary_cost(i,diagonal);return;}std::vector<Cost> rowDiffs(ki-1),colDiffs(kj-1);std::vector<std::vector<Cost>> mixed(ki-1,std::vector<Cost>(kj-1));for(Int a=1;a<ki;++a)rowDiffs[a-1]=detail::psSub(costs[a][0],costs[a-1][0]);for(Int b=1;b<kj;++b)colDiffs[b-1]=detail::psSub(costs[0][b],costs[0][b-1]);for(Int a=1;a<ki;++a)for(Int b=1;b<kj;++b){Cost left=detail::psSub(costs[a][b],costs[a][b-1]),right=detail::psSub(costs[a-1][b],costs[a-1][b-1]);if(left>right)throw std::invalid_argument("pair table is not Monge");mixed[a-1][b-1]=detail::psSub(right,left);}binary.add_unary_cost(0,costs[0][0],costs[0][0]);for(Int a=1;a<ki;++a)binary.add_unary_cost(ge(i,a),Cost(0),rowDiffs[a-1]);for(Int b=1;b<kj;++b)binary.add_unary_cost(ge(j,b),Cost(0),colDiffs[b-1]);for(Int a=1;a<ki;++a)for(Int b=1;b<kj;++b){Cost w=mixed[a-1][b-1];if(w!=0)binary.add_pair_cost(ge(i,a),ge(j,b),Cost(0),Cost(0),Cost(0),-w);}}
    void add_cost_if_ge_lt(Int i,Int lower_i,Int j,Int lower_j,Cost w){binary.add_cost_if_true_false(ge(i,lower_i),ge(j,lower_j),w);}
    void add_gain_if_all_ge(std::span<const std::pair<Int,Int>> conditions,Cost w){std::vector<Int> ids;for(auto [i,t]:conditions)ids.push_back(ge(i,t));binary.add_gain_if_all(ids,true,w);}
    void add_gain_if_all_le(std::span<const std::pair<Int,Int>> conditions,Cost w){std::vector<Int> ids;for(auto [i,t]:conditions)ids.push_back(gt(i,t));binary.add_gain_if_all(ids,false,w);}
    void set_min(Int i,Int lower){binary.force(ge(i,lower),true);}
    void set_max(Int i,Int upper){binary.force(gt(i,upper),false);}
    void force(Int i,Int value){checkIndex(i);if(value<0||value>=sizes[i])throw std::invalid_argument("forced value out of range");set_min(i,value);set_max(i,value);}
    void imply(Int i,Int lower_i,Int j,Int lower_j){binary.imply(ge(i,lower_i),ge(j,lower_j));}
    KProjectSelectionResult<Cost> solve()const{auto a=binary.solve();if(!a.feasible)return {};std::vector<Int> assignment(sizes.size());for(Int i=0;i<Int(sizes.size());++i)for(Int t=1;t<sizes[i];++t)if(a.assignment[ge(i,t)])++assignment[i];return {true,a.min_cost,std::move(assignment)};}
};
template<class Cost=Int> auto initKProjectSelection(std::span<const Int> sizes){return KProjectSelection<Cost>(sizes);}
template<class Cost=Int> auto initKProjectSelection(Int n,Int k){if(n<0||k<=0)throw std::invalid_argument("invalid variable count or domain size");return KProjectSelection<Cost>(std::vector<Int>(n,k));}
template<class Cost> auto solve(const KProjectSelection<Cost>& p){return p.solve();}
}
