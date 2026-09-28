#pragma once
#include <cplib/graph/maxflow.hpp>
#include <stdexcept>
#include <tuple>
namespace cplib {
namespace detail {
template<class C> C psAdd(C a,C b){if((b>0&&a>std::numeric_limits<C>::max()-b)||(b<0&&a<std::numeric_limits<C>::lowest()-b))throw std::overflow_error("ProjectSelection addition overflow");return a+b;}
template<class C> C psSub(C a,C b){if((b>0&&a<std::numeric_limits<C>::lowest()+b)||(b<0&&a>std::numeric_limits<C>::max()+b))throw std::overflow_error("ProjectSelection subtraction overflow");return a-b;}
}
template<class Cost> struct ProjectSelectionResult{bool feasible=false;Cost min_cost{};std::vector<bool> assignment;};
template<class Cost=Int> class ProjectSelection {
    enum Kind{Unary,Pair,AllGain,Force,Imply};
    struct Term{Kind kind;Int i=0,j=0;bool value=false;std::array<Cost,4> costs{};std::vector<Int> ids;};
    Int n;std::vector<Term> terms;
    void checkIndex(Int i)const{if(i<0||i>=n)throw std::invalid_argument("variable index out of range");}
    static void checkWeight(Cost w){if(w<0)throw std::invalid_argument("weight must be nonnegative");}
public:
    explicit ProjectSelection(Int count):n(count){static_assert(std::is_signed_v<Cost>);if(n<0)throw std::invalid_argument("negative variable count");}
    void add_unary_cost(Int i,Cost c0,Cost c1){checkIndex(i);terms.push_back({Unary,i,0,false,{c0,c1,Cost(0),Cost(0)},{}});}
    void add_cost(Int i,bool value,Cost w){if(value)add_unary_cost(i,Cost(0),w);else add_unary_cost(i,w,Cost(0));}
    void add_gain(Int i,bool value,Cost w){add_cost(i,value,detail::psSub(Cost(0),w));}
    void add_pair_cost(Int i,Int j,Cost c00,Cost c01,Cost c10,Cost c11){checkIndex(i);checkIndex(j);if(i==j){add_unary_cost(i,c00,c11);return;}Cost a=detail::psSub(c01,c00),b=detail::psSub(c11,c10);if(a<b)throw std::invalid_argument("pair cost is not submodular");(void)detail::psSub(a,b);terms.push_back({Pair,i,j,false,{c00,c01,c10,c11},{}});}
    void add_cost_if_true_false(Int i,Int j,Cost w){checkWeight(w);add_pair_cost(i,j,Cost(0),Cost(0),w,Cost(0));}
    void add_cost_if_different(Int i,Int j,Cost w){checkWeight(w);add_pair_cost(i,j,Cost(0),w,w,Cost(0));}
    void add_gain_if_all(std::span<const Int> ids,bool value,Cost w){checkWeight(w);for(Int i:ids)checkIndex(i);if(w==0)return;if(ids.size()==1)add_gain(ids[0],value,w);else if(ids.size()==2){if(value)add_pair_cost(ids[0],ids[1],Cost(0),Cost(0),Cost(0),-w);else add_pair_cost(ids[0],ids[1],-w,Cost(0),Cost(0),Cost(0));}else terms.push_back({AllGain,0,0,value,{w,Cost(0),Cost(0),Cost(0)},std::vector<Int>(ids.begin(),ids.end())});}
    void force(Int i,bool value){checkIndex(i);terms.push_back({Force,i,0,value,{},{}});}
    void imply(Int i,Int j){checkIndex(i);checkIndex(j);terms.push_back({Imply,i,j,false,{},{}});}
    void equal(Int i,Int j){checkIndex(i);checkIndex(j);imply(i,j);imply(j,i);}
    // 補助頂点込みの最小カット。O(V²E)、再実行可。加減算は容量型の溢れを検出する。
    ProjectSelectionResult<Cost> solve()const{
        Int source=n,sink=n+1,vertexCount=n+2;std::vector<std::tuple<Int,Int,Cost>> edges;std::vector<std::pair<Int,Int>> hardEdges;Cost offset=0,total=0;
        auto edge=[&](Int src,Int dst,Cost cap){if(src!=dst&&cap>0){total=detail::psAdd(total,cap);edges.emplace_back(src,dst,cap);}};
        auto unary=[&](Int i,Cost c0,Cost c1){offset=detail::psAdd(offset,std::min(c0,c1));if(c0<=c1)edge(i,sink,detail::psSub(c1,c0));else edge(source,i,detail::psSub(c0,c1));};
        for(const auto& t:terms){auto c=t.costs;switch(t.kind){
            case Unary:unary(t.i,c[0],c[1]);break;
            case Pair:offset=detail::psAdd(offset,c[0]);unary(t.i,Cost(0),detail::psSub(c[3],c[1]));unary(t.j,Cost(0),detail::psSub(c[1],c[0]));edge(t.i,t.j,detail::psSub(detail::psSub(c[1],c[0]),detail::psSub(c[3],c[2])));break;
            case AllGain:if(c[0]==0)break;if(t.ids.empty())offset=detail::psSub(offset,c[0]);else{Int aux=vertexCount++;if(t.value){unary(aux,Cost(0),-c[0]);for(Int i:t.ids)hardEdges.emplace_back(aux,i);}else{unary(aux,-c[0],Cost(0));for(Int i:t.ids)hardEdges.emplace_back(i,aux);}}break;
            case Force:if(t.value)hardEdges.emplace_back(source,t.i);else hardEdges.emplace_back(t.i,sink);break;
            case Imply:if(t.i!=t.j)hardEdges.emplace_back(t.i,t.j);break;
        }}Cost infinity=hardEdges.empty()?Cost(0):detail::psAdd(total,Cost(1));auto graph=initMaxFlow<Cost>(vertexCount);for(auto [src,dst,cap]:edges)graph.add_edge(src,dst,cap);for(auto [src,dst]:hardEdges)graph.add_edge(src,dst,infinity);Cost flow=graph.flow(source,sink,hardEdges.empty()?total:infinity);if(!hardEdges.empty()&&flow==infinity)return {};auto cut=graph.min_cut(source);cut.resize(n);return {true,detail::psAdd(offset,flow),std::move(cut)};
    }
};
template<class Cost=Int> auto initProjectSelection(Int n){return ProjectSelection<Cost>(n);}
template<class Cost> auto solve(const ProjectSelection<Cost>& p){return p.solve();}
}
