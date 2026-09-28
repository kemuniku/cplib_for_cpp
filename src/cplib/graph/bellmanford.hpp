#pragma once
#include <cplib/graph/dijkstra.hpp>
namespace cplib {
// Bellman–Ford法で最短距離を求める。負閉路の影響先は-inf。O(VE)。
template<GraphTypes G,class Start,class T=typename G::cost_type> RestoredDistances<T> restore_bellmanford(const G& g,const Start& start,T zero=T(0),T inf=detail::distance_inf<T>()){
    RestoredDistances<T> result{std::vector<T>(g.len,inf),std::vector<Int>(g.len,-1)};auto& costs=result.costs;auto& prev=result.prev;
    if constexpr(std::is_integral_v<Start>)costs[start]=zero;else for(Int s:start)costs[s]=zero;
    bool changed=false;
    for(Int step=0;step<g.len;++step){changed=false;for(Int i=0;i<g.len;++i){if(costs[i]==inf)continue;for(auto [j,c]:g.to_and_cost(i)){T temp=costs[i]+c;if(temp<costs[j]){prev[j]=i;costs[j]=temp;changed=true;}}}if(!changed)break;}
    if(changed)for(Int step=0;step<g.len;++step)for(Int i=0;i<g.len;++i){if(costs[i]==inf)continue;for(auto [j,c]:g.to_and_cost(i)){T temp=costs[i]+c;if(temp<costs[j]){costs[j]=-inf;prev[j]=-1;}}}
    return result;
}
template<GraphTypes G,class Start,class T=typename G::cost_type> auto bellmanford(const G& g,const Start& start,T zero=T(0),T inf=detail::distance_inf<T>()){return restore_bellmanford(g,start,zero,inf).costs;}
template<GraphTypes G,class T=typename G::cost_type> ShortestPath<T> shortest_path_bellmanford(const G& g,Int start,Int goal,T zero=T(0),T inf=detail::distance_inf<T>()){auto result=restore_bellmanford(g,start,zero,inf);return {restore_shortest_path_from_prev(result.prev,goal),result.costs[goal]};}
}
