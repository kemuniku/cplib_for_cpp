#pragma once
#include <cplib/graph/dijkstra.hpp>
namespace cplib {
namespace detail {
// 重み[0,k]用の循環バケット。前駆が不要な版は前駆配列を持たない。O(Vk+E)。
template<bool Restore,GraphTypes G,class Start,class T> auto maxk_dijkstra_core(const G& g,const Start& start,T k,T zero,T inf){
    static_assert(std::is_integral_v<T>);Int bcnt=Int(k+1),cnt=0,pos=0;assert(k>=0);std::vector<std::vector<std::int32_t>> queue(bcnt);T cur=zero;std::vector<T> costs(g.len,inf);std::vector<Int> prev;if constexpr(Restore)prev.assign(g.len,-1);auto init=[&](Int s){queue[0].push_back(std::int32_t(s));costs[s]=zero;};if constexpr(std::is_integral_v<Start>)init(start);else for(Int s:start)init(s);
    while(cnt<bcnt){if(queue[pos].empty()){if(++pos==bcnt)pos=0;++cur;++cnt;continue;}while(!queue[pos].empty()){Int i=queue[pos].back();queue[pos].pop_back();T cost=costs[i];if(cost!=cur)continue;for(auto [j,c]:g.to_and_cost(i)){T temp=cost+c;if(temp<costs[j]){if constexpr(Restore)prev[j]=i;costs[j]=temp;Int pn=pos+Int(c);if(pn>=bcnt)pn-=bcnt;queue[pn].push_back(std::int32_t(j));}}}cnt=0;if(++pos==bcnt)pos=0;++cur;}
    if constexpr(Restore)return RestoredDistances<T>{std::move(costs),std::move(prev)};else return costs;
}
}
template<GraphTypes G,class Start,class T=typename G::cost_type> auto maxk_dijkstra(const G& g,const Start& start,std::type_identity_t<T> k,T zero=T(0),T inf=detail::distance_inf<T>()){return detail::maxk_dijkstra_core<false>(g,start,k,zero,inf);}
template<GraphTypes G,class Start,class T=typename G::cost_type> auto restore_maxk_dijkstra(const G& g,const Start& start,std::type_identity_t<T> k,T zero=T(0),T inf=detail::distance_inf<T>()){return detail::maxk_dijkstra_core<true>(g,start,k,zero,inf);}
template<GraphTypes G,class T=typename G::cost_type> ShortestPath<T> shortest_path_maxk_dijkstra(const G& g,Int start,Int goal,std::type_identity_t<T> k,T zero=T(0),T inf=detail::distance_inf<T>()){auto out=restore_maxk_dijkstra(g,start,k,zero,inf);return {restore_shortest_path_from_prev(out.prev,goal),out.costs[goal]};}
}
