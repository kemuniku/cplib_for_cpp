#pragma once
#include <cplib/graph/maxk_dijkstra.hpp>
namespace cplib {
namespace detail {
template<class T> auto solve_tsp(std::vector<std::vector<T>> dist,const std::vector<T>& start,T inf,bool floyd){Int n=dist.size();if(floyd)for(Int k=0;k<n;++k)for(Int i=0;i<n;++i)for(Int j=0;j<n;++j)dist[i][j]=std::min(dist[i][j],dist[i][k]+dist[k][j]);std::vector<std::vector<T>> dp(n,std::vector<T>(Int(1)<<n,inf));for(Int i=0;i<n;++i)dp[i][Int(1)<<i]=start[i];for(Int bit=0;bit<(Int(1)<<n);++bit)for(Int i=0;i<n;++i)if((bit&(Int(1)<<i))&&dp[i][bit]!=inf)for(Int j=0;j<n;++j)if(!(bit&(Int(1)<<j)))dp[j][bit|(Int(1)<<j)]=std::min(dp[j][bit|(Int(1)<<j)],dp[i][bit]+dist[i][j]);return std::pair{std::move(dp),std::move(dist)};}
}
// 全頂点訪問の部分集合DP。O(N²2^N)、floyd指定時は先に距離閉包を取る。
template<class T> T tspPathCostFrom(const std::vector<std::vector<T>>& dist,Int start_v,T zero,T inf,bool floydwarshall=true){std::vector<T> start(dist.size(),inf);start[start_v]=zero;auto [dp,d]=detail::solve_tsp(dist,start,inf,floydwarshall);T out=inf;for(const auto& row:dp)out=std::min(out,row.back());return out;}
template<class T> T tspPathCostFrom(const std::vector<std::vector<T>>& dist,Int start_v,bool floydwarshall=true){return tspPathCostFrom(dist,start_v,T(0),detail::distance_inf<T>(),floydwarshall);}
// 元実装どおり、最後の訪問点からgoalへ移動する費用も加える（goalの再訪可）。
template<class T> T tspPathCostFromTo(const std::vector<std::vector<T>>& dist,Int start_v,Int goal_v,T zero,T inf,bool floydwarshall=true){std::vector<T> start(dist.size(),inf);start[start_v]=zero;auto [dp,d]=detail::solve_tsp(dist,start,inf,floydwarshall);T out=inf;for(Int i=0;i<Int(dist.size());++i)out=std::min(out,dp[i].back()+d[i][goal_v]);return out;}
template<class T> T tspPathCostFromTo(const std::vector<std::vector<T>>& dist,Int start_v,Int goal_v,bool floydwarshall=true){return tspPathCostFromTo(dist,start_v,goal_v,T(0),detail::distance_inf<T>(),floydwarshall);}
template<class T> T tspPathAnyStart(const std::vector<std::vector<T>>& dist,T zero,T inf,bool floydwarshall=true){std::vector<T> start(dist.size(),zero);auto [dp,d]=detail::solve_tsp(dist,start,inf,floydwarshall);T out=inf;for(const auto& row:dp)out=std::min(out,row.back());return out;}
template<class T> T tspPathAnyStart(const std::vector<std::vector<T>>& dist,bool floydwarshall=true){return tspPathAnyStart(dist,T(0),detail::distance_inf<T>(),floydwarshall);}
// 指定頂点間の最短距離を辺重みにする。重みなし版は循環バケットBFSを使用。
template<WeightedGraph G,class T=typename G::cost_type> auto toContractionGraph(const G& g,std::span<const Int> vertices,T zero=T(0),T inf=detail::distance_inf<T>()){auto out=initWeightedDirectedGraph<T>(vertices.size());for(Int i=0;i<Int(vertices.size());++i){auto res=dijkstra(g,vertices[i],zero,inf);for(Int j=0;j<Int(vertices.size());++j)out.add_edge(i,j,res[vertices[j]]);}return out;}
template<UnWeightedGraph G> auto toContractionGraph(const G& g,std::span<const Int> vertices){auto out=initWeightedDirectedGraph<Int>(vertices.size());for(Int i=0;i<Int(vertices.size());++i){auto res=maxk_dijkstra(g,vertices[i],Int(1));for(Int j=0;j<Int(vertices.size());++j)out.add_edge(i,j,res[vertices[j]]);}return out;}
// 自己辺がなければ対角もnoneのまま。多重辺は最小費用を採用。
template<DirectedGraph G,class T=typename G::cost_type> requires WeightedGraph<G> auto to_adjacency_matrix(const G& g,T none=detail::distance_inf<T>()){std::vector<std::vector<T>> out(g.len,std::vector<T>(g.len,none));for(Int i=0;i<g.len;++i)for(auto [j,c]:g.to_and_cost(i))out[i][j]=std::min(out[i][j],T(c));return out;}
}
