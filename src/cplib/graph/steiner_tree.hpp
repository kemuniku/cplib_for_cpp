#pragma once
#include <cplib/graph/dijkstra.hpp>
#include <cplib/utils/bititers.hpp>
namespace cplib {
// 部分集合の併合と多始点Dijkstra。O(V3^K+(V+E)2^K log V)。
template<GraphTypes G,class T=typename G::cost_type> std::vector<std::vector<T>> steiner_tree_dp(const G& g,std::span<const Int> terminal,T zero,T inf){Int k=terminal.size(),n=g.len;std::vector<std::vector<T>> dp(Int(1)<<k,std::vector<T>(n,inf));std::fill(dp[0].begin(),dp[0].end(),zero);for(Int i=0;i<k;++i)dp[Int(1)<<i][terminal[i]]=zero;for(Int bit=1;bit<(Int(1)<<k);++bit){for(Int u=0;u<n;++u)for(Int bn:bitsubset(bit))dp[bit][u]=std::min(dp[bit][u],dp[bn][u]+dp[bit^bn][u]);std::priority_queue<std::pair<T,Int>,std::vector<std::pair<T,Int>>,std::greater<>> q;for(Int u=0;u<n;++u)q.emplace(dp[bit][u],u);while(!q.empty()){auto [d,u]=q.top();q.pop();if(dp[bit][u]!=d)continue;for(auto [v,cost]:g.to_and_cost(u))if(dp[bit][v]>d+cost){dp[bit][v]=d+cost;q.emplace(dp[bit][v],v);}}}return dp;}
template<GraphTypes G,class T=typename G::cost_type> auto steiner_tree_dp(const G& g,std::span<const Int> terminal,T inf=detail::distance_inf<T>()){return steiner_tree_dp(g,terminal,T(0),inf);}
template<GraphTypes G,class T=typename G::cost_type> T steiner_tree_mincost(const G& g,std::span<const Int> terminal,T zero,T inf){if(terminal.empty())return zero;auto dp=steiner_tree_dp(g,terminal,zero,inf);return dp.back()[terminal[0]];}
template<GraphTypes G,class T=typename G::cost_type> T steiner_tree_mincost(const G& g,std::span<const Int> terminal,T inf=detail::distance_inf<T>()){return steiner_tree_mincost(g,terminal,T(0),inf);}
}
