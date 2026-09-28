#include <cplib/graph/k_shortest_walk.hpp>
#include <random>
using namespace cplib;
template<class G> std::vector<Int> naive(const G& graph,Int s,Int t,Int k){std::priority_queue<std::pair<Int,Int>,std::vector<std::pair<Int,Int>>,std::greater<>> queue;queue.emplace(0,s);std::vector<Int> count(graph.len),result;while(!queue.empty()&&Int(result.size())<k){auto [cost,u]=queue.top();queue.pop();if(count[u]>=k)continue;++count[u];if(u==t)result.push_back(cost);for(auto [v,c]:graph.to_and_cost(u))queue.emplace(cost+c,v);}result.resize(k,INF64);return result;}
int main(){std::mt19937 rng(81673);for(int trial=0;trial<1800;++trial){Int n=1+rng()%20,m=rng()%100,k=rng()%30,s=rng()%n,t=rng()%n;auto g=initWeightedDirectedGraph<Int>(n);auto sg=initWeightedDirectedStaticGraph<Int>(n);for(Int i=0;i<m;++i){Int u=rng()%n,v=rng()%n,w=rng()%6;g.add_edge(u,v,w);sg.add_edge(u,v,w);}sg.build();auto expected=naive(g,s,t,k);assert(k_shortest_walk(g,s,t,k)==expected);assert(k_shortest_walk(sg,s,t,k)==expected);}
 for(int trial=0;trial<300;++trial){Int n=1+rng()%20;auto g=initUnWeightedUnDirectedGraph(n);for(Int i=0;i<n;++i)g.add_edge(rng()%n,rng()%n);Int s=rng()%n,t=rng()%n;assert(k_shortest_walk(g,s,t,20)==naive(g,s,t,20));}
 auto g=initWeightedDirectedGraph<std::int32_t>(3);g.add_edge(0,1,3);g.add_edge(0,1,3);g.add_edge(1,2,4);assert((k_shortest_walk(g,0,2,4)==std::vector<std::int32_t>{7,7,std::int32_t(INF32),std::int32_t(INF32)}));
 auto chain=initWeightedDirectedGraph<Int>(100000);for(Int i=0;i<99999;++i){chain.add_edge(i,i+1,1);chain.add_edge(i,i+1,2);}auto result=k_shortest_walk(chain,0,99999,100);assert(result[0]==99999);for(Int i=1;i<100;++i)assert(result[i]==100000);
}
