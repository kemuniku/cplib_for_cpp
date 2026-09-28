#include <cplib/graph/general_matching.hpp>
#include <fstream>
#include <random>
using namespace cplib;
template<class G> void check(const G& graph,Int expected){auto matching=maximum_matching(graph);assert(Int(matching.size())==expected);std::vector<bool> used(graph.len);for(auto [u,v]:matching){assert(u!=v&&!used[u]&&!used[v]);used[u]=used[v]=true;bool present=false;for(auto [w,c]:graph.to_and_cost(u))if(w==v)present=true;assert(present);}}
Int brute(const std::vector<UInt>& adj){std::vector<int> memo(std::size_t(1)<<adj.size(),-1);auto solve=[&](auto&& self,UInt mask)->int{if(!mask)return 0;auto& ans=memo[mask];if(ans>=0)return ans;Int u=std::countr_zero(mask);UInt rest=mask^(UInt(1)<<u);ans=self(self,rest);for(UInt neighbors=adj[u]&rest;neighbors;neighbors&=neighbors-1){Int v=std::countr_zero(neighbors);ans=std::max(ans,1+self(self,rest^(UInt(1)<<v)));}return ans;};return solve(solve,(UInt(1)<<adj.size())-1);}
int main(int argc,char** argv){std::mt19937 rng(87234);for(int trial=0;trial<2000;++trial){Int n=rng()%17,m=n?rng()%100:0;auto g=initUnWeightedUnDirectedGraph(n);std::vector<UInt> adj(n);for(Int i=0;i<m;++i){Int u=rng()%n,v=rng()%n;g.add_edge(u,v);if(u!=v){adj[u]|=UInt(1)<<v;adj[v]|=UInt(1)<<u;}}check(g,brute(adj));}
 if(argc>1){std::ifstream input(argv[1]);assert(input);Int n,m,expected;while(input>>n>>m>>expected){auto g=initWeightedUnDirectedStaticGraph<Int>(n);for(Int i=0;i<m;++i){Int u,v;input>>u>>v;g.add_edge(u,v,Int(rng()%100)-50);}g.build();check(g,expected);}}
 auto large=initUnWeightedUnDirectedGraph(10000);for(Int i=0;i<10000;i+=2)large.add_edge(i,i+1);for(Int i=0;i<50000;++i)large.add_edge(rng()%10000,rng()%10000);check(large,5000);
 auto isolated=initUnWeightedUnDirectedGraph(200000);isolated.add_edge(191235,181234);isolated.add_edge(5,18);check(isolated,2);
}
