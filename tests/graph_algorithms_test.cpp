#include <cplib/graph/dijkstra.hpp>
#include <cplib/graph/bellmanford.hpp>
#include <cplib/graph/warshall_floyd.hpp>
#include <cplib/graph/topologicalsort.hpp>
#include <cplib/graph/count_topologicalsort.hpp>
#include <cplib/graph/reverse_edge.hpp>
#include <cplib/graph/bipartite_graph.hpp>
#include <cplib/graph/kruskal.hpp>
#include <cplib/graph/SCC.hpp>
#include <cplib/graph/two_sat.hpp>
#include <cplib/graph/euler_tour.hpp>
#include <random>
using namespace cplib;
int main(){
    std::mt19937 rng(35);
    for(int trial=0;trial<150;++trial){
        Int n=1+rng()%12;auto g=initWeightedDirectedGraph(n);auto ug=initUnWeightedDirectedGraph(n);
        for(Int i=0;i<n;++i)for(Int j=0;j<n;++j)if(rng()%4==0){g.add_edge(i,j,rng()%30);ug.add_edge(i,j);}
        auto all=warshall_floyd(g);for(Int i=0;i<n;++i){assert(dijkstra(g,i)==all[i]);assert(bellmanford(g,i)==all[i]);}
        auto mult=dijkstra(g,std::vector<Int>{0,n-1});for(Int i=0;i<n;++i)assert(mult[i]==std::min(all[0][i],all[n-1][i]));
        auto rev=reverse_edge(g);auto dist=warshall_floyd(rev);for(Int i=0;i<n;++i)for(Int j=0;j<n;++j)assert(dist[i][j]==all[j][i]);
        auto [condensed,map,groups]=SCCG(ug);assert(isDAG(condensed));for(Int i=0;i<n;++i)for(Int j=0;j<n;++j)assert((map[i]==map[j])==(all[i][j]!=INF64&&all[j][i]!=INF64));
        for(auto e:condensed.edge_info)assert(e.src<e.dst);
        for(auto& e:g.edge_info)e.cost=0;
    }
    for(int trial=0;trial<100;++trial){
        Int n=1+rng()%7;auto g=initUnWeightedDirectedGraph(n);for(Int i=0;i<n;++i)for(Int j=i+1;j<n;++j)if(rng()%2)g.add_edge(i,j);
        std::vector<Int> perm(n);std::iota(perm.begin(),perm.end(),0);Int count=0;
        do{std::vector<Int> pos(n);for(Int i=0;i<n;++i)pos[perm[i]]=i;bool ok=true;for(auto e:g.edge_info)if(pos[e.src]>=pos[e.dst])ok=false;count+=ok;}while(std::next_permutation(perm.begin(),perm.end()));assert(count_topologicalsort(g)==count);
    }
    auto neg=initWeightedDirectedGraph(4);neg.add_edge(0,1,1);neg.add_edge(1,2,-2);neg.add_edge(2,1,1);neg.add_edge(2,3,5);assert((bellmanford(neg,0)==std::vector<Int>{0,-INF64,-INF64,-INF64}));
    auto mst=initWeightedUnDirectedStaticGraph(4);mst.add_edge(0,1,4);mst.add_edge(1,2,2);mst.add_edge(0,2,3);mst.add_edge(2,3,1);mst.build();assert(get_MST_cost(mst)==6);assert(to_MST_Graph(mst).edge_count()==3);assert(!is_bipartite_graph(mst));
    for(int trial=0;trial<500;++trial){
        Int n=1+rng()%7;auto p=initTwoSat(n);std::vector<std::tuple<Int,bool,Int,bool>> clauses;
        for(int k=0;k<20;++k){Int i=rng()%n,j=rng()%n;bool f=rng()%2,g=rng()%2;p.add_clause(i,f,j,g);clauses.emplace_back(i,f,j,g);}
        bool exists=false;for(Int mask=0;mask<(Int(1)<<n);++mask){bool ok=true;for(auto [i,f,j,g]:clauses)if(bool(mask>>i&1)!=f&&bool(mask>>j&1)!=g)ok=false;if(ok){exists=true;break;}}
        assert(p.solve()==exists);assert(p.solve()==exists);if(exists){auto answer=p.answer();for(auto [i,f,j,g]:clauses)assert(answer[i]==f||answer[j]==g);}
    }
    auto p=initTwoSat(3);p+=p[0]==true;p+=implies(p[0],p[1]);p+=p[1]^p[2];assert(p.solve());assert(p[0].get()&&p[1].get()&&!p[2].get());auto shared=p;shared+=p[0]==false;assert(!p.solve());
    auto tour=initUnWeightedUnDirectedGraph(3);tour.add_edge(0,0);tour.add_edge(0,1);tour.add_edge(1,2);tour.add_edge(2,0);auto ids=euler_tour(tour,0);assert(ids.size()==4);Int now=0;std::vector<bool> seen(4);for(Int id:ids){assert(!seen[id]);seen[id]=true;auto e=tour.edge_info[id];assert(e.src==now||e.dst==now);now=e.src==now?e.dst:e.src;}assert(now==0);
    auto trail=initUnWeightedDirectedStaticGraph(3);trail.add_edge(0,1);trail.add_edge(1,2);trail.build();assert(euler_tour(trail).empty());assert((euler_trail(trail)==std::vector<Int>{0,1}));assert(euler_trail(trail,1).empty());
}
