#include <cplib/graph/maxk_dijkstra.hpp>
#include <cplib/graph/steiner_tree.hpp>
#include <cplib/graph/tsp.hpp>
#include <cplib/collections/unionfind.hpp>
#include <cplib/graph/grid_to_graph.hpp>
#include <cplib/graph/range_edge_graph.hpp>
#include <cplib/graph/merge_tree.hpp>
#include <cplib/graph/dynamic_bipartite.hpp>
#include <cplib/graph/namori_graph.hpp>
#include <cplib/graph/namori_forest.hpp>
#include <cplib/graph/graph_debug.hpp>
#include <random>
using namespace cplib;
void structures_test();
void namori_test();
int main(){
    structures_test();
    namori_test();
    std::mt19937 rng(28);
    for(int trial=0;trial<600;++trial){Int n=1+rng()%40,k=rng()%15;auto g=initWeightedDirectedGraph(n);for(int i=0;i<120;++i)g.add_edge(rng()%n,rng()%n,rng()%(k+1));std::vector<Int> sources(rng()%8);for(auto& v:sources)v=rng()%n;auto want=dijkstra(g,sources);auto a=maxk_dijkstra(g,sources,k);auto b=restore_maxk_dijkstra(g,sources,k);assert(want==a&&want==b.costs);Int start=rng()%n,goal=rng()%n;auto p=shortest_path_maxk_dijkstra(g,start,goal,k);assert(p.cost==dijkstra(g,start)[goal]);for(Int v=0;v<n;++v)if(b.prev[v]>=0){bool good=false;for(auto [to,c]:g.to_and_cost(b.prev[v]))if(to==v&&a[b.prev[v]]+c==a[v])good=true;assert(good);}auto contracted=toContractionGraph(g,sources);auto matrix=to_adjacency_matrix(contracted);for(Int i=0;i<Int(sources.size());++i){auto dist=dijkstra(g,sources[i]);for(Int j=0;j<Int(sources.size());++j)assert(matrix[i][j]==dist[sources[j]]);}}
    {auto g=initWeightedUnDirectedStaticGraph<std::int32_t>(3);g.add_edge(0,1,2);g.add_edge(1,2,1);g.build();auto v=maxk_dijkstra(g,0,2);static_assert(std::is_same_v<decltype(v),std::vector<std::int32_t>>);assert((v==std::vector<std::int32_t>{0,2,3}));}
    {auto g=initUnWeightedUnDirectedGraph(3);g.add_edge(0,1);g.add_edge(1,2);auto v=maxk_dijkstra(g,0,1);static_assert(std::is_same_v<decltype(v),std::vector<Int>>);assert((v==std::vector<Int>{0,1,2}));assert(toContractionGraph(g,std::vector<Int>{0,2}).edge_info[1].cost==2);}
    for(int trial=0;trial<400;++trial){Int n=1+rng()%7,m=rng()%10;auto g=initWeightedUnDirectedGraph(n);for(Int i=0;i<m;++i)g.add_edge(rng()%n,rng()%n,rng()%10);std::vector<Int> terminal(rng()%5);for(auto& v:terminal)v=rng()%n;Int want=terminal.empty()?0:INF64;for(Int bits=0;bits<(Int(1)<<m);++bits){auto uf=initUnionFind(n);Int cost=0;for(Int i=0;i<m;++i)if(bits>>i&1){auto e=g.edge_info[i];uf.unite(e.src,e.dst);cost+=e.cost;}bool good=true;for(Int t:terminal)good&=uf.issame(t,terminal[0]);if(good)want=std::min(want,cost);}assert(steiner_tree_mincost(g,terminal)==want);auto dp=steiner_tree_dp(g,terminal);for(Int v=0;v<n;++v){auto ts=terminal;ts.push_back(v);assert(dp.back()[v]==steiner_tree_mincost(g,ts));}}
    for(int trial=0;trial<250;++trial){Int n=1+rng()%7;std::vector<std::vector<Int>> dist(n,std::vector<Int>(n));for(Int i=0;i<n;++i)for(Int j=0;j<n;++j)dist[i][j]=i==j?0:1+rng()%30;Int start=rng()%n,goal=rng()%n;for(bool floyd:{false,true}){auto d=dist;if(floyd)for(Int k=0;k<n;++k)for(Int i=0;i<n;++i)for(Int j=0;j<n;++j)d[i][j]=std::min(d[i][j],d[i][k]+d[k][j]);std::vector<Int> perm(n);std::iota(perm.begin(),perm.end(),0);Int any=INF64,from=INF64,to=INF64;do{Int sum=0;for(Int i=1;i<n;++i)sum+=d[perm[i-1]][perm[i]];any=std::min(any,sum);if(perm[0]==start){from=std::min(from,sum);to=std::min(to,sum+d[perm.back()][goal]);}}while(std::next_permutation(perm.begin(),perm.end()));assert(tspPathAnyStart(dist,floyd)==any);assert(tspPathCostFrom(dist,start,floyd)==from);assert(tspPathCostFromTo(dist,start,goal,floyd)==to);}}
    {std::vector<std::vector<double>> d={{0,1.5},{2.5,0}};assert(tspPathCostFromTo(d,0,0)==4.0&&tspPathAnyStart(d,0.0,1e100)==1.5);}
}

void structures_test(){
    std::mt19937 rng(29);
    for(int trial=0;trial<150;++trial){Int h=rng()%8,w=rng()%8;std::vector<std::string> grid(h,std::string(w,'.'));for(auto& row:grid)for(auto& c:row)if(rng()%2)c='#';auto d=grid_to_graph(grid);auto s=grid_to_graph<true>(grid);s.build();assert(d.len==h*w&&s.len==d.len&&s.edge_count()==d.edge_count());for(Int i=0;i<h*w;++i)for(Int j=0;j<h*w;++j){bool want=(grid[i/w][i%w]=='.'&&grid[j/w][j%w]=='.'&&(std::abs(i/w-j/w)+std::abs(i%w-j%w)==1)),got=false;for(Int to:d[i])got|=to==j;assert(want==got);}}
    {auto e=initWeightedRangeGraph(0);assert(e.len()==0);e.add_edge(0,0,0,0,1);assert(e.len()==0);}
    for(int trial=0;trial<150;++trial){Int n=1+rng()%30;auto g=initWeightedRangeGraph(n);assert(g.len()==3*n-2&&g.base_len()==3*n-2);auto naive=initWeightedDirectedGraph(n);Int rangeCount=0;for(int q=0;q<70;++q){Int l=rng()%(n+1),r=rng()%(n+1),a=rng()%(n+1),b=rng()%(n+1),cost=rng()%10,u=rng()%n,v=rng()%n;if(l>r)std::swap(l,r);if(a>b)std::swap(a,b);switch(rng()%4){case 0:g.add_edge(l,r,a,b,cost);rangeCount+=l<r&&a<b;for(Int i=l;i<r;++i)for(Int j=a;j<b;++j)naive.add_edge(i,j,cost);break;case 1:g.add_point_to_range_edge(u,l,r,cost);for(Int i=l;i<r;++i)naive.add_edge(u,i,cost);break;case 2:g.add_range_to_point_edge(l,r,v,cost);for(Int i=l;i<r;++i)naive.add_edge(i,v,cost);break;case 3:g.add_edge(u,v,cost);naive.add_edge(u,v,cost);break;}}assert(g.len()==g.base_len()+2*rangeCount);for(Int start=0;start<n;++start){auto dist=dijkstra(g.graph(),start),want=dijkstra(naive,start);dist.resize(n);assert(dist==want);}}
    for(int trial=0;trial<300;++trial){Int n=1+rng()%30;std::vector<std::pair<Int,Int>> queries(70);for(auto& [u,v]:queries){u=rng()%n;v=rng()%n;}auto t=initMergeTree(n,queries);auto uf=initUnionFind(n);std::vector<Int> values(n);for(Int i=0;i<n;++i)values[i]=rng();assert(t.restore_seq(t.make_seq(values))==values);for(auto [u,v]:queries){t.unite(u,v);uf.unite(u,v);for(Int x=0;x<n;++x){auto range=t.get_range(x);assert(range.b-range.a+1==uf.siz(x));for(Int i=range.a;i<=range.b;++i)assert(uf.issame(x,t.et[i]));assert(t.et[t.index(x)]==x);}}}
    for(int trial=0;trial<400;++trial){Int n=1+rng()%8;auto d=initDynamicBipartite(n);std::vector<std::pair<Int,Int>> edges;auto uf=initUnionFind(n);auto best=[&](){Int out=-1;for(Int mask=0;mask<(Int(1)<<n);++mask){bool good=true;for(auto [u,v]:edges)good&=((mask>>u&1)!=(mask>>v&1));if(good)out=std::max(out,Int(std::popcount(UInt(mask))));}return out;};for(int q=0;q<30;++q){Int u=rng()%n,v=rng()%n;bool before=d.is_bipartite();edges.emplace_back(u,v);Int want=best();assert(d.can_unite(u,v)==(want>=0));d.unite(u,v);if(before&&want>=0)uf.unite(u,v);assert(d.is_bipartite()==(want>=0)&&d.cnt_sum==want);if(want>=0)for(Int a=0;a<n;++a)for(Int b=0;b<n;++b)assert(d.issame(a,b)==uf.issame(a,b));}}
}

void namori_test(){
    std::mt19937 rng(30);
    for(int trial=0;trial<400;++trial){Int components=1+rng()%4,n=0;auto g=initUnWeightedUnDirectedGraph(0);std::vector<Int> comp,cycle,cycleRoot;for(Int c=0;c<components;++c){Int count=3+rng()%8,cyc=rng()%2?3+rng()%(count-2):0;g.len+=count;g.edges.resize(g.len);comp.resize(g.len,c);cycle.resize(g.len,false);cycleRoot.resize(g.len);for(Int i=0;i<cyc;++i){g.add_edge(n+i,n+(i+1)%cyc);cycle[n+i]=true;cycleRoot[n+i]=n+i;}for(Int i=std::max(Int(1),cyc);i<count;++i){Int p=rng()%i;g.add_edge(n+i,n+p);cycleRoot[n+i]=cyc?cycleRoot[n+p]:n;}if(!cyc)cycleRoot[n]=n;n+=count;}auto forest=initNamoriForest(g);for(Int u=0;u<n;++u){assert(forest.incycle(u)==bool(cycle[u]));assert(forest.root(u)==cycleRoot[u]);for(Int v=0;v<n;++v){assert(forest.same_component(u,v)==(comp[u]==comp[v]));assert(forest.same_tree(u,v)==(cycleRoot[u]==cycleRoot[v]));std::vector<Int> paths;std::vector<bool> seen(n);auto dfs=[&](auto&& self,Int x,Int d)->void{if(x==v){paths.push_back(d);return;}seen[x]=true;for(Int y:g[x])if(!seen[y])self(self,y,d+1);seen[x]=false;};dfs(dfs,u,0);std::sort(paths.begin(),paths.end());std::pair<Int,Int> want={INF64,INF64};if(!paths.empty())want.first=paths.front();if(paths.size()>1)want.second=paths.back();assert(forest.dist(u,v)==want);}}if(components==1&&std::count(cycle.begin(),cycle.end(),1)){auto single=initNamoriGraph(g);for(Int u=0;u<n;++u)for(Int v=0;v<n;++v){assert(single.dist(u,v)==forest.dist(u,v));assert(single.root(u)==forest.root(u));}}}
    {auto empty=initNamoriForest(initUnWeightedUnDirectedGraph(0));(void)empty;}
    {auto g=initWeightedUnDirectedStaticGraph<>(2);g.add_edge(1,0,3);g.add_edge(1,1,4);g.build();std::ostringstream out;dump_graph(g,1,out);assert(out.str()=="2 3\n1 2 3\n2 2 4\n2 2 4\n");assert(to_graph_graph(g,true).ends_with("data=2+3%0A1+2+3%0A2+2+4%0A2+2+4"));}
    {auto g=initUnWeightedDirectedGraph(2);g.add_edge(1,0);std::ostringstream out;dump_graph(g,out);assert(out.str()=="2 1\n1 0\n");assert(to_graph_graph(g).find("weighted=false&directed=true")!=std::string::npos);}
}
