#include <cplib/graph/lowlink.hpp>
#include <cplib/graph/biconnected_components.hpp>
#include <cplib/graph/two_edge_connected_components.hpp>
#include <cplib/graph/block_cut_tree.hpp>
#include <cplib/graph/round_square_tree.hpp>
#include <cplib/collections/unionfind.hpp>
#include <cplib/graph/cycle_detection.hpp>
#include <random>
#include <set>
using namespace cplib;
int main(){
    std::mt19937 rng(27);
    for(int trial=0;trial<1500;++trial){Int n=rng()%10,m=n?rng()%35:0;auto g=initUnWeightedUnDirectedGraph(n);auto s=initWeightedUnDirectedStaticGraph<double>(n);for(Int i=0;i<m;++i){Int a=rng()%n,b=rng()%n;g.add_edge(a,b);s.add_edge(a,b,0.5);}s.build();auto cycle=cycle_detection(g);auto verts=restore_cycle_vertices(g,cycle);std::set<Int> seen(verts.begin(),verts.end());assert(seen.size()==cycle.size());for(std::size_t i=0;i<cycle.size();++i){auto e=g.edge_info[cycle[i]];Int a=verts[i],b=verts[(i+1)%cycle.size()];assert((e.src==a&&e.dst==b)||(e.src==b&&e.dst==a));}assert(cycle==cycle_detection(s));auto ll=initLowLink(g),sl=initLowLink(s);assert(ll.ord==sl.ord&&ll.low==sl.low&&ll.bridges==sl.bridges&&ll.articulation==sl.articulation);auto bc=initBiconnectedComponents(g),bl=initBiconnectedComponents(ll),bs=initBiconnectedComponents(s);assert(bc.groups==bl.groups&&bc.groups==bs.groups&&bc.belong==bl.belong&&bc.articulation==ll.articulation);
        auto component=[&](Int removedVertex,Int removedEdge){std::vector<Int> comp(n,-1);Int count=0;for(Int root=0;root<n;++root)if(root!=removedVertex&&comp[root]==-1){std::vector<Int> stack{root};comp[root]=count++;while(!stack.empty()){Int v=stack.back();stack.pop_back();for(auto e:g.adjacency(v))if(e.id!=removedEdge&&e.dst!=removedVertex&&comp[e.dst]==-1){comp[e.dst]=comp[v];stack.push_back(e.dst);}}}return std::pair{count,comp};};auto [base,comp]=component(-1,-1);assert(cycle.empty()==(m==n-base));std::set<std::pair<Int,Int>> wantBridges;std::vector<bool> bridge(m);for(Int i=0;i<m;++i)if(component(-1,i).first>base){auto e=g.edge_info[i];wantBridges.emplace(std::min(e.src,e.dst),std::max(e.src,e.dst));bridge[i]=true;}std::set<std::pair<Int,Int>> gotBridges;for(auto [a,b]:ll.bridges)gotBridges.emplace(std::min(a,b),std::max(a,b));assert(wantBridges==gotBridges);
        for(Int v=0;v<n;++v){Int delta=component(v,-1).first-base;assert(ll.is_articulation[v]==(delta>0));assert(bc.component_count_delta_after_removal(v)==delta);assert(!bc.belong[v].empty());}
        // 同一ブロックに属する2頂点は、他の任意の1頂点を消しても分離しない。
        for(Int u=0;u<n;++u)for(Int v=u+1;v<n;++v){bool want=comp[u]==comp[v],got=false;for(Int w=0;w<n&&want;++w)if(w!=u&&w!=v){auto c=component(w,-1).second;want=c[u]==c[v];}for(Int a:bc.belong[u])for(Int b:bc.belong[v])got|=a==b;assert(want==got);}
        auto te=initTwoEdgeConnectedComponents(ll);for(Int u=0;u<n;++u)for(Int v=0;v<n;++v){bool want=comp[u]==comp[v];for(Int i=0;i<m&&want;++i){auto c=component(-1,i).second;want=c[u]==c[v];}assert(want==(te.component[u]==te.component[v]));}
        auto block=initBlockCutTree(bc);auto round=initRoundSquareTree(bc);auto checkForest=[&,base=base](const auto& f){auto uf=initUnionFind(f.len);for(auto e:f.edge_info){assert(!uf.issame(e.src,e.dst));uf.unite(e.src,e.dst);}Int roots=0;for(Int v=0;v<f.len;++v)roots+=uf.root(v)==v;assert(roots==base);};checkForest(te.forest);checkForest(block.forest);checkForest(round);for(Int v=0;v<n;++v){if(ll.is_articulation[v])assert(block.id[v]>=Int(bc.groups.size()));else assert(block.id[v]==bc.belong[v][0]);}for(auto e:round.edge_info)assert((e.src<n)!=(e.dst<n));
    }
    // DFSを再帰に置き換えていないことを深いパスで確認する。
    Int n=200000;auto g=initUnWeightedUnDirectedGraph(n);for(Int i=1;i<n;++i)g.add_edge(i-1,i);auto ll=initLowLink(g);assert(Int(ll.bridges.size())==n-1&&Int(ll.articulation.size())==n-2);auto bc=initBiconnectedComponents(g);assert(Int(bc.groups.size())==n-1);assert(cycle_detection(g).empty());g.add_edge(n-1,0);assert(Int(cycle_detection(g).size())==n);
    for(int trial=0;trial<500;++trial){Int n=1+rng()%15;auto d=initUnWeightedDirectedStaticGraph(n);std::vector<std::vector<bool>> reach(n,std::vector<bool>(n));for(int i=0;i<30;++i){Int a=rng()%n,b=rng()%n;d.add_edge(a,b);reach[a][b]=true;}d.build();for(Int k=0;k<n;++k)for(Int i=0;i<n;++i)for(Int j=0;j<n;++j)reach[i][j]=reach[i][j]||(reach[i][k]&&reach[k][j]);bool has=false;for(Int i=0;i<n;++i)has=has||reach[i][i];auto c=cycle_detection(d);assert(has==!c.empty());auto v=restore_cycle_vertices(d,c);std::set<Int> unique(v.begin(),v.end());assert(unique.size()==v.size());for(std::size_t i=0;i<c.size();++i){auto e=d.edge_info[c[i]];assert(e.src==v[i]&&e.dst==v[(i+1)%c.size()]);}}
}
