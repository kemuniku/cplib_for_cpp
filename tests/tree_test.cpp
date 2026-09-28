#include <cplib/tree/lca.hpp>
#include <cplib/tree/rerooting.hpp>
#include <cplib/tree/tree_hash.hpp>
#include <cplib/tree/heavylightdecomposition.hpp>
#include <cplib/tree/diameter.hpp>
#include <cplib/tree/cartesiantree.hpp>
#include <cplib/tree/prufer.hpp>
#include <cplib/tree/centroid_decomposition.hpp>
#include <cplib/graph/warshall_floyd.hpp>
#include <random>
using namespace cplib;
int main(){
    std::mt19937 rng(11);
    for(int trial=0;trial<300;++trial){
        Int n=2+rng()%35;std::vector<Int> code(n-2);for(Int& x:code)x=rng()%n;auto g=prufer_decode(code);assert(g.len==n&&g.edge_count()==n-1);
        auto hashes=tree_hash(g);for(Int r=0;r<n;++r){auto h=tree_hash(g,r);assert(h.all_roots==hashes.all_roots);assert(h.subtree[r]==hashes.all_roots[r]);assert(subtree_hash(g,r)==h.subtree);}
        using DP=std::pair<Int,Int>;
        auto sums=solve_Rerooting(g,[](DP a,DP b){return DP{a.first+b.first,a.second+b.second};},DP{0,0},[](DP a,Int,Int){return DP{a.first,a.second+a.first};},[](DP a,Int){return DP{a.first+1,a.second};});
        auto dist=warshall_floyd(g);for(Int v=0;v<n;++v){assert(sums[v].first==n);assert(sums[v].second==std::accumulate(dist[v].begin(),dist[v].end(),Int(0)));}Int want=0;for(auto& row:dist)for(Int d:row){assert(d!=INF64);want=std::max(want,d);}assert(diameter(g)==want);auto [d,path]=diameter_path(g);assert(d==want&&Int(path.size())==d+1);for(std::size_t i=1;i<path.size();++i)assert(dist[path[i-1]][path[i]]==1);
        auto cd=initCentroidDecomposition(g);auto ct=initCentroidDecompositionTree(g);assert(cd.root==ct.root);assert(cd.depth==ct.depth);assert(ct.size[ct.root]==n);assert(ct.tree.edge_count()==n-1);
        auto walk=[&](auto&& self,Int v)->Int{Int size=1;for(Int c:cd.children[v]){assert(cd.parent[c]==v);assert(cd.depth[c]==cd.depth[v]+1);Int sub=self(self,c);assert(sub<=ct.size[v]/2);size+=sub;}assert(size==ct.size[v]);return size;};walk(walk,cd.root);
        std::vector<Int> a(n);for(Int& x:a)x=rng()%30;auto cart=cartesian_tree_tuple(a);Int root=-1;for(Int i=0;i<n;++i)if(cart[i].p==-1){assert(root==-1);root=i;}std::vector<Int> order;auto inorder=[&](auto&& self,Int v)->void{if(v==-1)return;for(Int child:{cart[v].l,cart[v].r})if(child!=-1){assert(cart[child].p==v);assert(a[v]<=a[child]);}self(self,cart[v].l);order.push_back(v);self(self,cart[v].r);};inorder(inorder,root);for(Int i=0;i<n;++i)assert(order[i]==i);
    }
    // 順序付き・番号シャッフル・鎖・スター・長い枝付き経路を親走査と照合。
    for(int mode=0;mode<6;++mode)for(int trial=0;trial<50;++trial){
        Int n=1+rng()%1000;std::vector<Int> p(n,-1),perm(n);std::iota(perm.begin(),perm.end(),0);
        for(Int i=1;i<n;++i)p[i]=mode==0?i-1:mode==1?0:mode>=4?(i<n*3/4?i-1:rng()%i):rng()%i;
        if(mode%2)std::shuffle(perm.begin(),perm.end(),rng);
        std::vector<Int> parent(n,-1);for(Int i=1;i<n;++i)parent[perm[i]]=perm[p[i]];
        Int root=perm[0];auto t=initLCAFromParent(parent,root),no=initLCAFromParent(parent,root,true);
        auto dg=initUnWeightedDirectedGraph(n);std::vector<std::vector<Int>> adj(n);
        for(Int v=0;v<n;++v)if(v!=root){dg.add_edge(v,parent[v]);adj[v].push_back(parent[v]);}
        auto gt=initLCA(dg,root),at=initLCA(adj,root);auto copy=t;auto hp=initHldFromParent(parent,root),hg=initHld(dg,root),ha=initHld(adj,root);
        auto up=[&](Int v,Int k){if(k<0)return Int(-1);while(k--&&v!=-1)v=parent[v];return v;};
        auto dep=[&](Int v){Int d=0;while(parent[v]!=-1){++d;v=parent[v];}return d;};
        auto naive=[&](Int u,Int v){Int a=dep(u),b=dep(v);while(a>b){u=parent[u];--a;}while(b>a){v=parent[v];--b;}while(u!=v){u=parent[u];v=parent[v];}return u;};
        assert(t.numVertices()==n);
        for(Int u=0;u<n;++u){assert(hp.toVtx(hp.toSeq(u))==u);assert(hp.depth(u)==dep(u));Int size=0;for(Int v:hp.subtreeV(u)){assert(naive(u,v)==u);++size;}auto [a,b]=hp.subtree(u);assert(size==b-a);Int cnt=0;for(Int v:hp.children(u)){assert(parent[v]==u);++cnt;}assert(cnt==std::count(parent.begin(),parent.end(),u));Int hc=hp.heavyChildOf(u);if(hc!=-1){assert(parent[hc]==u);for(Int c:hp.children(u))assert(hp.rangeR[c]-hp.rangeL[c]<=hp.rangeR[hc]-hp.rangeL[hc]);}}
        auto reordered=hp.toSeq(parent);for(Int i=0;i<n;++i)assert(reordered[hp.toSeq(i)]==parent[i]);
        for(int q=0;q<2000;++q){Int u=rng()%n,v=rng()%n,w=rng()%n,k=Int(rng()%(n+2))-1,l=naive(u,v);assert(t.parentOf(u)==parent[u]);assert(t.depth(u)==dep(u));assert(t.lca(u,v)==l&&no.lca(u,v)==l&&gt.lca(u,v)==l&&at.lca(u,v)==l&&copy.lca(u,v)==l);assert(t.la(u,k)==up(u,k));Int d=dep(u)+dep(v)-2*dep(l);assert(t.dist(u,v)==d);Int want=k<0||k>d?-1:k<=dep(u)-dep(l)?up(u,k):up(v,d-k);assert(t.la(u,v,k)==want);Int med=naive(u,v)^naive(v,w)^naive(u,w);assert(t.median(u,v,w)==med);assert(hp.lca(u,v)==l&&hg.lca(u,v)==l&&ha.lca(u,v)==l);assert(hp.la(u,v,k)==want);assert(hp.median(u,v,w)==med);
            if(q<10){std::vector<Int> path;for(auto [a,b,rev]:hp.pathWithDirection(u,v))for(Int i=a;i<b;++i)path.push_back(hp.toVtx(rev?n-1-i:i));assert(Int(path.size())==d+1);for(Int i=0;i<=d;++i)assert(path[i]==t.la(u,v,i));std::vector<Int> vertices;for(auto [a,b]:hp.path(u,v))for(Int i=a;i<b;++i)vertices.push_back(hp.toVtx(i));std::sort(vertices.begin(),vertices.end());std::sort(path.begin(),path.end());assert(path==vertices);
                std::vector<Int> ids={u,v,w};auto aux=initAuxiliaryWeightedTree(hp,ids);for(Int x:ids)assert(aux.toi.count(x));assert(aux.graph.edge_count()==Int(aux.v.size())-1);for(auto e:aux.graph.edge_info){assert(e.cost==t.dist(aux.v[e.src],aux.v[e.dst]));assert(naive(aux.v[e.src],aux.v[e.dst])==aux.v[e.src]);}
            }}
    }
    {auto g=initUnWeightedUnDirectedGraph(7);g.add_edge(0,2);g.add_edge(2,4);g.add_edge(1,3);g.add_edge(5,6);auto t=initLCAFromForest(g);assert(t.numVertices()==8);assert(t.lca(4,3)==7);assert(t.parentOf(0)==7&&t.parentOf(1)==7&&t.parentOf(5)==7);assert(t.lca(4,2)==2);assert(t.la(4,6,3)==7);}
    {std::vector<std::vector<Int>> e;auto t=initLCAFromForest(e);assert(t.numVertices()==1&&t.lca(0,0)==0&&t.la(0,1)==-1);}
    {auto g=initUnWeightedUnDirectedGraph(5);g.add_edge(0,1);g.add_edge(0,2);g.add_edge(1,3);g.add_edge(1,4);auto merge=[](std::string a,const std::string& b){return a+b;};auto edge=[](std::string a,Int,Int){return a;};auto vertex=[](std::string a,Int v){return "("+std::to_string(v)+a+")";};auto result=solve_Rerooting(g,merge,std::string{},edge,vertex);assert(result[0]=="(0(1(3)(4))(2))");assert(result[1]=="(1(3)(4)(0(2)))");assert(result[4]=="(4(1(3)(0(2))))");}
    for(int i=0;i<10000;++i){UInt a=((UInt(rng())<<32)|rng())%TREE_HASH_MOD,b=((UInt(rng())<<32)|rng())%TREE_HASH_MOD;assert(detail::treeHashMul(a,b)==UInt((__uint128_t(a)*b)%TREE_HASH_MOD));}
    auto empty=initUnWeightedUnDirectedGraph(0);assert(initCentroidDecomposition(empty).root==-1);assert(initCentroidDecompositionTree(empty).root==-1);
}
