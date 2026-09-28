#include <cplib/graph/functional_graph_with_op.hpp>
#include <cplib/graph/functional_graph_with_lazy_op.hpp>
#include <random>
using namespace cplib;
struct Sum {Int sum=0,size=0;bool operator==(const Sum&)const=default;};
struct Affine {Int a=1,b=0;};
Sum sum_op(Sum a,Sum b){return {a.sum+b.sum,a.size+b.size};}
Sum sum_e(){return {};}
Sum sum_mapping(Affine f,Sum x){return {f.a*x.sum+f.b*x.size,x.size};}
Affine sum_composition(Affine f,Affine g){return {f.a*g.a,f.a*g.b+f.b};}
Affine sum_id(){return {};}
struct ExplicitTree:LazySegmentTree<Sum,Affine>{explicit ExplicitTree(const std::vector<Sum>& v):LazySegmentTree(v,sum_op,sum_e(),sum_mapping,sum_composition,sum_id()){}};
void test_lazy(){
 std::mt19937 rng(8973);
 for(int trial=0;trial<150;++trial){Int n=1+rng()%80;std::vector<Int> next(n),plain(n);std::vector<Sum> values(n);for(auto& x:next)x=rng()%n;for(Int i=0;i<n;++i){plain[i]=rng()%10;values[i]={plain[i],1};}
  auto graph=initFunctionalGraph_with_lazy_op(next,values,sum_op,sum_e,sum_mapping,sum_composition,sum_id);
  for(int q=0;q<120;++q){Int start=rng()%n,k=rng()%300;bool include=rng()%2;
   if(q%3==0){Affine f{Int(rng()%2),Int(rng()%5)};graph.apply(start,k,f,include);Int now=start;for(Int j=0;j<=k;++j){if(include||j)plain[now]=f.a*plain[now]+f.b;now=next[now];}}
   else if(q%7==0){plain[start]=rng()%20;graph[start]=Sum{plain[start],1};}
   Int expected=0,now=start;for(Int j=0;j<=k;++j){if(include||j)expected+=plain[now];now=next[now];}assert(graph.prod(start,k,include).sum==expected);
   Int l=rng()%60,r=l+rng()%40;auto products=graph.prod_range(start,l,r,include);for(Int j=l;j<r;++j)assert(products[j-l]==graph.prod(start,j,include));auto folded=graph.prod_range_fold(start,l,r,[](Int acc,Sum v){return acc+v.sum;},Int(0),include);Int wanted=0;for(auto v:products)wanted+=v.sum;assert(folded==wanted);
   Int cap=rng()%5000,limit=rng()%500,total=0,stop=limit;now=start;for(Int j=0;j<=limit;++j){total+=plain[now];if(total>cap){stop=j;break;}now=next[now];}assert(graph.move_while([&](Sum x){return x.sum<=cap;},start,limit)==stop);
   for(Int x=0;x<n;++x)assert(graph.get(x)==(Sum{plain[x],1}));

  }
  auto maxop=[](Int a,Int b){return std::max(a,b);};auto maxe=[](){return Int(-1000000000000LL);};auto map=[](Int a,Int b){return a+b;};auto add=[](Int a,Int b){return a+b;};auto id=[](){return Int(0);};auto maximum=initFunctionalGraph_with_lazy_op(next,plain,maxop,maxe,map,add,id);maximum.apply(0,n*3,Int(5));Int now=0;for(Int j=0;j<=n*3;++j){plain[now]+=5;now=next[now];}auto all=maximum.prod_reachable_idempotent_all();for(Int x=0;x<n;++x){Int best=maxe();now=x;for(Int j=0;j<n;++j){best=std::max(best,plain[now]);now=next[now];}assert(all[x]==best);}
 }
 std::vector<Int> edges{1,2,3,1,0,4};std::vector<std::string> text{"a","b","c","d","e","f"};
 auto concat=initFunctionalGraph_with_lazy_op(edges,text,std::plus<std::string>{},[](){return std::string{};},[](int c,std::string x){if(c>=0)std::fill(x.begin(),x.end(),char(c));return x;},[](int a,int b){return a>=0?a:b;},[](){return -1;});
 concat.apply(5,8,int('x'));Int at=5;for(Int i=0;i<=8;++i){text[at]="x";at=edges[at];}concat[2]=std::string("z");text[2]="z";
 for(Int start=0;start<6;++start){std::string want;at=start;for(Int k=0;k<50;++k){want+=text[at];assert(concat.prod(start,k)==want);at=edges[at];}auto products=concat.prod_range(start,0,50);for(Int k=0;k<50;++k)assert(products[k]==concat.prod(start,k));}
 auto graph=initFunctionalGraph_with_lazy_op<ExplicitTree>(std::vector<Int>{1,2,0},std::vector<Sum>{{1,1},{2,1},{3,1}});graph.apply(0,999999999999LL,Affine{1,1});assert(graph.get(0).sum==333333333335LL&&graph.get(1).sum==333333333335LL&&graph.get(2).sum==333333333336LL);
}
int main(){
 test_lazy();
 std::mt19937 rng(85325);
 auto empty=initFunctionalGraph(std::vector<Int>{});assert(empty.cycle.empty()&&empty.tree.N==1);
 for(int trial=0;trial<300;++trial){Int n=1+rng()%60;std::vector<Int> a(n),values(n);for(auto& x:a)x=rng()%n;for(auto& x:values)x=rng()%10;auto f=initFunctionalGraph(a);auto g=initFunctionalGraph_with_op(a,values,std::plus<Int>{},Int(0));
  std::vector<std::vector<Int>> distance(n,std::vector<Int>(n,-1));
  for(Int x=0;x<n;++x){Int v=x,steps=0;while(distance[x][v]<0){distance[x][v]=steps++;v=a[v];}assert(f.canmove_size(x)==steps);Int incoming=0;for(Int u=0;u<n;++u){Int now=u;for(Int k=0;k<n;++k){if(now==x){++incoming;break;}now=a[now];}}assert(f.reachable_to_size(x)==incoming);}
  for(Int x=0;x<n;++x)for(Int y=0;y<n;++y)assert(f.dist(x,y)==distance[x][y]);
  for(Int k:{0,1,2,17,60,200}){std::vector<Int> counts(n);for(Int x=0;x<n;++x){Int now=x;std::vector<Int> path{x};for(Int j=0;j<k;++j){now=a[now];path.push_back(now);}++counts[now];assert(f.movekth(x,k)==now&&f.walk(x,k)==path);}for(Int x=0;x<n;++x)assert(f.count_kth(x,k)==counts[x]);}
  auto [forest,compressed,roots]=f.compressed_forest();assert(forest.edge_count()==forest.len-Int(roots.size()));for(Int x=0;x<n;++x)assert(compressed[f.root(x)]==roots[f.cycle_number[f.root(x)]]);
  for(int query=0;query<80;++query){Int x=rng()%n,k=rng()%200;bool include=rng()%2;Int expected=0,now=x;for(Int j=0;j<=k;++j){if(include||j)expected+=values[now];now=a[now];}assert(g.prod(x,k,include)==expected);
   Int l=rng()%40,r=l+rng()%20;auto products=g.prod_range(x,l,r,include);for(Int j=l;j<r;++j)assert(products[j-l]==g.prod(x,j,include));auto sum=g.prod_range_fold(x,l,r,std::plus<Int>{},Int(0),include);assert(sum==std::accumulate(products.begin(),products.end(),Int(0)));
   Int cap=rng()%300,limit=rng()%200,stop=limit,total=0;now=x;for(Int j=0;j<=limit;++j){total+=values[now];if(total>cap){stop=j;break;}now=a[now];}assert(g.move_while([&](Int v){return v<=cap;},x,limit)==stop);
   if(query%5==0){Int vertex=rng()%n;values[vertex]=rng()%10;g[vertex]=values[vertex];}
  }
  auto maximum=initFunctionalGraph_with_op(a,values,[](Int x,Int y){return std::max(x,y);},Int(-1));auto all=maximum.prod_reachable_idempotent_all();for(Int x=0;x<n;++x){Int best=-1;for(Int y=0;y<n;++y)if(distance[x][y]>=0)best=std::max(best,values[y]);assert(all[x]==best);}
  std::vector<std::string> text(n);for(Int x=0;x<n;++x)text[x]=char('a'+x%26);auto concat=initFunctionalGraph_with_op(a,text,std::plus<std::string>{},std::string{});for(Int x=0;x<n;++x){std::string expected;Int now=x;for(Int k=0;k<150;++k){expected+=text[now];assert(concat.prod(x,k)==expected);now=a[now];}}
 }
 auto loop=initFunctionalGraph_with_op(std::vector<Int>{1,2,0},std::vector<Int>{1,2,3},std::plus<Int>{},Int(0));assert(loop.prod(0,999999999999LL)==1999999999999LL);assert(loop.move_while([](Int x){return x<=1000000000000LL;},0,1000000000000LL)==500000000000LL);
 std::vector<Int> chain(200000);std::iota(chain.begin(),chain.end(),Int(1));chain.back()=chain.size()-1;auto deep=initFunctionalGraph(chain);assert(deep.depth(0)==199999&&deep.count_kth(199999,200000)==200000&&deep.movekth(0,1000000000000LL)==199999);
 auto graph=initUnWeightedDirectedGraph(2);graph.add_edge(0,1);graph.add_edge(1,1);assert(initFunctionalGraph(graph).movekth(0,100)==1);
}
