#include <cplib/collections/combined.hpp>
#include <cplib/graph/combined.hpp>
#include <random>
#include <bitset>
#include <iostream>
using cplib::Int;namespace old=cplib;
std::mt19937_64 rng(238421);
template<int N> void bits(){
 for(int rep=0;rep<40;++rep){std::array<bool,N> values{};std::bitset<N> a,b;for(int i=0;i<N;++i){a[i]=values[i]=rng()%2;b[i]=rng()%2;}
 auto x=cplib::initBitSet<N>(values);cplib::BitSet<N> y;for(int i=0;i<N;++i)y[i]=bool(b[i]);
 assert(x.str()==a.to_string());assert(x.popcount()==Int(a.count()));assert((~x).str()==(~a).to_string());assert((x&y).str()==(a&b).to_string());assert((x|y).str()==(a|b).to_string());assert((x^y).str()==(a^b).to_string());assert(x.andpopcount(y)==Int((a&b).count()));assert(x.orpopcount(y)==Int((a|b).count()));assert(x.xorpopcount(y)==Int((a^b).count()));
 std::vector<Int> indices;for(int i=0;i<N;++i)if(a[i])indices.push_back(i);std::vector<Int> actual;for(Int i:x)actual.push_back(i);assert(actual==indices);assert(x.lowestBit()==(indices.empty()?-1:indices[0]));assert(cplib::BitSet<N>::fromIndexes(indices).str()==a.to_string());
 for(int shift:{0,1,7,8,31,32,63,64,65,127,128,255,N,N+1}){assert((x<<shift).str()==(a<<shift).to_string());assert((x>>shift).str()==(a>>shift).to_string());}
 x&=y;assert(x.str()==(a&b).to_string());x|=y;assert(x.str()==b.to_string());x^=y;assert(!x.popcount());
 }
}
template<bool S,bool D,bool W> void graph(){old::BasicGraph<Int,S,D,W> g(10);std::vector<std::vector<std::pair<std::int32_t,Int>>> expected(10);for(int i=0;i<100;++i){Int u=rng()%10,v=rng()%10,w=W?Int(rng()%100):1;if constexpr(W)g.add_edge(u,v,w);else g.add_edge(u,v);expected[u].emplace_back(v,w);if(!D)expected[v].emplace_back(u,w);}if constexpr(S)g.build();for(int i=0;i<10;++i){auto range=g.to_and_cost(i);assert(std::vector(range.begin(),range.end())==expected[i]);Int j=0;for(auto e:g[i]){if constexpr(W)assert(e==std::make_pair(Int(expected[i][j].first),expected[i][j].second));else assert(e==expected[i][j].first);++j;}}}
void functional(){
 for(int trial=0;trial<200;++trial){int n=1+rng()%80;std::vector<Int> next(n),values(n);for(int i=0;i<n;++i){next[i]=rng()%n;values[i]=rng()%10+1;}
 auto f=old::initFunctionalGraph(next);auto op=old::initFunctionalGraph_with_op(next,values,std::plus<Int>{},Int(0));auto mx=old::initFunctionalGraph_with_op(next,values,[](Int a,Int b){return std::max(a,b);},Int(0));
 for(int rep=0;rep<100;++rep){Int x=rng()%n,k=rng()%200,target=x,product=values[x];for(Int i=0;i<k;++i){target=next[target];product+=values[target];}assert(f.movekth(x,k)==target);assert(op.prod(x,k)==product);Int incoming=0,reachable=0;for(Int start=0;start<n;++start){Int z=start;for(Int i=0;i<k;++i)z=next[z];incoming+=z==x;z=start;bool reaches=false;for(Int i=0;i<n;++i){reaches|=z==x;z=next[z];}reachable+=reaches;}assert(f.count_kth(x,k)==incoming&&f.reachable_to_size(x)==reachable);std::vector<bool> seen(n);Int z=x,moves=0;while(!seen[z]){seen[z]=true;++moves;z=next[z];}assert(f.canmove_size(x)==moves);auto cycle=f.get_cycle(x);assert(!cycle.empty());for(std::size_t i=0;i<cycle.size();++i)assert(next[cycle[i]]==cycle[(i+1)%cycle.size()]);Int y=rng()%n,distance=-1;z=x;for(Int i=0;i<n;++i){if(z==y){distance=i;break;}z=next[z];}assert(f.dist(x,y)==distance);
 Int bound=rng()%500,limit=rng()%150,sum=values[x],used=0,now=x;if(sum<=bound)while(used<limit&&sum+values[next[now]]<=bound){now=next[now];sum+=values[now];++used;}assert(op.move_while([&](Int z){return z<=bound;},x,limit)==used);
 if(rep%11==0){Int v=rng()%n,w=rng()%10+1;values[v]=w;op[v]=w;mx[v]=w;}
 }
 auto maximum=mx.prod_reachable_idempotent_all();for(int x=0;x<n;++x){Int ans=0,now=x;for(int i=0;i<n;++i){ans=std::max(ans,values[now]);now=next[now];}assert(maximum[x]==ans);}
 }
 // 旧版は失敗直前の移動回数を返す（現行版の最初の失敗位置とは異なる）。
 auto f=old::initFunctionalGraph_with_op(std::vector<Int>{1,2,2},std::vector<Int>{2,3,4},std::plus<Int>{},Int(0));assert(f.move_while([](Int s){return s<=5;},0,100)==1);
}
void hld(){for(int trial=0;trial<100;++trial){Int n=2+rng()%200;auto g=cplib::initUnWeightedUnDirectedStaticGraph(n);std::vector<std::vector<Int>> adj(n);for(Int i=1;i<n;++i){Int p=rng()%i;g.add_edge(i,p);adj[i].push_back(p);adj[p].push_back(i);}g.build();Int root=rng()%n;auto h=cplib::initHld(g,root);std::vector<Int> parent(n,-1),depth(n),order{root};for(std::size_t i=0;i<order.size();++i)for(Int v:adj[order[i]])if(v!=parent[order[i]]){parent[v]=order[i];depth[v]=depth[order[i]]+1;order.push_back(v);}assert(h.P==parent&&h.D==depth);for(int i=0;i<100;++i){Int a=rng()%n,b=rng()%n,d=rng()%n,u=a,v=b;std::vector<Int> left,right;while(depth[u]>depth[v]){left.push_back(u);u=parent[u];}while(depth[v]>depth[u]){right.push_back(v);v=parent[v];}while(u!=v){left.push_back(u);right.push_back(v);u=parent[u];v=parent[v];}assert(h.lca(a,b)==u);left.push_back(u);std::reverse(right.begin(),right.end());left.insert(left.end(),right.begin(),right.end());assert(h.la(a,b,d)==(d<Int(left.size())?left[d]:-1));std::vector<Int> actual;for(auto [l,r,reverse]:h.pathWithDirection(a,b))for(Int j=l;j<r;++j)actual.push_back(h.toVtx(reverse?n-1-j:j));assert(actual==left);}std::vector<Int> chosen{0,n-1,n/2};auto aux=cplib::initAuxiliaryWeightedTree(h,chosen);for(auto u:chosen)assert(aux.toi.contains(u));}}
int main(){bits<0>();bits<1>();bits<63>();bits<64>();bits<65>();bits<255>();bits<256>();bits<257>();bits<4097>();graph<0,0,0>();graph<0,0,1>();graph<0,1,0>();graph<0,1,1>();graph<1,0,0>();graph<1,0,1>();graph<1,1,0>();graph<1,1,1>();functional();hld();auto tab=old::initWeightedDirectedTableGraph(std::vector<std::string>{"a","b"});tab.add_edge("a","b",9);for(auto [v,w]:tab["a"])assert(v=="b"&&w==9);auto st=old::initSegmentTree(std::vector<Int>{1,2,3},std::plus<Int>{},Int(0));assert(st.get_all()==6&&st.min_left(3,[](Int x){return x<=5;})==1);st[1]=9;assert(st.max_right(0,[](Int x){return x<=10;})==2);std::cout<<"combined tests passed\n";}
