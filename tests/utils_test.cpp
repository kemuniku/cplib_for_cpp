#include <cplib/utils/random_helper.hpp>
#include <cplib/collections/unionfind.hpp>
#include <cplib/utils/game.hpp>
#include <cplib/utils/grid_searcher.hpp>
#include <cplib/utils/implicit_dijkstra.hpp>
#include <set>
#include <random>
#include <cplib/utils/inversion_number.hpp>
#include <cplib/utils/seqidx.hpp>
#include <cplib/utils/bititers.hpp>
#include <cplib/utils/memo.hpp>
#include <cplib/utils/list_procs.hpp>
#include <cplib/utils/area_of_union_of_rectangles.hpp>
#include <cplib/utils/binary_search.hpp>
#include <cplib/utils/cumsum2d.hpp>
#include <cplib/utils/imos2d.hpp>
#include <cplib/utils/gridutils.hpp>
#include <cplib/utils/kth_element.hpp>
#include <cplib/utils/knapsack.hpp>
#include <cplib/utils/lis.hpp>
#include <cplib/utils/mo.hpp>
#include <cplib/utils/monotone_minima.hpp>
#include <cplib/utils/sequtils2D.hpp>
#include <cplib/utils/smawk.hpp>
#include <random>
using namespace cplib;
void utility_more_test();
void grid_implicit_test();
void game_test();
void random_test();
int main(){random_test();game_test();grid_implicit_test();
    utility_more_test();
    assert(meguru_bisect(Int(0),Int(100),[](Int x){return x*x<=200;})==14);
    assert(std::abs(meguru_bisect(0.0,10.0,[](double x){return x*x<=2;})-std::sqrt(2.0))<1e-9);
    auto imos=initImos2D(9,7);std::vector<std::vector<Int>> a(9,std::vector<Int>(7));std::mt19937 rng(19);
    for(int trial=0;trial<100;++trial){Int l=rng()%10,r=rng()%10,u=rng()%8,v=rng()%8,x=Int(rng()%30)-15;if(l>r)std::swap(l,r);if(u>v)std::swap(u,v);imos.rectangle_add(l,r,u,v,x);for(Int i=l;i<r;++i)for(Int j=u;j<v;++j)a[i][j]+=x;}
    assert(imos.build()==a);assert(imos.build()==a);auto sum=toCumSum2D(a);
    for(Int l=0;l<9;++l)for(Int r=l;r<=9;++r)for(Int u=0;u<7;++u)for(Int v=u;v<=7;++v){Int want=0;for(Int i=l;i<r;++i)for(Int j=u;j<v;++j)want+=a[i][j];assert(sum.query(l,r,u,v)==want);}
    for(Int h=0;h<40;++h)for(Int w=0;w<40;++w){
        auto cost=[&](Int r,Int c){return (3*r-c)*(3*r-c)+c;};auto better=[&](Int r,Int old,Int col){return cost(r,col)<cost(r,old);};
        std::vector<Int> want(h,-1);for(Int r=0;r<h;++r)for(Int c=0;c<w;++c)if(want[r]==-1 || better(r,want[r],c))want[r]=c;
        assert(smawk(h,w,better)==want);assert(monotoneMinima(h,w,better)==want);
    }
    for(int trial=0;trial<200;++trial){
        std::vector<Int> xs(1+rng()%30);for(Int& x:xs)x=Int(rng()%20)-10;auto sorted=xs;std::sort(sorted.begin(),sorted.end());
        for(Int k=0;k<Int(xs.size());++k){assert(kth_element(xs,k)==sorted[k]);auto copy=xs;assert(kth_element_break(copy,k)==sorted[k]);}
        std::vector<Int> dp(xs.size(),1);for(std::size_t i=0;i<xs.size();++i)for(std::size_t j=0;j<i;++j)if(xs[j]<xs[i])dp[i]=std::max(dp[i],dp[j]+1);
        Int length=*std::max_element(dp.begin(),dp.end());assert(lis(xs)==length);auto indices=restore_lis_index(xs);auto values=restore_lis(xs);assert(Int(indices.size())==length);for(Int i=0;i<length;++i){assert(xs[indices[i]]==values[i]);if(i)assert(indices[i-1]<indices[i] && values[i-1]<values[i]);}
        Int W=rng()%25;std::vector<KnapsackItem> items;std::vector<BoundedKnapsackItem> bounded;for(int i=0;i<8;++i){Int v=rng()%15,w=1+rng()%8,m=rng()%5;items.emplace_back(v,w);bounded.emplace_back(v,w,m);}
        Int best=0;for(int mask=0;mask<256;++mask){Int v=0,w=0;for(int i=0;i<8;++i)if(mask>>i&1){v+=items[i].first;w+=items[i].second;}if(w<=W)best=std::max(best,v);}
        assert(solve_01knapsack_NW(items,W)==best);assert(solve_01knapsack_NV(items,W)==best);assert(solve_01knapsack_meet_in_middle(items,W)==best);
        std::vector<Int> bdp(W+1);for(auto [v,w,m]:bounded){auto next=bdp;for(Int capacity=0;capacity<=W;++capacity)for(Int count=0;count<=m && count*w<=capacity;++count)next[capacity]=std::max(next[capacity],bdp[capacity-count*w]+count*v);bdp=next;}
        assert(solve_BoundedKnapsack(bounded,W)==bdp[W]);
        std::vector<Int> udp(W+1);for(Int capacity=0;capacity<=W;++capacity)for(auto [v,w]:items)if(w<=capacity)udp[capacity]=std::max(udp[capacity],udp[capacity-w]+v);assert(solve_UBknapsack_NW(items,W)==udp[W]);
    }
    std::vector<Int> xs(200);for(Int& x:xs)x=rng()%100;auto mo=initMo(200,300);std::vector<std::pair<Int,Int>> queries;std::vector<Int> answers(300);Int total=0;
    for(int i=0;i<300;++i){Int l=rng()%201,r=rng()%201;if(l>r)std::swap(l,r);mo.insert(l,r);queries.emplace_back(l,r);}
    auto add=[&](Int i){total+=xs[i];};auto del=[&](Int i){total-=xs[i];};mo.run(add,add,del,del,[&](Int i){answers[i]=total;});
    for(int i=0;i<300;++i){auto [l,r]=queries[i];assert(answers[i]==std::accumulate(xs.begin()+l,xs.begin()+r,Int(0)));}
    std::vector<std::string> grid{"abc","bca"};assert((gridfind(grid,'c')==std::pair<Int,Int>{0,2}));assert(gridfinds(grid,'a').size()==2);assert((to_pos(grid,4)==std::pair<Int,Int>{1,1}));
    Int neighbors=0;for(auto [i,j]:griditer(0,0,2,3)){assert(i+j==1);++neighbors;}assert(neighbors==2);
    assert((maxIndex(grid)==std::pair<Int,Int>{0,2}));assert((minIndex(grid)==std::pair<Int,Int>{0,0}));
}
void utility_more_test(){
    using namespace cplib;std::mt19937 rng(23);
    for(int q=0;q<500;++q){std::vector<Int> a(rng()%80);for(auto& x:a)x=Int(rng()%50)-25;Int want=0;for(Int i=0;i<Int(a.size());++i)for(Int j=i+1;j<Int(a.size());++j)want+=a[i]>a[j];assert(inversion_number(a)==want);auto x=newSeqWithIdx(Int(a.size()),[](Int i){return i;});assert(allItidx(x,[](Int i,Int v){return i==v;}));assert(!anyItIdx(x,[](Int i,Int v){return i!=v;}));auto y=mapItIdx(a,[](Int i,Int v){return i+v;});for(Int i=0;i<Int(a.size());++i)assert(y[i]==a[i]+i);}
    auto collect=[](auto range){std::vector<Int> v;for(Int x:range)v.push_back(x);return v;};
    for(Int n=0;n<=10;++n){for(Int r=0;r<=n;++r){std::vector<Int> want;for(Int b=0;b<(Int(1)<<n);++b)if(std::popcount(UInt(b))==r)want.push_back(b);assert(collect(bitcomb(n,r))==want);}for(Int bits=0;bits<(Int(1)<<n);++bits){std::vector<Int> sub,sup,single,standing;for(Int b=0;b<(Int(1)<<n);++b){if((b&bits)==b)sub.push_back(b);if((b&bits)==bits)sup.push_back(b);}assert(collect(bitsubseteq(bits))==sub);auto proper=sub;if(bits)proper.pop_back();assert(collect(bitsubset(bits))==proper);std::reverse(sub.begin(),sub.end());std::reverse(proper.begin(),proper.end());assert(collect(bitsubseteq_descending(bits))==sub&&collect(bitsubset_descending(bits))==proper);assert(collect(bitsuperseteq(bits,n))==sup);sup.erase(sup.begin());assert(collect(bitsuperset(bits,n))==sup);for(Int i=0;i<n;++i)if(bits&(Int(1)<<i)){single.push_back(Int(1)<<i);standing.push_back(i);}assert(collect(bitsingleton(bits))==single&&collect(standingbits(bits))==standing);}}
    assert(collect(bitcomb(63,63))==std::vector<Int>{std::numeric_limits<Int>::max()});assert(collect(standingbits(std::numeric_limits<Int>::min()))==std::vector<Int>{63});
    Int calls=0;auto fib=memo<Int(Int)>([&](auto& self,Int n)->Int{++calls;return n<2?n:self(n-1)+self(n-2);});assert(fib(40)==102334155&&calls==41);auto copy=fib;assert(copy(39)==63245986&&calls==41);auto choose=memo<Int(Int,Int)>([](auto& self,Int n,Int k)->Int{return k==0||k==n?1:self(n-1,k-1)+self(n-1,k);});assert(choose(20,10)==184756);auto square=memo<Int(Int)>([](Int x){return x*x;});assert(square(7)==49);
    DoublyLinkedNode<Int> a{1},b{2},c{3},d{4};DoublyLinkedList<Int> list{&a,&a};insert(list,&a,&c);insertPrev(list,&c,&b);insertPrev(list,&a,&d);std::vector<Int> forward,back;for(auto p=list.head;p;p=p->next)forward.push_back(p->value);for(auto p=list.tail;p;p=p->prev)back.push_back(p->value);assert((forward==std::vector<Int>{4,1,2,3}));std::reverse(back.begin(),back.end());assert(back==forward);
    for(int q=0;q<1000;++q){std::vector<std::tuple<Int,Int,Int,Int>> rectangles;bool grid[20][20]{};for(int i=0;i<int(rng()%30);++i){Int l=rng()%21,r=rng()%21,b=rng()%21,t=rng()%21;if(l>r)std::swap(l,r);if(b>t)std::swap(b,t);rectangles.emplace_back(l-10,b-10,r-10,t-10);for(Int x=l;x<r;++x)for(Int y=b;y<t;++y)grid[x][y]=true;}Int want=0;for(auto& row:grid)for(bool x:row)want+=x;assert(area_of_union_of_rectangles(rectangles)==want);}
    std::vector<std::tuple<Int,Int,Int,Int>> wide={{-3000000000LL,0,3000000000LL,2},{0,-1,4000000000LL,1}};assert(area_of_union_of_rectangles(wide)==17000000000LL);
}
void grid_implicit_test(){
    using P=std::pair<Int,Int>;std::mt19937 rng(829);auto grid=initGridSearcher();std::multiset<P> walls;
    for(int q=0;q<5000;++q){P p{rng()%30,rng()%30};if(rng()%2){grid.incl(p);walls.insert(p);}else{grid.excl(p);auto it=walls.find(p);if(it!=walls.end())walls.erase(it);}assert(grid.len()==Int(walls.size())&&grid.contains(p)==walls.contains(p));P at{rng()%30,rng()%30};std::array<std::optional<P>,4> expected;for(auto [i,j]:walls){if(j==at.second&&i<at.first&&(!expected[0]||expected[0]->first<i))expected[0]=P{i,j};if(j==at.second&&i>at.first&&(!expected[1]||expected[1]->first>i))expected[1]=P{i,j};if(i==at.first&&j<at.second&&(!expected[2]||expected[2]->second<j))expected[2]=P{i,j};if(i==at.first&&j>at.second&&(!expected[3]||expected[3]->second>j))expected[3]=P{i,j};}assert(grid.updownleftright(at)==expected);std::vector<P> filtered;for(auto p:expected)if(p)filtered.push_back(*p);assert(grid.updownleftright_get(at)==filtered);if(expected[0])++expected[0]->first;if(expected[1])--expected[1]->first;if(expected[2])++expected[2]->second;if(expected[3])--expected[3]->second;assert(grid.updownleftright_move(at)==expected);filtered.clear();for(auto p:expected)if(p)filtered.push_back(*p);assert(grid.updownleftright_move_get(at)==filtered);}
    for(int trial=0;trial<200;++trial){Int n=1+rng()%25;std::vector<std::vector<std::pair<std::string,Int>>> edges(n);std::vector<std::vector<Int>> distance(n,std::vector<Int>(n,INF64));for(Int i=0;i<n;++i){distance[i][i]=0;for(Int j=0;j<n;++j)if(rng()%4==0){Int w=rng()%20;edges[i].emplace_back(std::to_string(j),w);distance[i][j]=std::min(distance[i][j],w);}}for(Int k=0;k<n;++k)for(Int i=0;i<n;++i)for(Int j=0;j<n;++j)distance[i][j]=std::min(distance[i][j],distance[i][k]+distance[k][j]);auto adjacent=[&](const std::string& v){return edges[std::stoll(v)];};auto result=restore_implicit_dijkstra(std::string("0"),adjacent);auto costs=implicit_dijkstra(std::string("0"),adjacent);assert(costs==result.costs);for(Int i=0;i<n;++i){auto goal=std::to_string(i);assert(costs.contains(goal)==(distance[0][i]!=INF64));auto path=shortest_path_implicit_dijkstra(std::string("0"),goal,adjacent);assert(path.cost==distance[0][i]);assert(implicit_dijkstra_until(std::string("0"),adjacent,[&](const auto& v){return v==goal;})==distance[0][i]);if(distance[0][i]==INF64){assert(path.path.empty());continue;}assert(path.path.front()=="0"&&path.path.back()==goal&&costs.at(goal)==distance[0][i]);Int sum=0;for(std::size_t j=1;j<path.path.size();++j){Int w=INF64;for(auto [v,c]:adjacent(path.path[j-1]))if(v==path.path[j])w=std::min(w,c);sum+=w;}assert(sum==path.cost);}}
}
void game_test(){
    auto next=[](Int n){std::vector<Int> out;for(Int k=1;k<=2&&k<=n;++k)out.push_back(n-k);return out;};auto g=init_grundy(next);auto w=init_can_win(next);auto p=init_optimal_play(next);auto misere=init_can_win(next,true);for(Int i=0;i<100;++i){assert(g(i)==i%3&&w(i)==bool(i%3)&&misere(i)==(i%3!=1));assert(grundy(i,next)==g(i)&&can_win(i,next)==w(i));}assert((p(4)==OptimalPlay<Int>{true,{4,3,2,0}}));assert((p(3)==OptimalPlay<Int>{false,{3,2,0}}));
    std::mt19937 rng(45);for(bool terminal:{false,true})for(int trial=0;trial<100;++trial){constexpr Int N=30;std::array<std::vector<std::vector<Int>>,2> edges;for(auto& list:edges){list.resize(N);for(Int i=0;i<N;++i)for(Int j=0;j<i;++j)if(rng()%4==0)list[i].push_back(j);}auto next=[&](Int state,bool first){return edges[first][state];};auto score=[](Int x){return -x;};auto optimal=init_optimal_play(next,score,terminal);auto win=init_can_win(next,terminal);std::array<std::array<bool,N>,2> winning{};std::array<std::array<Int,N>,2> turns{},choice{};for(Int i=0;i<N;++i)for(int t=0;t<2;++t){winning[t][i]=terminal;turns[t][i]=0;choice[t][i]=i;bool first=true;for(Int j:edges[t][i]){bool w=!winning[!t][j];Int len=turns[!t][j]+1;if(first||(w&&!winning[t][i])||(w==winning[t][i]&&((w&&len<turns[t][i])||(!w&&len>turns[t][i])||(len==turns[t][i]&&j<choice[t][i])))){winning[t][i]=w;turns[t][i]=len;choice[t][i]=j;}first=false;}}for(Int i=0;i<N;++i){auto play=optimal(i);assert(play.is_win==winning[1][i]&&win(i)==winning[1][i]&&Int(play.states.size())==turns[1][i]+1);Int state=i,t=1;for(std::size_t k=1;k<play.states.size();++k){state=choice[t][state];t=!t;assert(play.states[k]==state);}assert(optimal_play(i,next,score,terminal)==play);}}
    auto cycle=[](Int x){return std::vector<Int>{x};};auto cyclic=init_grundy(cycle);for(int i=0;i<2;++i){bool threw=false;try{cyclic(0);}catch(const std::invalid_argument&){threw=true;}assert(threw);}bool threw=false;try{can_win(Int(0),cycle);}catch(const std::invalid_argument&){threw=true;}assert(threw);threw=false;try{optimal_play(Int(0),cycle);}catch(const std::invalid_argument&){threw=true;}assert(threw);
    Int calls=0;auto shared=init_grundy([&](Int n){++calls;return next(n);});shared(20);auto copy=shared;Int before=calls;copy(20);assert(calls==before);
}
void random_test(){
    cplib::seed(923);auto first=randomseq(100,closed_slice(-100,100));cplib::seed(923);assert(first==randomseq(100,closed_slice(-100,100)));for(Int n:{0,1,2,10,50,100,200}){auto v=randomseq(n,closed_slice(-100,100),true);std::set<Int> unique(v.begin(),v.end());assert(Int(v.size())==n&&Int(unique.size())==n);for(Int x:v)assert(-100<=x&&x<=100);auto sum=randomseq_from_sum(n,n?27:0);assert(Int(sum.size())==n&&std::accumulate(sum.begin(),sum.end(),Int(0))==(n?27:0));for(Int x:sum)assert(x>=0);auto ps=random_parenthesis_sequence(2*n);Int balance=0;for(Int x:ps){balance+=x;assert(balance>=0);}assert(balance==0);auto str=random_parenthesis_string(2*n);balance=0;for(char c:str){balance+=c=='('?1:-1;assert(balance>=0);}assert(balance==0);if(n){auto binary=random_binary_tree(n),tree=random_tree(n);for(const auto* g:{&binary,&tree}){assert(g->edge_count()==n-1&&g->len==n);auto uf=initUnionFind(n);for(auto e:g->edge_info){assert(!uf.issame(e.src,e.dst));uf.unite(e.src,e.dst);}}auto bits=random_01sequence(n,n/2);assert(Int(bits.size())==n&&std::accumulate(bits.begin(),bits.end(),Int(0))==n/2);}}
    auto primes=random_prime_sequence(30,closed_slice(2,1000),true);assert(std::set<Int>(primes.begin(),primes.end()).size()==30);for(Int x:primes)assert(isprime(x));for(Int x:random_prime_sequence(100,closed_slice(100,200)))assert(100<=x&&x<=200&&isprime(x));
    for(Int n:{1,2,3,10,50,3200}){Int m=std::min(n*(n-1)/2,n+20);auto simple=random_simple_graph(n,m),connected=random_connected_graph(n,m);for(const auto* g:{&simple,&connected}){std::set<std::pair<Int,Int>> edges;for(auto e:g->edge_info){assert(e.src!=e.dst);edges.emplace(std::min(e.src,e.dst),std::max(e.src,e.dst));}assert(Int(edges.size())==m&&g->edge_count()==m);}auto uf=initUnionFind(n);for(auto e:connected.edge_info)uf.unite(e.src,e.dst);for(Int i=0;i<n;++i)assert(uf.issame(0,i));}
    for(Int n=0;n<20;++n){Int m=n<3?n*(n-1)/2:3*n-6;auto g=random_planar_graph(n,m);assert(g.edge_count()==m&&is_planar_graph(g));}
    for(char c:random_string(1000,closed_slice('a','z'))){assert('a'<=c&&c<='z');}for(char c:random_string(1000,std::string("xyz")))assert(c=='x'||c=='y'||c=='z');
}
