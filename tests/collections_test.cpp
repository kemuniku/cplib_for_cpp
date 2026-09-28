#include <cplib/collections/root_rangesum.hpp>
#include <cplib/collections/range_linear_add_range_min.hpp>
#include <cplib/collections/hashtable.hpp>
#include <cplib/collections/rangeset.hpp>
#include <cplib/collections/tatyamset.hpp>
#include <cplib/collections/avlset.hpp>
#include <cplib/collections/convex_hull_trick.hpp>
#include <cplib/collections/slopetrick.hpp>
#include <cplib/collections/lichaotree.hpp>
#include <cplib/collections/convex_hull_trick_monotone.hpp>
#include <cplib/collections/convex_hull_trick_monotone_slope.hpp>
#include <cplib/collections/binary_trie.hpp>
#include <cplib/collections/persistent_binary_trie.hpp>
#include <cplib/collections/deletable_heapqueue.hpp>
#include <cplib/collections/topk_sum_heapq.hpp>
#include <cplib/collections/defaultdict.hpp>
#include <set>
#include <cplib/collections/unionfind.hpp>
#include <cplib/collections/weightedunionfind.hpp>
#include <cplib/collections/rollback_unionfind.hpp>
#include <cplib/collections/ppunionfind.hpp>
#include <cplib/collections/rootvalue_unionfind.hpp>
#include <cplib/collections/fenwick.hpp>
#include <cplib/collections/fenwick2d.hpp>
#include <cplib/collections/SWAG.hpp>
#include <cplib/collections/persistent_array.hpp>
#include <cplib/collections/persistent_unionfind.hpp>
#include <cplib/collections/parallel_unionfind.hpp>
#include <deque>
#include <set>
#include <random>
using namespace cplib;
void persistent_test();
void heap_dict_test();
void binary_trie_test();
void convex_test();
void avl_test();
void range_bucket_test();
void hash_test();
void root_linear_test();
int main(){root_linear_test();hash_test();range_bucket_test();avl_test();convex_test();heap_dict_test();binary_trie_test();
    persistent_test();
    std::mt19937 rng(781);Int n=35;
    auto u=initUnionFind(n);auto weighted=initWeightedUnionFind(n);auto partial=initPartialPersistentUnionFind(n);auto rollback=initRollbackUnionFind(n);
    std::vector<UnionFind> states{u};std::vector<Int> potentials(n);for(Int& p:potentials)p=Int(rng()%100)-50;
    auto values=initRootValueUnionFind<Int>(n,[](Int& x,Int& y){x+=y;},Int(1));
    auto factory=initRootValueUnionFind<std::set<Int>>(n,[](auto& x,auto& y){x.insert(y.begin(),y.end());},[](){return std::set<Int>{};});
    std::vector<std::pair<Int,Int>> edges;
    for(Int t=0;t<100;++t){
        Int x=rng()%n,y=rng()%n;edges.emplace_back(x,y);u.unite(x,y);values.unite(x,y);factory.unite(x,y);partial.unite(x,y,t);rollback.unite(x,y);assert(weighted.unite(x,y,potentials[y]-potentials[x]));states.push_back(u);
        assert(u.count==rollback.count());
        for(Int v=0;v<n;++v){assert(u.siz(v)==values.get(v));assert(u.siz(v)==partial.size(v));}
        for(Int a=0;a<n;++a)for(Int b=0;b<n;++b)if(u.issame(a,b)){assert(weighted.diff(a,b)==potentials[b]-potentials[a]);assert(!weighted.unite(a,b,potentials[b]-potentials[a]+1));}
    }
    for(Int t=-1;t<100;++t)for(Int x=0;x<n;++x){assert(partial.size(x,t)==states[t+1].siz(x));for(Int y=0;y<n;++y)assert(partial.issame(x,y,t)==states[t+1].issame(x,y));}
    for(Int x=0;x<n;++x)for(Int y=0;y<n;++y){Int when=-2;for(Int t=-1;t<100;++t)if(states[t+1].issame(x,y)){when=t;break;}assert(partial.when_unite(x,y)==when);}
    for(Int t=99;t>=0;--t){rollback.undo();assert(rollback.count()==states[t].count);for(Int x=0;x<n;++x)for(Int y=0;y<n;++y)assert(rollback.issame(x,y)==states[t].issame(x,y));}
    rollback.unite(0,1);rollback.snapshot();rollback.unite(1,2);rollback.rollback();assert(rollback.issame(0,1)&&!rollback.issame(0,2));rollback.clear_snapshot();rollback.rollback();assert(rollback.count()==n);
    std::vector<Int> a(5000);for(Int& x:a)x=Int(rng()%100)-50;auto f=initFenwickTree(a);
    assert(Int(f[from_end(1)])==a.back());f[from_end(2)]=Int(123);a[a.size()-2]=123;
    assert(f[closed_slice(Int(0),from_end(1))]==std::accumulate(a.begin(),a.end(),Int(0)));
    for(int trial=0;trial<1000;++trial){Int p=rng()%a.size(),delta=Int(rng()%40)-20;f.add(p,delta);a[p]+=delta;Int l=rng()%(a.size()+1),r=rng()%(a.size()+1);if(l>r)std::swap(l,r);assert(f.get(l,r)==std::accumulate(a.begin()+l,a.begin()+r,Int(0)));assert(f.prefix(r)==std::accumulate(a.begin(),a.begin()+r,Int(0)));f[p]=delta;a[p]=delta;assert(Int(f[p])==a[p]);}
    std::vector<std::vector<Int>> positions(12,std::vector<Int>{-5,0,3,9});auto f2=initFenwick2D(positions);std::vector<std::tuple<Int,Int,Int>> updates;
    for(int trial=0;trial<300;++trial){Int x=rng()%12,y=positions[x][rng()%4],v=Int(rng()%30)-15;f2.add(x,y,v);updates.emplace_back(x,y,v);Int l=rng()%13,r=rng()%13,yl=Int(rng()%20)-7,yr=Int(rng()%20)-7;if(l>r)std::swap(l,r);if(yl>yr)std::swap(yl,yr);Int want=0;for(auto [px,py,pv]:updates)if(l<=px&&px<r&&yl<=py&&py<yr)want+=pv;assert(f2.get(l,r,yl,yr)==want);}
    auto op=[](const std::string& a,const std::string& b){return a+b;};auto swag=initSWAG(op,std::string());std::deque<std::string> dq,qq;
    for(int trial=0;trial<5000;++trial){
        int action=rng()%4;std::string x(1,char('a'+rng()%26));if(dq.empty()||action==0){swag.addFirst(x);dq.push_front(x);}else if(action==1){swag.addLast(x);dq.push_back(x);}else if(action==2){assert(swag.popFirst()==dq.front());dq.pop_front();}else{assert(swag.popLast()==dq.back());dq.pop_back();}
        std::string want;for(auto& x:dq)want+=x;assert(swag.fold()==want);assert(swag.len()==Int(dq.size()));for(Int i=0;i<swag.len();++i)assert(swag[i]==dq[i]);
        if(!dq.empty())assert(swag[from_end(1)]==dq.back());
    }
    auto sum=[](Int x,Int y){return x+y;};std::vector<Int> input{1,2,3,4,5,6};auto ok=[](Int x){return x<=7;};std::vector<Int> want;for(Int l=0;l<Int(input.size());++l){Int r=l,s=0;while(r<Int(input.size())&&s+input[r]<=7)s+=input[r++];want.push_back(r);}assert(get_maxrights(input,sum,Int(0),ok)==want);
}
void persistent_test(){
    std::mt19937 rng(18);
    for(Int n:{0,1,2,31,32,33,1000}){
        std::vector<Int> a(n);std::iota(a.begin(),a.end(),0);auto initial=initPersistentArray(a);assert(initial.toseq()==a);std::vector<decltype(initial)> versions{initial};std::vector<std::vector<Int>> refs{a};
        for(int q=0;n&&q<200;++q){Int v=rng()%versions.size(),i=rng()%n,x=rng();auto b=refs[v];b[i]=x;versions.push_back(versions[v].change_value(i,x));refs.push_back(b);v=rng()%versions.size();assert(versions[v].toseq()==refs[v]);for(int j=0;j<20;++j){i=rng()%n;assert(versions[v][i]==refs[v][i]);}}
        auto binary=initPersistentArray<Int,1>(a);assert(binary.toseq()==a);if(n){auto b=binary.change_value(n-1,Int(-1));assert(b[n-1]==-1&&binary[n-1]==n-1);}
    }
    for(Int n:{0,1,17,257,1025}){
        std::vector<PersistentUnionFind> versions{initPersistentUnionFind(n)};std::vector<UnionFind> refs{initUnionFind(n)};
        for(int q=0;n&&q<1000;++q){Int v=rng()%versions.size(),x=rng()%n,y=rng()%n;versions.push_back(versions[v].unite(x,y));auto ref=refs[v];ref.unite(x,y);refs.push_back(ref);v=rng()%versions.size();assert(versions[v].count==refs[v].count);for(int k=0;k<20;++k){x=rng()%n;y=rng()%n;assert(versions[v].issame(x,y)==refs[v].issame(x,y));assert(versions[v].siz(x)==refs[v].siz(x));}}
    }
    for(Int n:{0,1,3,4,5,16,17,65,500}){auto p=initParallelUnionFind(n);auto ref=initUnionFind(n);for(int q=0;q<500;++q){Int a=rng()%(n+1),b=rng()%(n+1),length=rng()%(n-std::max(a,b)+1),before=ref.count,calls=0;Int got=p.unite(a,b,length,[&](Int x,Int y){assert(p.root(x)==x&&p.root(y)==y&&x!=y);assert(p.siz(x)>=p.siz(y));++calls;});for(Int i=0;i<length;++i)ref.unite(a+i,b+i);assert(got==before-ref.count&&calls==got&&p.count()==ref.count);for(int k=0;n&&k<30;++k){Int x=rng()%n,y=rng()%n;assert(p.issame(x,y)==ref.issame(x,y));assert(p.siz(x)==ref.siz(x));}assert(Int(p.roots().size())==p.count());}auto copy=p.copy();if(n>1){copy.unite(0,1);assert(p.count()==ref.count);}}
    auto one=initParallelUnionFind(3);assert(unite(one,0,1));assert(!unite(one,0,1));assert(unite(one,Int(0),Int(1),Int(2))==1);
}

void heap_dict_test(){
    std::mt19937 rng(40);auto heap=initDeletableHeapQueue<Int>();std::multiset<Int> values;
    for(int i=0;i<20000;++i){if(values.empty()||rng()%2==0){Int x=Int(rng()%100)-50;heap.push(x);values.insert(x);}else if(rng()%2){auto it=values.begin();std::advance(it,rng()%values.size());heap.erase(*it);values.erase(it);}else{assert(heap.pop()==*values.begin());values.erase(values.begin());}assert(heap.len()==Int(values.size()));if(!values.empty())assert(heap[0]==*values.begin());}
    std::vector<Int> initial(100);for(auto& x:initial)x=Int(rng()%100)-50;auto top=initTopKHeapq(initial,40);values={initial.begin(),initial.end()};for(int i=0;i<10000;++i){if(values.empty()||rng()%3==0){Int x=Int(rng()%100)-50;top.incl(x);values.insert(x);}else if(rng()%2){if(Int(values.size())<=top.k)top.setK(values.size()-1);auto it=values.begin();std::advance(it,rng()%values.size());top.excl(*it);values.erase(it);}else top.setK(rng()%(values.size()+1));Int sum=0,large=0,j=0;for(auto it=values.rbegin();it!=values.rend();++it){sum+=*it;if(j++<top.k)large+=*it;}assert(top.sm==sum&&top.topk==large);}
    auto d=initDefaultDict<std::string,Int>(7);const auto& cd=d;assert(cd["absent"]==7&&d.len()==0);d["one"]+=3;assert(cd["one"]==10&&d.len()==1);auto other=initDefaultDict<std::string,Int>(99);other["one"]=10;assert(d==other&&d.hash()==other.hash());d["two"]=20;Int count=0;for(const auto& [k,v]:d.pairs()){assert(k=="one"||k=="two");assert(v==10||v==20);++count;}assert(count==2);for(auto& [k,v]:d.mpairs())v+=1;Int total=0;for(auto v:d.values())total+=v;assert(total==32);count=0;for(const auto& key:d.keys()){assert(d.hasKey(key));++count;}assert(count==2);Int value=55;assert(!d.pop("none",value)&&value==55);assert(d.take("one",value)&&value==11);d.del("two");assert(d.len()==0);d["x"]=1;d.clear();assert(d.len()==0&&std::as_const(d)["x"]==7);auto from=toDefaultDict(std::vector<std::pair<Int,Int>>{{1,3},{1,5},{2,7}},Int(9));assert(from[1]==5&&from.len()==2);
}

void binary_trie_test(){
    std::mt19937 rng(42);auto t=initBineryTrie(10);std::vector<Int> counts(1024);Int length=0;
    auto check=[&](const auto& trie,const std::vector<Int>& counts){std::vector<Int> sorted;for(Int i=0;i<1024;++i){assert(trie.count(i)==counts[i]);for(Int k=0;k<counts[i];++k)sorted.push_back(i);}Int x=rng()%1024;assert(trie.lowerBound(x)==std::lower_bound(sorted.begin(),sorted.end(),x)-sorted.begin());assert(trie.upperBound(x)==std::upper_bound(sorted.begin(),sorted.end(),x)-sorted.begin());assert(trie.get_kth(sorted.size())==-1);if(!sorted.empty()){Int k=rng()%sorted.size();assert(trie[k]==sorted[k]);Int mask=rng()%1024;std::sort(sorted.begin(),sorted.end(),[=](Int a,Int b){return (a^mask)<(b^mask);});assert(trie.get_kth(k,mask)==sorted[k]);}return sorted;};
    for(int q=0;q<3000;++q){Int x=rng()%1024;if(counts[x]&&rng()%2){t.excl(x);--counts[x];--length;}else{Int v=1+rng()%3;t.incl(x,v);counts[x]+=v;length+=v;}assert(t.len()==length);if(q%13==0){auto sorted=check(t,counts);if(length)assert(t[from_end(1)]==*std::max_element(sorted.begin(),sorted.end()));Int mask=rng()%1024,lower=0,upper=0;for(Int v:sorted){lower+=(v^mask)<x;upper+=(v^mask)<=x;}assert(t.lowerBound(x,mask)==lower&&t.upperBound(x,mask)==upper);}}
    std::vector<PersistentBinaryTrie> versions{initPersistentBineryTrie(10)};std::vector<std::vector<Int>> states{std::vector<Int>(1024)};for(int q=0;q<1000;++q){Int base=rng()%versions.size(),x=rng()%1024,v=rng()%5;auto c=states[base];PersistentBinaryTrie next=versions[base];if(q%3==0){next=next.set_value(x,v);c[x]=v;}else if(c[x]&&q%2){next=next.excl(x);--c[x];}else{next=next.incl(x,v);c[x]+=v;}versions.push_back(next);states.push_back(c);check(versions[base],states[base]);check(versions.back(),states.back());std::vector<std::pair<Int,Int>> rle;for(Int i=0;i<1024;++i)if(c[i])rle.emplace_back(i,c[i]);assert(next.RLE()==rle);}
    auto copy=t;Int before=t.count(0);copy.incl(0);assert(t.count(0)==before&&copy.count(0)==before+1);
}

void convex_test(){
    std::mt19937 rng(782);
    for(int repeat=0;repeat<100;++repeat){
        auto f=initSlopeTrick(3);constexpr int W=2000;std::vector<Int> values(2*W+1,3);
        for(int step=0;step<100;++step){int op=rng()%7;Int a=Int(rng()%21)-10;
            if(op<3){if(op==0)f.add_abs(a);if(op==1)f.add_a_minus_x(a);if(op==2)f.add_x_minus_a(a);for(int x=-W;x<=W;++x)values[x+W]+=op==0?std::abs(x-a):op==1?std::max(a-x,Int(0)):std::max(x-a,Int(0));}
            else if(op==3){f.add_all(a);for(auto& v:values)v+=a;}
            else if(op==4){f.clearL();for(int i=2*W;i-->0;)values[i]=std::min(values[i],values[i+1]);}
            else if(op==5){f.clearR();for(int i=1;i<=2*W;++i)values[i]=std::min(values[i],values[i-1]);}
            else {Int b=a+rng()%11;f.shift(a,b);auto old=values;for(int x=-W;x<=W;++x){values[x+W]=INF64;for(Int y=x-b;y<=x-a;++y)if(y>=-W&&y<=W)values[x+W]=std::min(values[x+W],old[y+W]);}}
            for(int x=-100;x<=100;++x){assert(f.get_value(x)==values[x+W]);}assert(f.min()==*std::min_element(values.begin()+W-1000,values.begin()+W+1001));
        }
    }
    for(int repeat=0;repeat<100;++repeat){std::vector<Int> xs;for(Int x=-50;x<=50;++x)if(rng()%2)xs.push_back(x);auto tree=initLiChaoTree(xs);std::vector<Int> best(xs.size(),INF64);for(int q=0;q<200;++q){Int a=Int(rng()%101)-50,b=Int(rng()%301)-150,l=-60,r=60;bool segment=rng()%2;if(segment){l=Int(rng()%121)-60;r=Int(rng()%121)-60;if(l>r)std::swap(l,r);tree.add_segment(a,b,l,r);}else tree.add_line(a,b);for(std::size_t i=0;i<xs.size();++i){if(l<=xs[i]&&xs[i]<r)best[i]=std::min(best[i],a*xs[i]+b);assert(tree.get_min(xs[i])==best[i]);}}}
    for(bool si:{false,true})for(bool xi:{false,true})for(int repeat=0;repeat<100;++repeat){auto a=initConvexHullTrickMonotone(si,xi);auto b=initConvexHullTrickMonotoneSlope(si);std::vector<CHTLine> lines;Int slope=si?-100:100,x=xi?-100:100;for(int q=0;q<100;++q){slope+=(si?1:-1)*Int(rng()%3);Int intercept=Int(rng()%1001)-500;lines.push_back({slope,intercept});a.add_line(slope,intercept);b.add_line(slope,intercept);x+=(xi?1:-1)*Int(rng()%3);Int expected=INF64;for(auto line:lines)expected=std::min(expected,chtAnswer(chtValue(line,x)));assert(a.get_min(x)==expected);Int y=Int(rng()%301)-150;expected=INF64;for(auto line:lines)expected=std::min(expected,chtAnswer(chtValue(line,y)));assert(b.get_min(y)==expected);}}
}
template<int Old> void avl_version_test(){
    std::mt19937 rng(19);auto m=initAvlSortedMultiSet<Int>();auto s=initAvlSortedSet<Int>();std::multiset<Int> expected;std::set<Int> unique;
    for(int q=0;q<30000;++q){Int x=Int(rng()%1000)-500;if(rng()%2){m.incl(x);s.incl(x);expected.insert(x);unique.insert(x);}else{auto it=expected.find(x);bool has=it!=expected.end();assert(m.excl(x)==has);if(has)expected.erase(it);assert(s.excl(x)==bool(unique.erase(x)));}assert(m.len()==Int(expected.size())&&s.len()==Int(unique.size()));assert(m.count(x)==Int(expected.count(x)));assert(m.lowerBound(x)==std::distance(expected.begin(),expected.lower_bound(x)));assert(m.upperBound(x)==std::distance(expected.begin(),expected.upper_bound(x)));if(q%100==0){std::vector<Int> v;for(Int y:m)v.push_back(y);assert(std::equal(v.begin(),v.end(),expected.begin(),expected.end()));for(Int i=0;i<m.len();++i){auto n=get(m.root,i);assert(n->key==v[i]&&index(n)==i&&rootOf(n)==m.root);assert(m[i]==v[i]);if(i+1<m.len())assert(next(n)->key==v[i+1]);if(i)assert(prev(n)->key==v[i-1]);}auto copy=m;copy.incl(10000);assert(!m.contains(10000));}}
    while(m.len()){assert(m.pop()==*expected.rbegin());expected.erase(std::prev(expected.end()));}
    for(int repeat=0;repeat<100;++repeat){auto cht=initConvexHullTrick();std::vector<CHTLine> lines;for(int q=0;q<300;++q){Int a=Int(rng()%201)-100,b=Int(rng()%2001)-1000;cht.add_line(a,b);lines.push_back({a,b});for(int k=0;k<4;++k){Int x=Int(rng()%2001)-1000,answer=INF64;for(auto line:lines)answer=std::min(answer,chtAnswer(chtValue(line,x)));assert(cht.get_min(x)==answer);}}}
}

void avl_test(){avl_version_test<false>();}
void range_bucket_test(){
    std::mt19937 rng(723);auto ranges=initRangeSet(Int(0));std::vector<Int> values(200);for(int q=0;q<10000;++q){Int l=rng()%200,r=rng()%200;if(l>r)std::swap(l,r);++r;Int value=rng()%5;ranges.update(l,r,value);std::fill(values.begin()+l,values.begin()+r,value);for(Int i=0;i<200;++i){auto [a,b,c]=ranges.get_segment(i);assert(a<=i&&i<b&&c==values[i]);}Int prev=std::numeric_limits<Int>::min(),last=-1;for(auto [a,b,c]:ranges.st){assert(a==prev&&a<b&&c!=last);prev=b;last=c;}assert(prev==std::numeric_limits<Int>::max());}
    std::vector<Int> initial(1000);for(auto& x:initial)x=rng()%500;auto buckets=initSortedMultiset(initial);std::multiset<Int> expected(initial.begin(),initial.end());for(int q=0;q<20000;++q){Int x=rng()%500;switch(rng()%3){case 0:buckets.incl(x);expected.insert(x);break;case 1:{auto it=expected.find(x);assert(buckets.excl(x)==(it!=expected.end()));if(it!=expected.end())expected.erase(it);break;}default:if(!expected.empty()){Int i=rng()%expected.size();auto it=expected.begin();std::advance(it,i);assert(buckets[i]==*it&&buckets[i-Int(expected.size())]==*it);assert(buckets.pop(i)==*it);expected.erase(it);}}
        assert(buckets.len()==Int(expected.size())&&buckets.count(x)==Int(expected.count(x)));assert(buckets.index(x)==std::distance(expected.begin(),expected.lower_bound(x)));auto lb=expected.lower_bound(x),ub=expected.upper_bound(x);assert(buckets.ge(x)==(lb==expected.end()?std::optional<Int>{}:*lb));assert(buckets.gt(x)==(ub==expected.end()?std::optional<Int>{}:*ub));assert(buckets.lt(x)==(lb==expected.begin()?std::optional<Int>{}:*std::prev(lb)));assert(buckets.le(x)==(ub==expected.begin()?std::optional<Int>{}:*std::prev(ub)));if(q%100==0){std::vector<Int> actual;for(auto y:buckets)actual.push_back(y);assert(std::equal(actual.begin(),actual.end(),expected.begin(),expected.end()));}}
}
struct CollisionHash {std::size_t operator()(Int x)const{return std::size_t(x%4);}};
void hash_test(){
    std::mt19937 rng(543);HashSet<Int,CollisionHash> set;HashTable<Int,Int,CollisionHash> table;std::set<Int> keys;std::unordered_map<Int,Int> values;
    for(int q=0;q<30000;++q){Int key=rng()%300,value=rng()%10000;if(rng()%2){set.incl(key);keys.insert(key);if(rng()%2)table[key]=value;else table.incl({key,value});values[key]=value;}else{set.excl(key);keys.erase(key);table.del(key);values.erase(key);}assert(set.len()==Int(keys.size())&&table.len()==Int(values.size()));assert(set.contains(key)==keys.contains(key)&&table.contains(key)==values.contains(key));if(values.contains(key)){assert(table[key]==values[key]);table[key]+=1;++values[key];assert(std::as_const(table)[key]==values[key]);}if(q%100==0){std::set<Int> actual;for(Int k:set)actual.insert(k);assert(keys==actual);std::unordered_map<Int,Int> actual_values;for(auto [k,v]:table.pairs())actual_values[k]=v;assert(values==actual_values);Int sum=0,expected_sum=0;for(auto v:table.values())sum+=v;for(auto [k,v]:values)expected_sum+=v;assert(sum==expected_sum);}}table.clear();assert(table.len()==0);for(auto [k,v]:table.pairs()){(void)k;(void)v;assert(false);}
}
void root_linear_test(){
    std::mt19937 rng(924);for(Int n:{0,1,2,3,15,16,17,100,513}){std::vector<Int> a(n);for(auto& x:a)x=rng()%100;auto linear=initRangeLinearAddRangeMin(a);for(int q=0;q<3000;++q){Int l=rng()%(n+1),r=rng()%(n+1);if(l>r)std::swap(l,r);if(rng()%2){Int b=Int(rng()%101)-50,c=Int(rng()%101)-50;linear.add(l,r,b,c);for(Int i=l;i<r;++i)a[i]+=b*i+c;}else{Int expected=l==r?std::numeric_limits<Int>::max():*std::min_element(a.begin()+l,a.begin()+r);assert(linear.prod(l,r)==expected);if(l<r)assert(linear[closed_slice(l,r-1)]==expected);}if(n){Int i=rng()%n;assert(linear[i]==a[i]&&linear[from_end(n-i)]==a[i]);}}
        for(Int bs:{0,1,5,31,1000}){for(auto& x:a)x=rng()%100;auto sums=initrangesum(a,bs);for(int q=0;q<1000;++q){if(n&&rng()%2){Int i=rng()%n,v=rng()%100;a[i]=v;sums[from_end(n-i)]=v;assert(Int(sums[i])==v);}Int l=rng()%(n+1),r=rng()%(n+1);if(l>r)std::swap(l,r);assert(sums.get(l,r)==std::accumulate(a.begin()+l,a.begin()+r,Int(0)));Int limit=rng()%1000,sm=0,right=l,left=r;while(right<n&&sm+a[right]<=limit)sm+=a[right++];sm=0;while(left>0&&sm+a[left-1]<=limit)sm+=a[--left];auto predicate=[=](Int x){return x<=limit;};assert(sums.max_right(l,predicate)==right&&sums.min_left(r,predicate)==left);}}
    }
}
