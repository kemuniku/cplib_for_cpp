#include <cplib/collections/range_reverse_array.hpp>
#include <cplib/collections/range_reverse_array_monoid.hpp>
#include <cplib/collections/range_reverse_dualsegtree.hpp>
#include <cplib/collections/range_reverse_lazysegtree.hpp>
#include <random>
using namespace cplib;
struct Action{Int a=1,b=0;};
struct Sum{Int sum=0,len=0;bool operator==(const Sum&)const=default;};
Action composition(Action f,Action g){return {f.a*g.a,f.a*g.b+f.b};}
template<class Tree,class S,class Merge,class Map,class Make> void check(Tree tree,std::vector<S> plain,S identity,Merge merge,Map map,Make make){
 std::mt19937 rng(85824);for(int q=0;q<6000;++q){Int n=plain.size(),l=rng()%(n+1),r=rng()%(n+1);if(l>r)std::swap(l,r);int choice=rng()%7;S value=make(rng()%20);
  if(choice==0){tree.reverse(closed_slice(l,r-1));std::reverse(plain.begin()+l,plain.begin()+r);}
  else if(choice==1&&n){Int p=rng()%n;tree[p]=value;plain[p]=value;}
  else if(choice==2){if constexpr(requires{tree.insert(l,value);}){tree.insert(l,value);plain.insert(plain.begin()+l,value);}}
  else if(choice==3){if constexpr(requires{tree.erase(l,r);}){tree.erase(closed_slice(l,r-1));plain.erase(plain.begin()+l,plain.begin()+r);}}
  else if(choice==4&&n){if constexpr(requires{tree.erase(l);}){Int p=rng()%n;tree.erase(p);plain.erase(plain.begin()+p);}}
  else if(choice==5){if constexpr(requires{tree.apply(l,r,Action{});}){Action action{Int(rng()%2),Int(rng()%5)};tree.apply(l,r,action);for(Int i=l;i<r;++i)plain[i]=map(action,plain[i]);}}
  else if(choice==6&&n){if constexpr(requires{tree.apply(l,Action{});}){Int p=rng()%n;Action action{0,3};tree.apply(p,action);plain[p]=map(action,plain[p]);}}
  n=plain.size();assert(tree.len()==n);l=rng()%(n+1);r=rng()%(n+1);if(l>r)std::swap(l,r);
  if constexpr(requires{tree.get(l,r);}){S want=identity,total=identity;for(Int i=l;i<r;++i)want=merge(want,plain[i]);for(const S& x:plain)total=merge(total,x);assert(tree.get(l,r)==want&&tree.fold(closed_slice(l,r-1))==want&&tree.get_all()==total&&tree.fold()==total);}
  if(n){Int p=rng()%n;assert(tree.get(p)==plain[p]);assert(S(tree[from_end(1)])==plain.back());}
  if(q%37==0){assert(tree.toSeq()==plain);std::vector<S> actual;for(S x:items(tree))actual.push_back(x);assert(actual==plain);auto copy=tree;if(n){copy.update(0,value);assert(tree.get(0)==plain[0]);}}
 }
}
template<class Tree> void check_bounds(Tree tree){std::vector<std::string> a{"ab","c","def","g","hi","j"};std::mt19937 rng(8577);for(int q=0;q<1000;++q){Int l=rng()%7,r=rng()%7;if(l>r)std::swap(l,r);tree.reverse(l,r);std::reverse(a.begin()+l,a.begin()+r);Int cap=rng()%15;auto pred=[&](const std::string& x){return Int(x.size())<=cap;};Int right=l,left=r,len=0;while(right<6&&len+Int(a[right].size())<=cap)len+=a[right++].size();len=0;while(left>0&&len+Int(a[left-1].size())<=cap)len+=a[--left].size();assert(tree.max_right(l,pred)==right&&tree.min_left(r,pred)==left);}}
int main(){
 detail::range_reverse_random.seed(982733);
 auto add=std::plus<Int>{};auto intmap=[](Action f,Int x){return f.a*x+f.b;};auto make=[](Int x){return x;};std::vector<Int> initial(25);std::iota(initial.begin(),initial.end(),Int(0));
 check(initRangeReverseArray(initial),initial,Int(0),add,intmap,make);
 check(initRangeReverseArrayMonoid(initial,add,Int(0)),initial,Int(0),add,intmap,make);
 check(initRangeReverseDualSegmentTree(initial,intmap,composition,Action{}),initial,Int(0),add,intmap,make);
 auto sumop=[](Sum x,Sum y){return Sum{x.sum+y.sum,x.len+y.len};};auto summap=[](Action f,Sum x){return Sum{f.a*x.sum+f.b*x.len,x.len};};auto summake=[](Int x){return Sum{x,1};};std::vector<Sum> sums;for(Int x:initial)sums.push_back({x,1});check(initRangeReverseLazySegmentTree(sums,sumop,Sum{},summap,composition,Action{}),sums,Sum{},sumop,summap,summake);
 auto concat=std::plus<std::string>{};auto textmap=[](Action f,std::string x){if(!f.a)std::fill(x.begin(),x.end(),char('a'+f.b%26));else for(char& c:x)c=char('a'+(c-'a'+f.b)%26);return x;};auto textmake=[](Int x){return std::string(1+x%3,char('a'+x));};std::vector<std::string> text{"ab","c","def","g","hi","j"};
 check(initRangeReverseArrayMonoid(text,concat,std::string{}),text,std::string{},concat,textmap,textmake);
 check(initRangeReverseLazySegmentTree(text,concat,std::string{},textmap,composition,Action{}),text,std::string{},concat,textmap,textmake);
 check_bounds(initRangeReverseArrayMonoid(text,concat,std::string{}));check_bounds(initRangeReverseLazySegmentTree(text,concat,std::string{},textmap,composition,Action{}));
 auto empty=initRangeReverseLazySegmentTree(Int(0),sumop,Sum{},summap,composition,Action{});assert(empty.len()==0&&empty.fold()==Sum{});empty.insert(0,{2,1});empty.apply(0,1,{1,5});assert(empty.get(0)==(Sum{7,1}));empty.erase(0);assert(empty.toSeq().empty());
 auto dual=initRangeReverseDualSegmentTree(Int(100),Int(2),intmap,composition,Action{});dual.apply(0,100,{1,3});dual.insert(20,Int(9));assert(dual.get(20)==9&&dual.get(21)==5);
 auto booleans=initRangeReverseArray(std::vector<bool>{true,false,false});booleans.reverse(0,3);assert((booleans.toSeq()==std::vector<bool>{false,false,true}));
 std::vector<Int> big(200000);std::iota(big.begin(),big.end(),Int(0));auto linear=initRangeReverseArray(big);linear.reverse(0,big.size());assert(linear.get(0)==199999&&linear.get(199999)==0);
}
