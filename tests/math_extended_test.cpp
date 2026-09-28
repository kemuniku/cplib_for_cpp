#include <cplib/math/stern_brocot_tree.hpp>
#include <cplib/math/int128.hpp>
#include <cplib/math/float128.hpp>
#include <cplib/math/fractions.hpp>
#include <cplib/math/generalized_floor_sum.hpp>
#include <cplib/math/sqrt_heuristic_for_floor_sum.hpp>
#include <cplib/math/xor_basis.hpp>
#include <set>
#include <random>
using namespace cplib;
void stern_test();
int main(){
    stern_test();
    for(std::string s:{"0","1","-1","170141183460469231731687303715884105727","-170141183460469231731687303715884105728"}) assert(to_string(parseInt128(s))==s);
    assert(pow(Int128(3),Int128(10))==59049);assert(pow(Int128(3),Int128(10),Int128(7))==4);
    assert(to_float(to_float128(1)/to_float128(2))==0.5);
    using F=Fraction<Int>;
    assert(F(2,4)==F(1,2));assert(F(1,2)+F(1,3)==F(5,6));assert(2/F(1,3)==F(6));
    assert(F(1,0)>F(100));assert(F(-1,0)<F(-100));assert((F(1,0)+F(-1,0)).isNaN());
    assert(inv(F(-2,3))==F(-3,2));assert(to_string(F(6,8,false))=="3/4");
    for(Int n=0;n<12;++n)for(Int m=1;m<12;++m)for(Int a=0;a<12;++a)for(Int b=0;b<12;++b){
        auto tab=generalizedFloorSumTable<Int>(n,m,a,b,2,2);
        for(Int p=0;p<=2;++p)for(Int q=0;q<=2;++q){Int want=0;for(Int i=0;i<n;++i){Int v=1;for(Int j=0;j<p;++j)v*=i;for(Int j=0;j<q;++j)v*=(a*i+b)/m;want+=v;}assert(tab[p][q]==want);}
        auto cat=[](const std::string& x,const std::string& y){return x+y;};
        auto got=monoidFloorSum(n,m,a,b,std::string("x"),std::string("y"),cat,std::string());
        std::string want(b/m,'y');for(Int i=1;i<=n;++i)want+="x"+std::string((a*i+b)/m-(a*(i-1)+b)/m,'y');assert(got==want);
    }
    for(Int n=0;n<50;++n)for(Int m=1;m<40;++m)for(Int a=0;a<m;++a){
        Int b=m/2;std::vector<Int> got(n,-1);
        for(auto [x,y,dx,dy,count]:sqrt_heuristic_for_floor_sum(a,b,n,m))for(Int i=0;i<count;++i){assert(x+i*dx<n);assert(got[x+i*dx]==-1);got[x+i*dx]=y+i*dy;}
        for(Int i=0;i<n;++i)assert(got[i]==(a*i+b)%m);
    }
    std::mt19937 rng(42);
    for(int trial=0;trial<200;++trial){
        std::vector<Int> input;for(int i=0;i<8;++i)input.push_back(rng()%256);
        auto basis=initXorBasis(input);auto incremental=initXorBasis();for(Int x:input)incremental.incl(x);assert(basis.basis==incremental.basis);
        std::set<Int> values{0};for(Int x:input){auto old=values;for(Int y:old)values.insert(x^y);}
        std::vector<Int> v(values.begin(),values.end());
        for(Int x=0;x<300;++x){assert(basis.can_make(x)==values.count(x));auto it=values.lower_bound(x);assert(basis.lt(x)==(it==values.begin()?-1:*std::prev(it)));auto jt=values.upper_bound(x);assert(basis.le(x)==*std::prev(jt));}
        for(Int k=0;k<Int(v.size());++k){assert(basis.kth_smallest(k)==v[k]);assert(basis.index(v[k])==k);}
        Int x=rng()%256;auto sorted=v;std::sort(sorted.begin(),sorted.end(),[&](Int a,Int b){return (a^x)<(b^x);});
        assert(basis.xor_min(x)==sorted[0]);for(Int k=0;k<Int(v.size());++k)assert(basis.xor_kth(x,k)==sorted[k]);
        assert(basis.kth_smallest(v.size())==-1);
    }
}
void stern_test(){
    for(cplib::Int a=1;a<=40;++a)for(cplib::Int b=1;b<=40;++b){using namespace cplib;auto x=to_SBTNode(a,b);Int g=std::gcd(a,b);assert(x.num()==a/g&&x.den()==b/g);assert(decode_path(encode_path(x))==x);std::vector<SBTNode<Int>> chain{sbt_root()};for(auto [c,d]:encode_path(x))for(Int i=0;i<d;++i)chain.push_back(c=='L'?move_left(chain.back(),Int(1)):move_right(chain.back(),Int(1)));assert(Int(chain.size())==x.depth+1);for(Int i=0;i<=x.depth;++i)assert(ancestor(x,i)==chain[i]);assert(!ancestor(x,x.depth+1));for(Int m=1;m<=10;++m){Fraction<Int> lower(0,1),upper(1,0),target(a,b);for(Int p=1;p<=m;++p)for(Int q=1;q<=m;++q){Fraction<Int> f(p,q);if(f<target&&lower<f)lower=f;if(target<f&&f<upper)upper=f;}assert(Fraction<Int>(max_less_with_den_at_most(x,m))==lower);assert(Fraction<Int>(min_greater_with_den_at_most(x,m))==upper);auto bounds=get_bounds([&](SBTNode<Int> n){return Fraction<Int>(n)<target;},m);auto [lo,hi]=get_range_fraction(bounds);assert(lo==lower);assert(hi==((a/g<=m&&b/g<=m)?target:upper));}}
    using namespace cplib;std::vector<SBTNode<Int>> v;for(Int a=1;a<12;++a)for(Int b=1;b<12;++b)if(std::gcd(a,b)==1)v.push_back(to_SBTNode(a,b));for(auto a:v)for(auto b:v){auto c=LCA(a,b);auto p=encode_path(a),q=encode_path(b);std::string x,y;for(auto [ch,d]:p)x+=std::string(d,ch);for(auto [ch,d]:q)y+=std::string(d,ch);Int len=0;while(len<Int(std::min(x.size(),y.size()))&&x[len]==y[len])++len;assert(c.depth==len);assert(ancestor(a,len)==c&&ancestor(b,len)==c);}auto aux=initAuxiliaryWeightedTree(v);assert(aux.graph.edge_count()==Int(aux.v.size())-1);for(auto e:aux.graph.edge_info){auto a=aux.v[e.src],b=aux.v[e.dst];assert(LCA(a,b)==a&&e.cost==b.depth-a.depth);}
}
