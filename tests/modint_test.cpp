#include <cplib/modint/modint.hpp>
#include <cplib/modint/exp_modint.hpp>
#include <cplib/collections/lazysegtree_template.hpp>
#include <cplib/math/combination_prefix_sum.hpp>
#include <cplib/math/modfast.hpp>
#include <random>
#include <unordered_set>
using namespace cplib;
template<class M> void test(){
    std::mt19937_64 rng(21);Int m=M::mod();
    for(int i=0;i<50000;++i){Int a=std::bit_cast<Int>(rng()),b=std::bit_cast<Int>(rng());M x=a,y=b;Int ar=a%m,br=b%m;if(ar<0)ar+=m;if(br<0)br+=m;assert(x.val()==ar&&y.val()==br);assert((x+y).val()==(ar+br)%m);assert((x-y).val()==(ar-br+m)%m);assert((x*y).val()==Int(__int128(ar)*br%m));assert((-x).val()==(m-ar)%m);assert((Int(1)+x).val()==(ar+1)%m);assert((Int(1)-x).val()==(1-ar+m)%m);assert((UInt(2)*x).val()==ar*2%m);if(br&&std::gcd(br,m)==1){assert((y*y.inv()).val()==1);assert((x/y*y)==x);assert((Int(1)/y)==y.inv());}UInt u=rng();assert(M(u).val()==Int(u%UInt(m)));assert((x+M(m))==x);assert(std::hash<M>{}(x+M(m))==std::hash<M>{}(x));Int e=rng()%50,w=1%m;for(Int j=0;j<e;++j)w=Int(__int128(w)*ar%m);assert(x.pow(e).val()==w);}
    std::unordered_set<M> values;values.insert(M(0));values.insert(M(m));values.insert(M(1));assert(values.size()==std::size_t(m>1?2:1));
}
void modfast_combination_test();
int main(){
    modfast_combination_test();
    static_assert(sizeof(modint998244353_montgomery)==4&&sizeof(modint998244353_barrett)==4);
    static_assert((modint998244353_montgomery(5)*7).val()==35);
    static_assert(divide_constant<2>(modint998244353_barrett(10)).val()==5);
    test<modint998244353_montgomery>();test<modint1000000007_montgomery>();test<modint998244353_barrett>();test<modint1000000007_barrett>();test<StaticMontgomeryModint<1>>();test<StaticBarrettModint<1>>();test<StaticMontgomeryModint<9>>();test<StaticBarrettModint<10>>();test<StaticBarrettModint<2147483647>>();
    for(Int m:{1,3,9,1000000007}){modint_montgomery::setMod(m);test<modint_montgomery>();}
    for(Int m:{1,2,10,1000000007,2147483647}){modint_barrett::setMod(m);test<modint_barrett>();}
    DynamicMontgomeryModint<2>::setMod(17);assert(DynamicMontgomeryModint<2>::mod()==17&&modint_montgomery::mod()==1000000007);
    DynamicBarrettModint<2>::setMod(19);assert(DynamicBarrettModint<2>::mod()==19&&modint_barrett::mod()==2147483647);
    using M=modint998244353_montgomery;assert(estimate_rational(M(3)/7,20)=="3/7");assert(estimate_rational(M(-4)/9,20)=="-4/9");assert(estimate_rational(M(0),5)=="0/1");
    auto s=initRangeAffineRangeSum(std::vector<M>{1,2,3});s.apply(0,3,RangeAffine<M>{2,1});assert(s.get_all().sum==M(15));
    // 元のexpmodintは法1の状態を常に0にする。そのphi鎖の規則をそのまま照合。
    for(Int a=0;a<100;++a)for(Int b=0;b<100;++b){auto x=tomodint<17>(a),y=tomodint<17>(b);assert((x+y).val()==(a+b)%17);assert((x*y).val()==a*b%17);auto p=pow(x,y);Int want=1;for(Int i=0;i<b;++i)want=want*a%17;assert(p.val()==want);}
}
void modfast_combination_test(){
    using M=modint998244353_montgomery;auto c=initCombination<M>(300);std::vector<std::vector<M>> pascal(151);for(Int n=0;n<=150;++n){pascal[n].resize(n+1);pascal[n][0]=pascal[n][n]=1;for(Int r=1;r<n;++r)pascal[n][r]=pascal[n-1][r-1]+pascal[n-1][r];for(Int r=0;r<=n;++r){assert(c.ncr(n,r)==pascal[n][r]);M product=1;for(Int i=0;i<r;++i)product*=n-i;assert(c.npr(n,r)==product);}}assert(c.nhr(0,0)==M(1));assert(c.ncr(-1,0)==M(0));
    std::mt19937 rng(22);std::vector<std::pair<Int,Int>> prefixes;std::vector<std::tuple<Int,Int,Int>> ranges;for(int i=0;i<10000;++i){Int n=rng()%151,y=Int(rng()%180)-10,l=Int(rng()%180)-10,r=Int(rng()%180)-10;prefixes.emplace_back(n,y);ranges.emplace_back(n,l,r);}auto p=combinationPrefixSum<M>(prefixes),s=combinationRangeSum<M>(ranges);for(Int i=0;i<Int(p.size());++i){auto [n,y]=prefixes[i];M want=0;for(Int j=0;j<=std::min(n,y);++j)want+=pascal[n][j];assert(p[i]==want);auto [nn,l,r]=ranges[i];want=0;for(Int j=std::max<Int>(0,l);j<std::min(nn+1,r);++j)want+=pascal[nn][j];assert(s[i]==want);}
    using D=DynamicBarrettModint<42>;
    for(Int modulus:{2,3,5,17,61,67,101,1009,65537,1000003}){D::setMod(modulus);auto t=initModFast<D>();D x=1;for(Int e=0;e<modulus-1;++e){if(modulus<10000||e%17==0){assert(t.log(x)==e);assert(t.powRoot(e)==x);assert(t.inv(x)*x==D(1));Int k=Int(rng()%2001)-1000;D want=k>=0?x.pow(k):x.inv().pow(-k);assert(t.pow(x,k)==want);}x*=t.root();}assert(t.pow(0,0)==D(1)&&t.pow(0,3)==D(0));assert(t.powRoot(-1)*t.root()==D(1));}
}
