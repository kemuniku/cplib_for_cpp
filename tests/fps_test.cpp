#include <cplib/fps/formal_power_series.hpp>
#include <cplib/fps/bmbm.hpp>
#include <cplib/fps/taylor_shift.hpp>
#include <cplib/convolution/relaxed_convolution.hpp>
#include <cplib/fps/product_tree.hpp>
#include <cplib/fps/polynomial_interpolation.hpp>
#include <cplib/fps/product_of_polynomial_sequence.hpp>
#include <cplib/fps/shift_of_sampling_points.hpp>
#include <cplib/fps/composition.hpp>
#include <cplib/fps/power_projection.hpp>
#include <cplib/fps/fps.hpp>
#include <cplib/math/many_factorials.hpp>
#include <random>
using namespace cplib;
template<class T> void check_fps(){std::mt19937 rng(35);for(Int n:{0,1,2,17,63,64,65,128,129,300,1024}){std::vector<T> a(n);for(auto& x:a)x=rng();if(n){a[0]=1;assert(inv(a)==invRelaxed(a));assert(log(a)==logRelaxed(a));assert(prefix(a*inv(a),n)==prefix(std::vector<T>{1},n));assert(derivative(integral(a))==a);a[0]=0;assert(exp(a)==expRelaxed(a));assert(log(exp(a))==a);a[0]=rng();auto sq=prefix(a*a,n);auto root=sqrt(sq);assert(root&&prefix((*root)*(*root),n)==sq);assert(root==sqrtRelaxed(sq));for(Int k:{0,1,2,5})assert(pow(a,k)==powRelaxed(a,k));if(n>3){a[0]=a[1]=0;assert(pow(a,3)==powRelaxed(a,3));sq=prefix(a*a,n);assert(sqrt(sq)==sqrtRelaxed(sq));}}
        assert(prefix(a,-1).empty());assert(coefficient(a,-1)==T(0)&&coefficient(a,n)==T(0));for(Int size:{0,1,5,100}){std::vector<T> b(size);for(auto& x:b)x=rng();if(!b.empty()){b.back()=1;auto [q,r]=divmod(a,b);assert(normalized(q*b+r)==normalized(a));assert(r.empty()||r.size()<normalized(b).size());assert(div(a,b)==q&&mod(a,b)==r);}auto c=a+b;c-=b;assert(normalized(c)==normalized(a));c=a;c+=b;assert(c==a+b);c=a;c+=2;assert(c==a+2&&c==2+a);c=a;c*=T(3);assert(c==a*3&&c==3*a);c/=T(3);assert(c==a);}
        T shift=rng();auto shifted=taylorShift(a,shift);for(int k=0;k<10;++k){T x=rng();assert(eval(shifted,x)==eval(a,x+shift));}assert(taylorShift(shifted,-shift)==a);
    }
    assert((exp(std::vector<T>{},5)==std::vector<T>{1,0,0,0,0}));assert(!sqrt(std::vector<T>{0,1}));assert(sqrt(std::vector<T>{},5)==std::optional<std::vector<T>>(std::vector<T>(5)));
    for(int trial=0;trial<150;++trial){Int d=1+rng()%20;std::vector<T> a(300),c(d);for(auto& x:c)x=rng();for(Int i=0;i<d;++i)a[i]=rng();for(Int i=d;i<Int(a.size());++i)for(Int j=0;j<d;++j)a[i]+=c[j]*a[i-j-1];auto first=std::vector<T>(a.begin(),a.begin()+2*d);auto inferred=berlekampMassey(first);for(Int i=inferred.size();i<Int(a.size());++i){T x=0;for(Int j=0;j<Int(inferred.size());++j)x+=inferred[j]*a[i-j-1];assert(x==a[i]);}for(Int k:{0,1,27,100,299}){assert(linearRecurrenceKth(std::vector<T>(a.begin(),a.begin()+d),c,k)==a[k]);assert(bmbm(first,k)==a[k]);}}
    assert(berlekampMassey(std::vector<T>{}).empty());assert(berlekampMassey(std::vector<T>(20)).empty());assert((berlekampMassey(std::vector<T>{1,0,0,0})==std::vector<T>{0}));assert(bmbm(std::vector<T>{},100)==T(0));
    {std::vector<T> p={1},q={1,-T(1)};assert(bostanMori(p,q,1000000000000000000LL)==T(1));}
}
void extra_fps_test();
void sparse_test();
void factorial_test();
int main(){extra_fps_test();sparse_test();factorial_test();check_fps<StaticMontgomeryModint<998244353>>();check_fps<StaticBarrettModint<998244353>>();check_fps<StaticBarrettModint<1000000007>>();using C=StaticBarrettModint<35>;assert((integral(std::vector<C>{1,2,3})==std::vector<C>{0,1,1,1}));}

template<class T> void extra_fps_check(){
    std::mt19937 rng(36);
    for(Int n:{0,1,2,17,61,64,65,128,129,300}){std::vector<T> f(n),xs(n);for(auto& x:f)x=rng();for(Int i=0;i<n;++i)xs[i]=T(i);auto tree=initPolynomialProductTree(xs);auto values=multipointEvaluation(f,xs);assert(tree.evaluate(f)==values);for(Int i=0;i<n;++i)assert(values[i]==eval(f,xs[i]));assert(polynomialInterpolation(xs,values)==f);for(int q=0;q<4;++q){T shift=rng();Int m=rng()%350;auto moved=shiftOfSamplingPoints(values,shift,m);for(Int i=0;i<m;++i)assert(moved[i]==eval(f,shift+T(i)));}if(n){f.resize(n*3);for(auto& x:f)x=rng();auto got=multipointEvaluation(f,xs);assert(tree.evaluate(f)==got);for(Int i=0;i<n;++i)assert(got[i]==eval(f,xs[i]));}}
    for(int trial=0;trial<150;++trial){std::vector<std::vector<T>> factors(rng()%60);std::vector<T> want={T(1)};for(auto& p:factors){p.resize(rng()%25);for(auto& x:p)x=rng();want=convolution_naive(want,p);}assert(productOfPolynomialSequence(factors)==want);}
    {std::vector<std::vector<T>> factors(200,std::vector<T>{1,1});auto got=productOfPolynomialSequence(factors);std::vector<T> want={1};for(auto& f:factors)want=convolution_naive(want,f);assert(got==want);factors.push_back({T(0)});assert(productOfPolynomialSequence(factors)==std::vector<T>(201));factors={{T(3)},{T(7)}};assert((productOfPolynomialSequence(factors)==std::vector<T>{21}));}
    for(Int n:{1,2,5,17,32,65,129}){std::vector<T> outer(n),inner(n),g(n);for(auto& x:outer)x=rng();for(auto& x:inner)x=rng();for(auto& x:g)x=rng();inner[0]=0;std::vector<T> want(n),p(n);p[0]=1;for(Int i=0;i<n;++i){for(Int j=0;j<n;++j)want[j]+=outer[i]*p[j];p=prefix(convolution_naive(p,inner),n);}assert(compose(outer,inner,n)==want);if(n>=2){inner[1]=1;auto inverse=compositionalInverse(inner);assert(compose(inner,inverse)==prefix(std::vector<T>{0,1},n));assert(compose(inverse,inner)==prefix(std::vector<T>{0,1},n));}for(Int zero:{0,1}){outer[0]=T(zero);Int m=n+3;auto fixed=powerProjection(outer,g,m),diagonal=powerProjectionDiagonal(outer,g,m);std::vector<T> power(m+1);power[0]=1;for(Int i=0;i<=m;++i){auto product=convolution_naive(power,g);assert(fixed[i]==coefficient(product,n-1));assert(diagonal[i]==coefficient(product,i));power=prefix(convolution_naive(power,outer),m+1);}}}
}
void extra_fps_test(){extra_fps_check<StaticBarrettModint<998244353>>();extra_fps_check<StaticMontgomeryModint<998244353>>();extra_fps_check<StaticBarrettModint<1000000007>>();using M=StaticBarrettModint<17>;std::vector<M> f={1,5,3,2},ys(17);for(Int i=0;i<17;++i)ys[i]=eval(f,M(i));for(Int t=0;t<17;++t){auto out=shiftOfSamplingPoints(ys,M(t),70);for(Int i=0;i<70;++i)assert(out[i]==eval(f,M(t+i)));}}

void sparse_test(){
    using T=StaticMontgomeryModint<998244353>;std::mt19937 rng(37);
    auto f=sfps<T>([](auto x){return 1+2*x-(x^3)+(x^5);});auto want=initSparseFPS<T>({{0,1},{1,2},{3,-T(1)},{5,1}});assert(f==want);assert(SFPS<T>([](auto x){return 2*(x^3)-3*(x^3)+1+x;})==initSparseFPS<T>({{0,1},{1,1},{3,-T(1)}}));assert(sfps<T>([](auto){return 5;})==initSparseFPS<T>({{0,5}}));
    for(int trial=0;trial<300;++trial){Int n=1+rng()%160;SparseFPS<T> terms;for(int i=0;i<15;++i)terms.push_back({Int(rng()%(n+20)),T(rng())});auto s=initSparseFPS<T>(terms);auto dense=toDense(s,n);for(Int i=-1;i<=n+20;++i){T want=0;for(auto term:terms)if(term.degree==i)want+=term.coefficient;assert(coefficient(s,i)==want);}auto plus=s+2;assert(toDense(plus,n)==dense+2);plus+=3;assert(plus==5+s);assert(isZero(s)==s.empty());assert(degree(s)==(s.empty()?-1:s.back().degree));
        std::vector<T> b(n);for(auto& x:b)x=rng();assert(mulPrefix(b,s,n)==prefix(convolution_naive(b,dense),n));auto c=b;c*=s;assert(c==mulPrefix(b,s,n));auto full=b*s;assert(full==convolution_naive(b,toDense(s,degree(s)+1)));assert(s*b==full);
        // 定数項を1にして除算・対数を比較。
        SparseFPS<T> unit={{0,T(1)}};for(auto term:s)if(term.degree>0)unit.push_back(term);auto d=toDense(unit,n);assert(inv(unit,n)==inv(d,n));assert(log(unit,n)==log(d,n));assert(divPrefix(b,unit,n)==prefix(b*inv(d,n),n));auto divided=b;divided/=unit;assert(divided==b/unit);for(Int k:{0,1,2,7})assert(pow(unit,k,n)==pow(d,k,n));assert(sqrt(unit,n)==sqrt(d,n));
        auto zero=unit;zero.erase(zero.begin());assert(exp(zero,n)==exp(toDense(zero,n),n));for(Int k:{0,1,2,4})assert(pow(s,k,n)==pow(dense,k,n));assert(sqrt(s,n)==sqrt(toDense(s,degree(s)+1),n));
    }
    for(int trial=0;trial<50;++trial){Int dim=1+rng()%5,n=20+rng()%300;SparseFPS<T> unit={{0,1}},zero;for(Int i=1;i<=dim;++i){T x=rng();unit.push_back({i,x});zero.push_back({i,x});}auto iv=inv(unit,n),lg=log(unit,n),ex=exp(zero,n);for(Int k:{Int(0),Int(1),n/2,n-1}){assert(invCoefficient(unit,k)==iv[k]);assert(logCoefficient(unit,k)==lg[k]);assert(expCoefficient(zero,k)==ex[k]);}for(Int exponent:{0,1,2,9,100}){auto power=pow(unit,exponent,n);for(Int k:{Int(0),Int(1),n/2,n-1})assert(powCoefficient(unit,exponent,k)==power[k]);}for(auto& term:unit)term.degree+=2;auto power=pow(unit,3,n);for(Int k:{Int(0),Int(6),n/2,n-1})if(k<n)assert(powCoefficient(unit,3,k)==power[k]);}
    assert((initSparseFPS<T>({{3,2},{3,-T(2)},{3,5},{0,0}})==initSparseFPS<T>({{3,5}})));
    {auto e=initSparseFPS<T>({{1,1}});T fact=1;Int n=10000;for(Int i=1;i<=n;++i)fact*=T(i);assert(expCoefficient(e,n)==fact.inv());auto u=initSparseFPS<T>({{0,1},{1,-T(1)}});assert(invCoefficient(u,1000000000000LL)==T(1));}
}

void factorial_test(){
    using T=StaticMontgomeryModint<998244353>;std::mt19937 rng(38);std::vector<T> facts(200001);facts[0]=1;for(Int i=1;i<Int(facts.size());++i)facts[i]=facts[i-1]*T(i);for(Int width:{1,2,32,256,1024}){auto table=initLargeFactorial<T>(200000,width);for(Int n=0;n<100;++n)assert(table.fact(n)==facts[n]);for(int i=0;i<1000;++i){Int n=rng()%facts.size();assert(table.fact(n)==facts[n]);}assert(table.fact(998244353)==T(0));}for(Int limit:{0,1,100,1024,1025,4000,200000}){std::vector<Int> queries(500);for(auto& n:queries)n=rng()%(limit+1);queries.push_back(998244353);queries.push_back(1000000000);auto result=manyFactorials<T>(queries);for(std::size_t i=0;i<result.size();++i)assert(result[i]==(queries[i]<998244353?facts[queries[i]]:T(0)));}using S=StaticBarrettModint<17>;auto all=initLargeFactorial<S>();S fact=1;for(Int i=0;i<17;++i){if(i)fact*=S(i);assert(all.fact(i)==fact);}assert(all.fact(17)==S(0));
    auto big=initLargeFactorial<T>();assert(big.fact(998244352)==T(-1)&&big.fact(998244351)==T(1));assert(manyFactorials<T>(std::vector<Int>{998244352,998244351,0})==std::vector<T>({-T(1),T(1),T(1)}));
}
