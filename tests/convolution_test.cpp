#include <cplib/convolution/convolution.hpp>
#include <cplib/convolution/ntt.hpp>
#include <cplib/convolution/bitwise_and_convolution.hpp>
#include <cplib/convolution/xor_convolution.hpp>
#include <cplib/convolution/gcd_convolution.hpp>
#include <cplib/convolution/lcm_convolution.hpp>
#include <cplib/convolution/min_plus_convolution.hpp>
#include <cplib/convolution/semi_relaxed_convolution.hpp>
#include <cplib/convolution/relaxed_convolution.hpp>
#include <random>
using namespace cplib;
template<class T,bool Old> void check_mod(){std::mt19937 rng(31);for(Int n:{0,1,2,17,60,61,64,127,128,129,200,511,700})for(Int m:{0,1,3,60,61,64,90,129,300}){std::vector<T> a(n),b(m);for(auto& x:a)x=T(Int(rng())-Int(rng()));for(auto& x:b)x=T(Int(rng())-Int(rng()));auto want=convolution_naive(a,b);assert(convolution(a,b)==want);Int size=std::bit_ceil(UInt(std::max({n,m,Int(1)})));auto cyc=convolutionCyclicPowerOfTwo(a,b,size);std::vector<T> c(size);for(std::size_t i=0;i<want.size();++i)c[i%size]+=want[i];assert(cyc==c);}
    if constexpr(Old)for(Int n=1;n<=4096;n*=2){std::vector<T> f(n);for(auto& x:f)x=T(rng());auto original=f;ntt(f);intt(f);assert(f==original);}
}
void extra_test();
void relaxed_test();
int main(){
    extra_test();
    relaxed_test();
    check_mod<StaticMontgomeryModint<998244353>,true>();check_mod<StaticBarrettModint<998244353>,true>();check_mod<StaticMontgomeryModint<754974721>,true>();check_mod<StaticBarrettModint<1000000007>,false>();check_mod<StaticBarrettModint<1000000000>,false>();check_mod<StaticBarrettModint<2147483647>,false>();check_mod<StaticBarrettModint<1>,false>();
    using D=DynamicMontgomeryModint<32>;for(int mod:{998244353,167772161,469762049}){D::setMod(mod);check_mod<D,true>();}using B=DynamicBarrettModint<32>;for(int mod:{998244353,1000000007}){B::setMod(mod);if(mod==998244353)check_mod<B,true>();else check_mod<B,false>();}
    std::mt19937 rng(32);for(int trial=0;trial<250;++trial){Int n=rng()%200,m=rng()%200;std::vector<Int> a(n),b(m);for(auto& x:a)x=Int(rng()%200001)-100000;for(auto& x:b)x=Int(rng()%200001)-100000;auto want=convolution_naive(a,b);assert(convolution_ll(a,b)==want);auto mod=convolution<1000000007>(a,b);for(std::size_t i=0;i<mod.size();++i)assert(mod[i]==detail::convolution_floor_mod(want[i],1000000007));}
    {std::vector<Int> a={std::numeric_limits<Int>::min(),std::numeric_limits<Int>::max()},b={1};assert(convolution_ll(a,b)==a);}
    for(Int n=0;n<100;++n){std::vector<Int> a(n),b(n),gcd(n),lcm(n);for(auto& x:a)x=Int(rng()%200)-100;for(auto& x:b)x=Int(rng()%200)-100;for(Int i=0;i<n;++i)for(Int j=0;j<n;++j){gcd[std::gcd(i,j)]+=a[i]*b[j];Int l=std::lcm(i,j);if(l<n)lcm[l]+=a[i]*b[j];}assert(gcdConvolution(a,b)==gcd&&lcmConvolution(a,b)==lcm);}
    for(Int n=1;n<=256;n*=2){std::vector<Int> a(n),b(n),x(n),an(n);for(auto& v:a)v=Int(rng()%200)-100;for(auto& v:b)v=Int(rng()%200)-100;for(Int i=0;i<n;++i)for(Int j=0;j<n;++j){x[i^j]+=a[i]*b[j];an[i&j]+=a[i]*b[j];}assert(xorConvolution(a,b)==x&&bitwiseAndConvolution(a,b)==an);using M=StaticMontgomeryModint<998244353>;std::vector<M> am(a.begin(),a.end()),bm(b.begin(),b.end()),xm(x.begin(),x.end());assert(xorConvolution(am,bm)==xm);}
    assert(bitwiseAndConvolution(std::vector<Int>{},std::vector<Int>{}).empty());
}

void extra_test(){
    std::mt19937 rng(33);
    auto naive=[](const auto& a,const auto& b){std::vector<Int> out;if(a.empty()||b.empty())return out;out.assign(a.size()+b.size()-1,std::numeric_limits<Int>::max());for(std::size_t i=0;i<a.size();++i)for(std::size_t j=0;j<b.size();++j)out[i+j]=std::min(out[i+j],a[i]+b[j]);return out;};
    for(int trial=0;trial<500;++trial){Int n=rng()%350,m=rng()%350;auto shape=[&](Int size){std::vector<Int> out(size),diff(size);for(auto& x:diff)x=Int(rng()%400)-200;std::sort(diff.begin(),diff.end());Int sum=Int(rng()%400)-200;for(Int i=0;i<size;++i){out[i]=sum;sum+=diff[i];}return out;};auto a=shape(n),b=shape(m);std::vector<Int> other(m);for(auto& x:other)x=Int(rng()%10000)-5000;assert(minPlusConvolutionConvexConvex(a,b)==naive(a,b));auto want=naive(a,other);assert(minPlusConvolutionConvexArbitraryMonotoneMinima(a,other)==want&&minPlusConvolutionConvexArbitrarySmawk(a,other)==want);for(auto& x:a)x=-x;for(auto& x:b)x=-x;assert(minPlusConvolutionConcaveConcave(a,b)==naive(a,b));assert(minPlusConvolutionConcaveArbitrary(a,other)==naive(a,other));}
    using M=StaticMontgomeryModint<998244353>;
    for(Int n:{0,1,17,61,300,1024}){std::vector<M> fixed(n),online(1100);for(auto& x:fixed)x=rng();for(auto& x:online)x=rng();auto s=initSemiRelaxedConvolution(fixed);auto want=convolution_naive(fixed,online);want.resize(online.size());for(Int i=0;i<Int(online.size());++i){M result=i%3==0?s.add(online[i]):i%3==1?s.append(online[i]):s.get(online[i]);assert(result==want[i]&&s.len()==i+1);if(i%127==0)assert(s.coefficients()==std::vector<M>(want.begin(),want.begin()+i+1));}assert(s.coefficients()==want);}
}

template<class T> void relaxed_check(){
    std::mt19937 rng(34);
    for(Int n:{0,1,2,17,64,65,128,200,1024}){std::vector<T> a(n),b(n);for(auto& x:a)x=rng();for(auto& x:b)x=rng();auto state=initRelaxedConvolution<T>(n);auto want=convolution_naive(a,b);want.resize(n);for(Int i=0;i<n;++i){T known=0;for(Int j=1;j<i;++j)known+=a[j]*b[i-j];assert(state.pendingCoefficient()==known);auto result=i%3==0?state.add(a[i],b[i]):i%3==1?state.append(a[i],b[i]):state.get(a[i],b[i]);assert(result==want[i]&&state.len()==i+1);}assert(state.coefficients()==want);if(!n)continue;
        a[0]=T(1+rng()%(T::umod()-1));auto inverse=invRelaxed(a);auto identity=convolution_naive(a,inverse);identity.resize(n);assert(identity[0]==T(1));for(Int i=1;i<n;++i)assert(identity[i]==T(0));auto invstate=initRelaxedInv<T>(n);for(Int i=0;i<n;++i)assert(invstate.get(a[i])==inverse[i]);assert(invstate.coefficients()==inverse);
        a[0]=0;auto exponential=expRelaxed(a);std::vector<T> expected(n);expected[0]=1;for(Int k=1;k<n;++k){T sum=0;for(Int i=1;i<=k;++i)sum+=T(i)*a[i]*expected[k-i];expected[k]=sum/T(k);}assert(expected==exponential);assert(logRelaxed(exponential)==a);auto ex=initRelaxedExp<T>(n);auto log=initRelaxedLog<T>(n);for(Int i=0;i<n;++i){assert(ex.append(a[i])==exponential[i]);assert(log.get(exponential[i])==a[i]);}assert(ex.coefficients()==exponential&&log.coefficients()==a);
        a[0]=T(1+rng()%(T::umod()-1));auto sq=convolution_naive(a,a);sq.resize(n);auto root=sqrtRelaxed(sq);assert(root);auto square=convolution_naive(*root,*root);square.resize(n);assert(square==sq);auto rt=initRelaxedSqrt<T>(n);for(Int i=0;i<n;++i)assert(rt.append(sq[i])==std::optional<T>((*root)[i]));assert(rt.coefficients()==root);
        if(n>2){a[0]=a[1]=0;sq=convolution_naive(a,a);sq.resize(n);root=sqrtRelaxed(sq);assert(root);square=convolution_naive(*root,*root);square.resize(n);assert(square==sq);}for(Int k:{0,1,2,5}){std::vector<T> power(n);power[0]=1;for(Int e=0;e<k;++e){power=convolution_naive(power,a);power.resize(n);}assert(powRelaxed(a,k)==power);auto ps=initRelaxedPow<T>(n,k);for(Int i=0;i<n;++i)assert(ps.append(a[i])==power[i]);assert(ps.coefficients()==power);}
    }
    assert(invRelaxed(std::vector<T>{T(1)},0).empty());assert((invRelaxed(std::vector<T>{T(2)},5)==std::vector<T>{T(2).inv(),0,0,0,0}));assert((expRelaxed(std::vector<T>{0},5)==std::vector<T>{1,0,0,0,0}));assert((logRelaxed(std::vector<T>{1},5)==std::vector<T>(5)));assert(sqrtRelaxed(std::vector<T>{},5)==std::optional<std::vector<T>>(std::vector<T>(5)));assert(!sqrtRelaxed(std::vector<T>{0,1},5));T nonresidue=2;while(nonresidue.pow((T::umod()-1)/2)==T(1))nonresidue+=T(1);auto rt=initRelaxedSqrt<T>(3);assert(!rt.add(nonresidue)&&!rt.add(T(1))&&!rt.add(T(0))&&rt.len()==3&&!rt.coefficients());
}
void relaxed_test(){relaxed_check<StaticMontgomeryModint<998244353>>();relaxed_check<StaticBarrettModint<1000000007>>();}
