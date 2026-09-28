#include <cplib/convolution/convolution_old.hpp>
#include <random>
using namespace cplib;
template<class T,bool Old> void check_mod(){std::mt19937 rng(31);for(Int n:{0,1,2,17,60,61,64,127,128,129,200,511,700})for(Int m:{0,1,3,60,61,64,90,129,300}){std::vector<T> a(n),b(m);for(auto& x:a)x=T(Int(rng())-Int(rng()));for(auto& x:b)x=T(Int(rng())-Int(rng()));auto want=convolution_naive(a,b);assert(convolution(a,b)==want);}

    if constexpr(Old)for(Int n=1;n<=4096;n*=2){std::vector<T> f(n);for(auto& x:f)x=T(rng());auto original=f;ntt(f);intt(f);assert(f==original);}
}
int main(){check_mod<StaticMontgomeryModint<998244353>,true>();check_mod<StaticBarrettModint<998244353>,true>();check_mod<StaticMontgomeryModint<754974721>,true>();std::mt19937 rng(32);for(int trial=0;trial<250;++trial){Int n=rng()%200,m=rng()%200;std::vector<Int> a(n),b(m);for(auto& x:a)x=Int(rng()%200001)-100000;for(auto& x:b)x=Int(rng()%200001)-100000;assert(convolution_ll(a,b)==convolution_naive(a,b));}}
