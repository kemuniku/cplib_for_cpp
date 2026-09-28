#include <cplib/math/baser.hpp>
#include <cplib/math/combination_int.hpp>
#include <cplib/math/divisor.hpp>
#include <cplib/math/eratosthenes.hpp>
#include <cplib/math/euler_phi.hpp>
#include <cplib/math/ext_gcd.hpp>
#include <cplib/math/floor_sum.hpp>
#include <cplib/math/inv_gcd.hpp>
#include <cplib/math/mex_naive.hpp>
#include <cplib/math/nearest_equiv.hpp>
#include <cplib/math/osa_k.hpp>
#include <cplib/math/primitive_root.hpp>
#include <iostream>
using namespace cplib;
int main() {
    for(Int a=-100;a<=100;++a) for(Int b=-100;b<=100;++b) {
        auto [x,y]=ext_gcd(a,b);assert(a*x+b*y==std::gcd(a,b));
        if(b>0) {auto [g,i]=inv_gcd(a,b);assert(g==std::gcd(a,b));assert((a*i-g)%b==0);assert(i>=0 && i<b/g);}
    }
    for(Int n=0;n<20;++n) for(Int m=1;m<20;++m) for(Int a=0;a<20;++a) for(Int b=0;b<20;++b) {
        Int want=0;for(Int i=0;i<n;++i) want+=(a*i+b)/m;
        assert(floor_sum(n,m,a,b)==want);
    }
    for(Int n=0;n<10000;++n) {Int x=isqrt(n);assert(x*x<=n && (x+1)*(x+1)>n);}
    assert(isqrt(std::numeric_limits<Int>::max())==3037000499LL);
    for(Int x=-20;x<=20;++x) for(Int l=-20;l<=20;++l) for(Int m=1;m<=20;++m) {
        Int y=nearest_equiv(x,l,m);assert(y>=l && y-m<l && (y-x)%m==0);assert(y==nearest_equiv(x,l,-m));
    }
    auto table=initPrimeFactorTable(20000); auto phi=euler_phi_list(20000);
    for(Int n=1;n<=20000;++n) {
        bool prime=n>=2;
        for(Int d=2;d<=n/d;++d) if(n%d==0){prime=false;break;}
        assert(isprime(n)==prime);auto factors=primefactor(n);assert(factors==table.primefactor(n));
        Int product=1;for(Int p:factors) product*=p;assert(product==n);
        assert(euler_phi(n)==phi[n]);
    }
    assert(isprime(UInt(18446744073709551557ULL)));
    assert(!isprime(UInt(18446744073709551615ULL)));
    assert(!isprime(Int(341550071728321LL)));
    for(Int n:{1000000007LL*1000000009LL,9223372036854775807LL,9999999967LL}) {
        auto f=primefactor(n);__int128 product=1;for(Int p:f){assert(isprime(p));product*=p;}assert(product==n);
    }
    auto sieve=initEratosthenes(2000000);assert(sieve.byte_size()==2000000/30+1);assert(sieve.count_primes()==148933);
    Int count=0,prev=0;for(Int p:sieve){assert(p>prev && isprime(p));prev=p;++count;}assert(count==sieve.count_primes());
    for(Int low=0;low<100;++low) for(Int high=low;high<150;++high) {
        auto seg=initSegmentedEratosthenes(low,high);Int want=0;
        for(Int n=low;n<=high;++n){assert(seg.is_prime(n)==isprime(n));want+=isprime(n);}
        assert(seg.count_primes()==want);auto primes=get_primes(low,high+1);assert(Int(primes.size())==want);
    }
    for(Int low:{999900LL,999999900LL,1000000000000LL}) {
        auto seg=initSegmentedEratosthenes(low,low+1000);
        for(Int n=low;n<=low+1000;++n) assert(seg.is_prime(n)==isprime(n));
    }
    for(Int p:{2,3,5,17,97,998244353}) {Int g=primitive_root(p);for(auto [q,k]:primefactor_tuple(p-1)) assert(powmod(g,(p-1)/q,p)!=1);}
    assert(ncr_int(5,2)==10);assert(ncr_int(100,50,1000000)==1000000);assert(ncr_int(-1,0)==0);
    assert((baser(10,2)==std::vector<Int>{0,1,0,1}));assert(baser(0,2).empty());
    assert(mex_naive(std::vector<Int>{0,2,3,0,-1})==1);
    assert((divisor(36)==std::vector<Int>{1,2,3,4,6,9,12,18,36}));
    assert(divisor(1000000007LL*1000000009LL).size()==4);
    assert(powmod(-2,3,5)==2);assert(powmod(1,0,1)==0);
    std::cout<<"math tests passed\n";
}
