#pragma once
#include <cplib/common.hpp>
namespace cplib {
// y^h(0) Π(i=1..n)(x y^(h(i)-h(i-1)))。O(log(m+1)log(n+a+b+2))回のop。
template<class T,class Op> T monoidFloorSum(Int n,Int m,Int a,Int b,T x,T y,Op op,T e) {
    assert(n>=0 && m>0 && a>=0 && b>=0);
    assert(n==0 || a<=(std::numeric_limits<Int>::max()-b)/n);
    auto power=[&](T value,Int exponent){
        T result=e;
        while(exponent>0){if(exponent&1)result=op(result,value);exponent>>=1;if(exponent>0)value=op(value,value);}
        return result;
    };
    T prefix=e,suffix=e;
    while(true) {
        prefix=op(prefix,power(y,b/m));x=op(x,power(y,a/m));a%=m;b%=m;
        Int height=(a*n+b)/m;
        if(height==0)return op(op(prefix,power(x,n)),suffix);
        Int last_x=(m*height-b-1)/a+1;
        prefix=op(prefix,x);suffix=op(op(y,power(x,n-last_x)),suffix);
        n=height-1;b=m-b-1;std::swap(a,m);std::swap(x,y);
    }
}
}
