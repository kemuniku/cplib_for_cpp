#pragma once
#include <cplib/common.hpp>
#include <concepts>
#include <functional>
namespace cplib {
template<class X,class V> struct GoldenSectionResult{X x;V fx;};
template<std::floating_point T,class F> GoldenSectionResult<T,T> golden_section_search(T l,T r,F f,Int iterations=100,bool maximize=false){
 assert(l<=r&&iterations>=0);auto better=[&](T a,T b){return maximize?b<a:a<b;};
 auto interp=[](T a,T b,T t){return a<=0&&0<=b?a*(1-t)+b*t:a+(b-a)*t;};
 T fl=f(l);GoldenSectionResult<T,T> out{l,fl};if(l==r)return out;T fr=f(r);if(better(fr,out.fx))out={r,fr};
 auto eval=[&](T x){T fx=x==l?fl:x==r?fr:T(f(x));if(better(fx,out.fx))out={x,fx};return fx;};
 constexpr T ratio=T(0.6180339887498948482L);T left=l,right=r,x1=interp(l,r,1-ratio),x2=interp(l,r,ratio),f1=eval(x1),f2=x1==x2?f1:eval(x2);
 for(Int i=0;i<iterations;++i){if(!(left<x1&&x1<x2&&x2<right))break;if(better(f2,f1)){left=x1;x1=x2;f1=f2;T next=interp(left,right,ratio);if(!(x1<next&&next<right))break;x2=next;f2=eval(x2);}else{right=x2;x2=x1;f2=f1;T next=interp(left,right,1-ratio);if(!(left<next&&next<x2))break;x1=next;f1=eval(x1);}}return out;
}
template<std::integral I,class F> auto golden_section_search(I l,I r,F f,bool maximize=false){
 static_assert(std::is_signed_v<I>&&sizeof(I)<=8);assert(l<=r);using V=std::invoke_result_t<F,I>;using U=std::make_unsigned_t<I>;
 auto pos=[&](UInt offset){return std::bit_cast<I>(U(U(l)+U(offset)));};UInt width=U(U(r)-U(l)),a=1,b=1;
 while(a-1<width-(b-1)){UInt t=a+b;a=b;b=t;}UInt offset=0;V f1=f(pos(a-1));if(a==b)return GoldenSectionResult<I,V>{l,f1};V f2=f(pos(b-1));bool valid=true;
 while(a!=b){bool right=valid&&(maximize?f1<f2:f2<f1);if(!right){UInt t=b-a;b=a;a=t;f2=f1;valid=true;if(a!=b)f1=f(pos(offset+a-1));}else{offset+=a;UInt t=b-a;b=a;a=t;f1=f2;valid=b-1<=width-offset;if(a!=b&&valid)f2=f(pos(offset+b-1));}}return GoldenSectionResult<I,V>{pos(offset),f1};
}
}
