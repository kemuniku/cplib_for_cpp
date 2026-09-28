#pragma once
#include <cplib/common.hpp>
#include <cmath>
#include <concepts>
namespace cplib {
// 境界を二分探索する。O(log |ok-ng|)回の判定。
template<class Predicate> Int meguru_bisect(Int ok,Int ng,Predicate is_ok){
    while(static_cast<__int128>(ok)-ng>1 || static_cast<__int128>(ng)-ok>1){Int mid=static_cast<Int>((static_cast<__int128>(ok)+ng)/2);if(is_ok(mid))ok=mid;else ng=mid;}return ok;
}
// 絶対・相対誤差のいずれかがeps以下になるまで二分探索する。
template<std::floating_point T,class Predicate> T meguru_bisect(T ok,T ng,Predicate is_ok,T eps=T(1e-10)){
    while(std::abs(ok-ng)>eps && std::abs(ok-ng)/std::max(std::abs(ok),std::abs(ng))>eps){T mid=(ok+ng)/2;if(is_ok(mid))ok=mid;else ng=mid;}return ok;
}
}
