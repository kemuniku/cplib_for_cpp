#pragma once
#include <cplib/collections/private/lazysegtree_base.hpp>
namespace cplib {
using detail::lazysegtree_base::apply;
using detail::lazysegtree_base::get;
using detail::lazysegtree_base::get_all;
using detail::lazysegtree_base::len;
using detail::lazysegtree_base::max_right;
using detail::lazysegtree_base::min_left;
using detail::lazysegtree_base::to_string;
using detail::lazysegtree_base::update;

// 単位作用を必要としないhasLazy版。演算を静的に特殊化する。
template<class S,class F,class Op,class Map,class Comp> using LazySegmentTree=detail::lazysegtree_base::LazySegmentTree<S,F,Op,Map,Comp,true>;
template<class S,class F,class V,class Op,class Map,class Comp> auto initLazySegmentTree(const V& v,Op op,S e,Map map,Comp comp,F id){static_assert(std::is_empty_v<Op>&&std::is_empty_v<Map>&&std::is_empty_v<Comp>);return LazySegmentTree<S,F,Op,Map,Comp>(v,op,e,map,comp,id);}
template<class V,class Op,class S,class Map,class Comp,class F> auto newLazySegWith(const V& v,Op op,S e,Map map,Comp comp,F id){return cplib::initLazySegmentTree(v,op,e,map,comp,id);}
template<auto Op,auto Map,auto Comp,class V,class S,class F> auto initLazySegmentTree(const V& v,S e,F id){return cplib::initLazySegmentTree(v,[](S a,S b){return Op(a,b);},e,[](F f,S x){return Map(f,x);},[](F f,F g){return Comp(f,g);},id);}
}
