#pragma once
#include <cplib/collections/private/range_reverse_treap.hpp>
namespace cplib {
template<class S,class F> using RangeReverseLazySegmentTree=detail::RangeReverseTreap<S,F,true,true>;
template<class S,class F,class V,class Op,class Map,class Comp> auto initRangeReverseLazySegmentTree(const V& values,Op op,S e,Map mapping,Comp composition,F id){if constexpr(std::integral<V>){assert(values>=0);return RangeReverseLazySegmentTree<S,F>(std::vector<S>(values,e),op,e,mapping,composition,id);}else return RangeReverseLazySegmentTree<S,F>(values,op,e,mapping,composition,id);}
template<class V,class Op,class S,class Map,class Comp,class F> auto newRangeReverseLazySegWith(const V& values,Op op,S e,Map mapping,Comp composition,F id){return initRangeReverseLazySegmentTree(values,op,e,mapping,composition,id);}
}
