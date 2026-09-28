#pragma once
#include <cplib/collections/private/range_reverse_treap.hpp>
namespace cplib {
template<class S,class F> using RangeReverseDualSegmentTree=detail::RangeReverseTreap<S,F,false,true>;
template<class S,class F,class Map,class Comp> auto initRangeReverseDualSegmentTree(std::span<const S> values,Map mapping,Comp composition,F id){return RangeReverseDualSegmentTree<S,F>(values,{},S{},mapping,composition,id);}
template<class S,class F,class Map,class Comp> auto initRangeReverseDualSegmentTree(const std::vector<S>& values,Map mapping,Comp composition,F id){return RangeReverseDualSegmentTree<S,F>(values,{},S{},mapping,composition,id);}
template<class S,class F,class Map,class Comp> auto initRangeReverseDualSegmentTree(Int n,S value,Map mapping,Comp composition,F id){assert(n>=0);return initRangeReverseDualSegmentTree(std::vector<S>(n,value),mapping,composition,id);}
template<class V,class Map,class Comp,class F> auto newRangeReverseDualSegWith(const V& values,Map mapping,Comp composition,F id){return initRangeReverseDualSegmentTree(values,mapping,composition,id);}
template<class S,class Map,class Comp,class F> auto newRangeReverseDualSegWith(Int n,S value,Map mapping,Comp composition,F id){return initRangeReverseDualSegmentTree(n,value,mapping,composition,id);}
}
