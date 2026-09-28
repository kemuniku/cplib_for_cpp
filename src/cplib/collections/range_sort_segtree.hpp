#pragma once
#include <cplib/collections/private/range_sort_impl.hpp>
namespace cplib {
template<class T> class RangeSortSegmentTree:public detail::RangeSortEngine<T,true>{
 using Base=detail::RangeSortEngine<T,true>;
public:
 template<class Op> RangeSortSegmentTree(std::span<const Int> keys,std::span<const T> values,Int limit,Op op,T e):Base(keys,values,limit,op,e){}
 using Base::len;using Base::get;using Base::sort;
 T operator[](Int i)const{return this->value(i);}T operator[](BackwardsIndex i)const{return this->value(len()-i.value);}
 void set(Int i,T v){this->update(i,std::move(v));}void set(BackwardsIndex i,T v){set(len()-i.value,std::move(v));}
 template<class L,class R> T get(ClosedSlice<L,R> s)const{return get(resolve_index(len(),s.a),resolve_index(len(),s.b)+1);}
 template<class L,class R> T operator[](ClosedSlice<L,R> s)const{return get(s);}
 template<class L,class R> void sort(ClosedSlice<L,R> s,SortOrder order=Ascending){sort(resolve_index(len(),s.a),resolve_index(len(),s.b)+1,order);}
};
template<class Values,class Op,class T> auto initRangeSortSegmentTree(std::span<const Int> keys,const Values& values,Int limit,Op op,T e){return RangeSortSegmentTree<T>(keys,values,limit,op,e);}
}
