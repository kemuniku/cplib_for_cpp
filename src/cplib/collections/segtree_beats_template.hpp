#pragma once
#include <cplib/collections/segtree_beats.hpp>
#include <cplib/utils/constants.hpp>
namespace cplib {
template<class T> struct S_rch {T max{},max2{},min{},min2{},sum{};Int sz=0,n_min=0,n_max=0;bool fail=false;};
template<class T> struct F_rch {T lb,ub,add;};
template<class T> S_rch<T> init_S(T value,T inf,Int size=1){return {value,-inf,value,inf,value*T(size),size,size,size,false};}
template<class T> std::ostream& operator<<(std::ostream& out,const S_rch<T>& x){return out<<"(max: "<<x.max<<", max2: "<<x.max2<<", min: "<<x.min<<", min2: "<<x.min2<<", sum: "<<x.sum<<", sz: "<<x.sz<<", n_min: "<<x.n_min<<", n_max: "<<x.n_max<<", fail: "<<(x.fail?"true":"false")<<")";}
// 区間chmin/chmax/addとsum/min/max。元と同じ第2極値による失敗判定を使う。
template<class T> class RangeChminChmaxRangeSumMaxMin {
 static T second_lowest(T a,T b,T c,T d){if(a==c)return std::min(b,d);if(b<=c)return b;if(d<=a)return d;return std::max(a,c);}
 static S_rch<T> merge(S_rch<T> l,S_rch<T> r){return {std::max(l.max,r.max),-second_lowest(-l.max,-l.max2,-r.max,-r.max2),std::min(l.min,r.min),second_lowest(l.min,l.min2,r.min,r.min2),l.sum+r.sum,l.sz+r.sz,l.n_min*Int(l.min<=r.min)+r.n_min*Int(r.min<=l.min),l.n_max*Int(l.max>=r.max)+r.n_max*Int(r.max>=l.max),true};}
public:
 SegmentTreeBeats<S_rch<T>,F_rch<T>> seg;T inf,zero;
 RangeChminChmaxRangeSumMaxMin(std::span<const T> values,T infinity,T zero):inf(infinity),zero(zero){
  auto identity=S_rch<T>{-inf,-inf,inf,inf,zero,0,0,0,false};
  auto composition=[](F_rch<T> f,F_rch<T> g){return F_rch<T>{std::max(std::min(g.lb+g.add,f.ub),f.lb)-g.add,std::min(std::max(g.ub+g.add,f.lb),f.ub)-g.add,f.add+g.add};};
  auto mapping=[infinity,identity](F_rch<T> f,S_rch<T> x){x.fail=false;if(!x.sz)return identity;if(x.min==x.max||f.lb==f.ub||f.lb>=x.max||f.ub<=x.min)return init_S(std::min(std::max(x.min,f.lb),f.ub)+f.add,infinity,x.sz);
   if(x.min2==x.max){x.max2=x.min=std::max(x.min,f.lb)+f.add;x.min2=x.max=std::min(x.max,f.ub)+f.add;x.sum=x.min*T(x.n_min)+x.max*T(x.n_max);return x;}
   if(f.lb<x.min2&&f.ub>x.max2){T next_min=std::max(x.min,f.lb),next_max=std::min(x.max,f.ub);x.sum+=(next_min-x.min)*T(x.n_min)-(x.max-next_max)*T(x.n_max)+f.add*T(x.sz);x.min=next_min+f.add;x.max=next_max+f.add;x.min2+=f.add;x.max2+=f.add;return x;}x.fail=true;return x;};
  std::vector<S_rch<T>> initial;initial.reserve(values.size());for(T x:values)initial.push_back(init_S(x,inf));initSegmentTreeBeatsInPlace(seg,initial,merge,identity,mapping,composition,F_rch<T>{-inf,inf,zero});
 }
 Int len()const{return seg.len();}void update(Int p,T value){seg.update(p,init_S(value,inf));}S_rch<T> get(Int p){return seg.get(p);}S_rch<T> get(Int l,Int r){return seg.get(l,r);}
 struct Reference:S_rch<T>{RangeChminChmaxRangeSumMaxMin* owner;Int index;Reference(RangeChminChmaxRangeSumMaxMin* tree,Int i):S_rch<T>(tree->get(i)),owner(tree),index(i){}Reference& operator=(T value){owner->update(index,value);static_cast<S_rch<T>&>(*this)=init_S(value,owner->inf);return *this;}};
 Reference operator[](Int p){return {this,p};}Reference operator[](BackwardsIndex p){return (*this)[len()-p.value];}
 template<class L,class R> S_rch<T> operator[](ClosedSlice<L,R> range){return seg.get(range);}
 void chmin(Int l,Int r,T value){seg.apply(l,r,{-inf,value,zero});}void chmax(Int l,Int r,T value){seg.apply(l,r,{value,inf,zero});}void add(Int l,Int r,T value){seg.apply(l,r,{-inf,inf,value});}
 template<class L,class R> void chmin(ClosedSlice<L,R> range,T value){seg.apply(range,{-inf,value,zero});}
 template<class L,class R> void chmax(ClosedSlice<L,R> range,T value){seg.apply(range,{value,inf,zero});}
 template<class L,class R> void add(ClosedSlice<L,R> range,T value){seg.apply(range,{-inf,inf,value});}
 std::string str(){return seg.str();}
};
template<class T> auto initRangeChminChmaxRangeSumMaxMin(std::span<const T> values,T inf,T zero){return RangeChminChmaxRangeSumMaxMin<T>(values,inf,zero);}
template<class T> auto initRangeChminChmaxRangeSumMaxMin(const std::vector<T>& values,T inf,T zero){return initRangeChminChmaxRangeSumMaxMin(std::span<const T>(values),inf,zero);}
template<class T> auto initRangeChminChmaxRangeSumMaxMin(std::span<const T> values){T inf;if constexpr(std::is_floating_point_v<T>)inf=T(1e100);else if constexpr(sizeof(T)>=8)inf=T(INF64);else inf=T(INF32);return initRangeChminChmaxRangeSumMaxMin(values,inf,T(0));}
template<class T> auto initRangeChminChmaxRangeSumMaxMin(const std::vector<T>& values){return initRangeChminChmaxRangeSumMaxMin(std::span<const T>(values));}
template<class T> void update(RangeChminChmaxRangeSumMaxMin<T>& tree,Int p,T value){tree.update(p,value);}
template<class T> Int len(const RangeChminChmaxRangeSumMaxMin<T>& tree){return tree.len();}
template<class T> std::string to_string(RangeChminChmaxRangeSumMaxMin<T>& tree){return tree.str();}
template<class T,class L,class R> void chmin(RangeChminChmaxRangeSumMaxMin<T>& tree,ClosedSlice<L,R> range,T value){tree.chmin(range,value);}
template<class T,class L,class R> void chmax(RangeChminChmaxRangeSumMaxMin<T>& tree,ClosedSlice<L,R> range,T value){tree.chmax(range,value);}
template<class T,class L,class R> void add(RangeChminChmaxRangeSumMaxMin<T>& tree,ClosedSlice<L,R> range,T value){tree.add(range,value);}
}
