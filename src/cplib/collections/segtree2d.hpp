#pragma once
#include <cplib/collections/segtree.hpp>
namespace cplib {
// 密な2次元セグメント木。構築O(HW)、点更新・矩形取得O(log H log W)。
template<class T> class SegmentTree2D {
    T identity;std::function<T(T,T)> merge;Int lastnode=1,H=0,W=0;
public:
    std::vector<SegmentTree<T>> arr;
    template<class Op> SegmentTree2D(const std::vector<std::vector<T>>& v,Op op,T e):identity(e),merge(op),H(v.size()),W(H?Int(v[0].size()):0){while(lastnode<H)lastnode*=2;arr.reserve(lastnode*2);for(Int i=0;i<lastnode*2;++i)arr.emplace_back(W,op,e);for(Int i=0;i<H;++i)arr[i+lastnode]=initSegmentTree(v[i],op,e);for(Int i=lastnode-1;i>0;--i){std::vector<T> tmp(W);for(Int j=0;j<W;++j)tmp[j]=merge(arr[2*i][j],arr[2*i+1][j]);arr[i]=initSegmentTree(tmp,op,e);}}
    T get(Int il,Int ir,Int jl,Int jr)const{assert(0<=il&&il<=ir&&ir<=H&&0<=jl&&jl<=jr&&jr<=W);il+=lastnode;ir+=lastnode;T left=identity,right=identity;while(il<ir){if(il&1)left=merge(left,arr[il++].get(jl,jr));if(ir&1)right=merge(arr[--ir].get(jl,jr),right);il>>=1;ir>>=1;}return merge(left,right);}
    void update(Int i,Int j,const T& value){assert(0<=i&&i<H&&0<=j&&j<W);i+=lastnode;arr[i].update(j,value);while(i>1){i>>=1;arr[i].update(j,merge(arr[2*i][j],arr[2*i+1][j]));}}
};
template<class T,class Op> auto initSegmentTree2D(const std::vector<std::vector<T>>& v,Op op,T e){return SegmentTree2D<T>(v,op,e);}
template<class T> T get(const SegmentTree2D<T>& s,Int il,Int ir,Int jl,Int jr){return s.get(il,ir,jl,jr);}
template<class T> void update(SegmentTree2D<T>& s,Int i,Int j,const T& v){s.update(i,j,v);}
}
