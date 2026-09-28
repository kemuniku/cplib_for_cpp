#pragma once
#include <cplib/utils/backwards_index.hpp>
namespace cplib {
class RangeLinearAddRangeMin {
    struct Point {Int x=0,y=0;};struct Node {Point left,right;Int slope=0,intercept=0;};
    Int length_;std::vector<Node> nodes_;
    static Point shifted(Point p,Int s,Int t){p.y+=s*p.x+t;return p;}
    static __int128 cross(Point a,Point b,Point c,Point d){return (__int128(b.x)-a.x)*(__int128(d.y)-c.y)-(__int128(b.y)-a.y)*(__int128(d.x)-c.x);}
    // 左右の下側凸包をたどり、共通接線を更新。O(log N)。
    void pull(Int k,Int border){Int l=k*2,r=k*2+1,ls=nodes_[l].slope,lc=nodes_[l].intercept,rs=nodes_[r].slope,rc=nodes_[r].intercept;
        for(;;){auto a=shifted(nodes_[l].left,ls,lc),b=shifted(nodes_[l].right,ls,lc),c=shifted(nodes_[r].left,rs,rc),d=shifted(nodes_[r].right,rs,rc);bool lLeaf=a.x==b.x,rLeaf=c.x==d.x;if(lLeaf&&rLeaf){nodes_[k].left=a;nodes_[k].right=c;return;}bool descendLeft=false;Int child;
            if(!lLeaf&&cross(a,b,a,c)<0){descendLeft=true;child=l*2;}else if(!rLeaf&&cross(b,c,b,d)<0)child=r*2+1;else if(lLeaf)child=r*2;else if(rLeaf){descendLeft=true;child=l*2+1;}else{auto c1=cross(a,b,c,d),c2=cross(a,b,c,b);descendLeft=c1==0&&c2==0?c.x<border:__int128(c.x-border)*c1+__int128(d.x-c.x)*c2<0;child=descendLeft?l*2+1:r*2;}
            if(descendLeft){l=child;ls+=nodes_[l].slope;lc+=nodes_[l].intercept;}else{r=child;rs+=nodes_[r].slope;rc+=nodes_[r].intercept;}}
    }
    void build(const std::vector<Int>& v,Int k,Int l,Int r){if(r-l==1){nodes_[k].left=nodes_[k].right={l,v[l]};return;}Int m=(l+r)/2;build(v,k*2,l,m);build(v,k*2+1,m,r);pull(k,m);}
    void push(Int k){for(Int child=k*2;child<=k*2+1;++child){nodes_[child].slope+=nodes_[k].slope;nodes_[child].intercept+=nodes_[k].intercept;}nodes_[k].slope=nodes_[k].intercept=0;}
    void addImpl(Int k,Int l,Int r,Int ql,Int qr,Int b,Int c){if(ql<=l&&r<=qr){nodes_[k].slope+=b;nodes_[k].intercept+=c;return;}push(k);Int m=(l+r)/2;if(ql<m)addImpl(k*2,l,m,ql,qr,b,c);if(m<qr)addImpl(k*2+1,m,r,ql,qr,b,c);pull(k,m);}
    Int subtreeMin(Int k,Int s,Int t)const{for(;;){s+=nodes_[k].slope;t+=nodes_[k].intercept;auto a=shifted(nodes_[k].left,s,t),b=shifted(nodes_[k].right,s,t);if(a.x==b.x)return a.y;k=a.y<b.y?k*2:k*2+1;}}
    Int prodImpl(Int k,Int l,Int r,Int ql,Int qr,Int s,Int t)const{if(ql<=l&&r<=qr)return subtreeMin(k,s,t);Int m=(l+r)/2;s+=nodes_[k].slope;t+=nodes_[k].intercept;Int result=std::numeric_limits<Int>::max();if(ql<m)result=prodImpl(k*2,l,m,ql,qr,s,t);if(m<qr)result=std::min(result,prodImpl(k*2+1,m,r,ql,qr,s,t));return result;}
public:
    // 構築 O(N)、一次式区間加算・区間最小値 O(log² N)。
    explicit RangeLinearAddRangeMin(const std::vector<Int>& v):length_(v.size()),nodes_(4*v.size()){if(length_)build(v,1,0,length_);}
    Int len()const{return length_;}void add(Int l,Int r,Int b,Int c){assert(0<=l&&l<=r&&r<=length_);if(l<r)addImpl(1,0,length_,l,r,b,c);}
    Int prod(Int l,Int r)const{assert(0<=l&&l<=r&&r<=length_);return l==r?std::numeric_limits<Int>::max():prodImpl(1,0,length_,l,r,0,0);}
    template<class L,class R> void add(ClosedSlice<L,R> s,Int b,Int c){add(resolve_index(len(),s.a),resolve_index(len(),s.b)+1,b,c);}template<class L,class R> Int prod(ClosedSlice<L,R> s)const{return prod(resolve_index(len(),s.a),resolve_index(len(),s.b)+1);}template<class L,class R> Int operator[](ClosedSlice<L,R> s)const{return prod(s);}
    Int operator[](Int i)const{assert(0<=i&&i<len());return prod(i,i+1);}Int operator[](BackwardsIndex i)const{return (*this)[len()-i.value];}
};
inline auto initRangeLinearAddRangeMin(const std::vector<Int>& v){return RangeLinearAddRangeMin(v);}
}
