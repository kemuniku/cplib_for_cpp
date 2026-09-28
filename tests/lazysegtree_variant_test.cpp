#include CPLIB_LAZY_HEADER
#include CPLIB_DUAL_HEADER
#include <cplib/collections/lazysegtree_template.hpp>
#include <random>
using namespace cplib;
void lazy_test(){
    std::mt19937 rng(13);
    using S=RangeSum<Int>;using F=RangeAffine<Int>;
    auto merge=[](S a,S b){return S{a.sum+b.sum,a.len+b.len};};auto map=[](F f,S x){return S{f.a*x.sum+f.b*x.len,x.len};};auto comp=[](F f,F g){return F{f.a*g.a,f.a*g.b+f.b};};auto point=[](F f,Int x){return f.a*x+f.b;};
    for(Int n:{0,1,2,5,16,17,65,257}){
        std::vector<Int> a(n);std::vector<S> nodes(n);for(Int i=0;i<n;++i){a[i]=rng()%10;nodes[i]={a[i],1};}
        auto s=initLazySegmentTree(nodes,merge,S{},map,comp,F{1,0});auto st=cplib::initLazySegmentTree(nodes,merge,S{},map,comp,F{1,0});auto d=initDualSegmentTree(a,point,comp,F{1,0});auto dt=cplib::initDualSegmentTree(a,point,comp,F{1,0});
        for(int q=0;q<2000;++q){Int l=rng()%(n+1),r=rng()%(n+1);if(l>r)std::swap(l,r);F f{Int(rng()%2),Int(rng()%5)};
            if(q%3){s.apply(l,r,f);st.apply(closed_slice(l,r-1),f);d.apply(l,r,f);dt.apply(closed_slice(l,r-1),f);for(Int i=l;i<r;++i)a[i]=point(f,a[i]);}
            else if(n){Int p=rng()%n,x=rng()%15;s[p]=S{x,1};st[p]=S{x,1};d[p]=x;dt[p]=x;a[p]=x;}
            Int sum=std::accumulate(a.begin()+l,a.begin()+r,Int(0));assert(s.get(l,r).sum==sum&&st.get(l,r).sum==sum);assert(s.get_all()==st.get_all());assert(s.get_all().sum==std::accumulate(a.begin(),a.end(),Int(0)));
            if(n){Int p=rng()%n;assert(s.get(p).sum==a[p]&&st.get(p).sum==a[p]&&d.get(p)==a[p]&&dt.get(p)==a[p]);}
            Int cap=rng()%80,rr=l,ll=r,sm=0;while(rr<n&&sm+a[rr]<=cap)sm+=a[rr++];sm=0;while(ll>0&&sm+a[ll-1]<=cap)sm+=a[--ll];auto pred=[&](S x){return x.sum<=cap;};assert(s.max_right(l,pred)==rr&&st.max_right(l,pred)==rr);assert(s.min_left(r,pred)==ll&&st.min_left(r,pred)==ll);
        }
        assert(d.toSeq()==a&&dt.toSeq()==a);Int i=0;for(Int x:dt)assert(x==a[i++]);
        if(n){s[from_end(1)]=S{8,1};st[from_end(1)]=S{8,1};d[from_end(1)]=8;dt[from_end(1)]=8;assert(S(s[from_end(1)]).sum==8&&S(st[from_end(1)]).sum==8&&Int(d[from_end(1)])==8&&Int(dt[from_end(1)])==8);}
    }
    std::vector<Int> a(39);for(auto& x:a)x=Int(rng()%30)-15;
    auto amin=initRangeAddRangeMin(a);auto amax=initRangeAddRangeMax(a);auto smin=initRangeAssignRangeMin(a);auto smax=initRangeAssignRangeMax(a);
    auto ami=initRangeAddRangeMinIndex(a);auto ama=initRangeAddRangeMaxIndex(a);auto smi=initRangeAssignRangeMinIndex(a);auto sma=initRangeAssignRangeMaxIndex(a);
    auto asum=initRangeAddRangeSum(a);auto ssum=initRangeAssignRangeSum(a);auto aff=initRangeAffineRangeSum(a);auto b=a,c=a;
    for(int q=0;q<1000;++q){Int l=rng()%40,r=rng()%40;if(l>r)std::swap(l,r);Int x=Int(rng()%11)-5;F f{Int(rng()%3)-1,x};
        amin.apply(l,r,x);amax.apply(l,r,x);ami.apply(l,r,x);ama.apply(l,r,x);asum.apply(l,r,x);smin.apply(l,r,x);smax.apply(l,r,x);smi.apply(l,r,x);sma.apply(l,r,x);ssum.apply(l,r,x);aff.apply(l,r,f);for(Int i=l;i<r;++i){a[i]+=x;b[i]=x;c[i]=point(f,c[i]);}
        l=rng()%40;r=rng()%40;if(l>r)std::swap(l,r);
        auto check=[&](auto& lo,auto& hi,auto& li,auto& ri,auto& sum,const auto& vec){assert(sum.get(l,r).sum==std::accumulate(vec.begin()+l,vec.begin()+r,Int(0)));if(l==r){assert(li.get(l,r).index==-1&&ri.get(l,r).index==-1);assert(lo.get(l,r)==std::numeric_limits<Int>::max()&&hi.get(l,r)==std::numeric_limits<Int>::lowest());}else{Int mn=std::min_element(vec.begin()+l,vec.begin()+r)-vec.begin(),mx=std::max_element(vec.begin()+l,vec.begin()+r)-vec.begin();assert(lo.get(l,r)==vec[mn]&&hi.get(l,r)==vec[mx]);assert(li.get(l,r).value==vec[mn]&&li.get(l,r).index==mn&&li.get(l,r).left==l);assert(ri.get(l,r).value==vec[mx]&&ri.get(l,r).index==mx);}};
        check(amin,amax,ami,ama,asum,a);check(smin,smax,smi,sma,ssum,b);assert(aff.get(l,r).sum==std::accumulate(c.begin()+l,c.begin()+r,Int(0)));
    }

#ifdef CPLIB_STATIC_LAZY
    {auto assign=[](Int f,Int){return f;};auto overwrite=[](Int f,Int){return f;};auto d=cplib::initDualSegmentTree(Int(3),Int(7),assign,overwrite,Int(0));assert(d.get(1)==7);d.apply(0,3,0);assert(d.toSeq()==std::vector<Int>(3,0));}
#endif
}
int main(){lazy_test();}
