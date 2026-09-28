#include <cplib/collections/segtree.hpp>
#include <cplib/collections/lazysegtree.hpp>
#include <cplib/collections/lazysegtree_template.hpp>
#include <cplib/collections/dualsegtree.hpp>
#include <cplib/collections/compressed_segtree.hpp>
#include <cplib/collections/compressed_lazysegtree.hpp>
#include <cplib/collections/compressed_fenwick2d.hpp>
#include <cplib/collections/segtree2d.hpp>
#include <cplib/collections/compressed_segtree2d.hpp>
#include <cplib/collections/dynamic_segtree.hpp>
#include <cplib/collections/dynamic_lazysegtree.hpp>
#include <random>
#include <map>
#include <set>
using namespace cplib;
void compressed_lazy_test();
void compressed_fenwick_test();
void dynamic_lazy_test();
void dynamic_persistent_test();
void compressed_test();
void lazy_test();
int main(){lazy_test();compressed_test();compressed_lazy_test();compressed_fenwick_test();dynamic_persistent_test();dynamic_lazy_test();}
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
void compressed_test(){
    std::mt19937_64 rng(14);
    for(Int n:{0,1,63,64,100,2047,2048,10000}){
        std::vector<Int> coords(n);for(auto& x:coords)x=std::bit_cast<Int>(rng());auto sorted=coords;std::sort(sorted.begin(),sorted.end());detail::sortCompressedCoordinates(coords);assert(coords==sorted);auto slots=detail::initCompressedCoordinateIndex<Int>(coords);for(Int i=0;i<n;++i)assert(detail::findCompressedCoordinate<Int>(coords,slots,coords[i])==i);for(int i=0;i<100;++i){Int x=std::bit_cast<Int>(rng());auto it=std::lower_bound(coords.begin(),coords.end(),x);Int want=it!=coords.end()&&*it==x?it-coords.begin():-1;assert(detail::findCompressedCoordinate<Int>(coords,slots,x)==want);}
    }
    {std::vector<UInt> a(10000);for(auto& x:a)x=rng();auto b=a;std::sort(b.begin(),b.end());detail::sortCompressedCoordinates(a);assert(a==b);}
    for(int trial=0;trial<100;++trial){std::vector<Int> input(100);for(auto& x:input)x=Int(rng()%300)-150;auto op=[](std::string a,const std::string& b){return a+b;};auto init=[](Int x){return std::to_string(x)+",";};auto a=initCompressedSegmentTree(input,op,std::string{},init);auto b=newCompressedSegWith(input,op,std::string{},init);std::map<Int,std::string> naive;for(Int x:input)naive[x]=init(x);assert(a.len()==Int(naive.size()));
        for(int q=0;q<100;++q){Int x=input[rng()%input.size()];auto value=std::to_string(q)+"!";naive[x]=value;a[x]=value;b[x]=value;Int l=Int(rng()%350)-175,r=Int(rng()%350)-175;if(l>r)std::swap(l,r);std::string want,closed,all;for(auto& [key,val]:naive){all+=val;if(l<=key&&key<r)want+=val;if(l<=key&&key<=r)closed+=val;}assert(a.get(l,r)==want&&b.get(l,r)==want);assert(a.get(closed_slice(l,r))==closed&&b.get(closed_slice(l,r))==closed);assert(a.get_all()==all&&b.get_all()==all);assert(std::string(a[x])==value);assert(std::string(a[1000]).empty());}
    }
    {std::vector<std::string> coords={"z","a","m","a"};auto s=initCompressedSegmentTree(coords,[](Int a,Int b){return a+b;},Int(0));s["m"]=4;s["z"]=8;assert(s.get("b","z")==4&&s.get(closed_slice(std::string("b"),std::string("z")))==12);}
}
void compressed_lazy_test(){
    using S=RangeSum<Int>;using F=RangeAffine<Int>;auto merge=[](S a,S b){return S{a.sum+b.sum,a.len+b.len};};auto map=[](F f,S x){return S{f.a*x.sum+f.b*x.len,x.len};};auto comp=[](F f,F g){return F{f.a*g.a,f.a*g.b+f.b};};
    std::mt19937 rng(16);std::vector<Int> coords;for(Int i=-100;i<=100;i+=5)coords.push_back(i);
    auto initial=[](Int x){return S{x,1};};auto p=initCompressedLazySegmentTree(coords,merge,S{},map,comp,F{1,0},initial);auto ps=newCompressedLazySegWith(coords,merge,S{},map,comp,F{1,0},initial);
    auto interval=[](Int l,Int r){return S{0,r-l};};auto s=initCompressedLazySegmentTree(coords,merge,S{},map,comp,F{1,0},interval);auto ss=newCompressedLazySegWith(coords,merge,S{},map,comp,F{1,0},interval);std::vector<Int> a=coords,b(200,0);
    assert(p.len()==41&&s.len()==40);assert(S(p[102]).len==0);
    for(int q=0;q<2000;++q){Int l=Int(rng()%241)-120,r=Int(rng()%241)-120;if(l>r)std::swap(l,r);F f{Int(rng()%3)-1,Int(rng()%9)-4};p.apply(l,r,f);ps.apply(l,r,f);for(Int i=0;i<41;++i)if(l<=coords[i]&&coords[i]<r)a[i]=f.a*a[i]+f.b;
        Int li=rng()%41,ri=rng()%41;if(li>ri)std::swap(li,ri);s.apply(coords[li],coords[ri],f);ss.apply(coords[li],coords[ri],f);for(Int i=coords[li]+100;i<coords[ri]+100;++i)b[i]=f.a*b[i]+f.b;
        if(q%5==0){Int i=rng()%41;S x{Int(rng()%11),1};a[i]=x.sum;p[coords[i]]=x;ps[coords[i]]=x;i=rng()%40;x={Int(rng()%11)*5,5};s[coords[i]]=x;ss[coords[i]]=x;std::fill(b.begin()+coords[i]+100,b.begin()+coords[i+1]+100,x.sum/5);}
        l=Int(rng()%241)-120;r=Int(rng()%241)-120;if(l>r)std::swap(l,r);Int sum=0,cl=0;for(Int i=0;i<41;++i){if(l<=coords[i]&&coords[i]<r)sum+=a[i];if(l<=coords[i]&&coords[i]<=r)cl+=a[i];}assert(p.get(l,r).sum==sum&&ps.get(l,r).sum==sum);assert(p.get(closed_slice(l,r)).sum==cl&&ps.get(closed_slice(l,r)).sum==cl);
        li=rng()%41;ri=rng()%41;if(li>ri)std::swap(li,ri);sum=std::accumulate(b.begin()+coords[li]+100,b.begin()+coords[ri]+100,Int(0));assert(s.get(coords[li],coords[ri]).sum==sum&&ss.get(coords[li],coords[ri]).sum==sum);assert(s.get(closed_slice(coords[li],coords[ri]-1)).sum==sum&&ss.get(closed_slice(coords[li],coords[ri]-1)).sum==sum);
        assert(s.get_all()==ss.get_all()&&p.get_all()==ps.get_all());assert(s.get_all().sum==std::accumulate(b.begin(),b.end(),Int(0)));assert(p.get_all().sum==std::accumulate(a.begin(),a.end(),Int(0)));
    }
}
void compressed_fenwick_test(){
    std::mt19937 rng(17);
    for(int trial=0;trial<80;++trial){std::vector<std::tuple<Int,Int,Int>> ps;std::vector<std::pair<Int,Int>> coords;std::map<std::pair<Int,Int>,Int> ref;for(int i=0;i<200;++i){Int x=Int(rng()%50)-25,y=Int(rng()%50)-25,w=Int(rng()%10)-5;ps.emplace_back(x,y,w);coords.emplace_back(x,y);ref[{x,y}]+=w;}auto op=[](Int a,Int b){return a+b;};auto s=initCompressedSegmentTree2D(ps,op,Int(0));auto ss=newCompressedSeg2DWith(ps,op,Int(0));auto zero=initCompressedSegmentTree2D(coords,op,Int(0));for(auto [p,w]:ref)zero.update(p.first,p.second,w);auto t=initCompressedFenwick2D(ps);auto z=initCompressedFenwick2D(coords);for(auto [p,w]:ref)z.add(p.first,p.second,w);assert(t.len()==Int(ref.size())&&z.len()==t.len());
        for(int q=0;q<500;++q){auto p=coords[rng()%coords.size()];Int w=Int(rng()%20)-10;if(q%2){t.add(p.first,p.second,w);z.add(p.first,p.second,w);ref[p]+=w;}else{t[p]=w;z[p]=w;ref[p]=w;}s[p]=ref[p];ss[p]=ref[p];zero[p]=ref[p];assert(Int(s[p])==ref[p]&&Int(ss[p])==ref[p]);assert(Int(t[p])==ref[p]);Int xl=Int(rng()%60)-30,xr=Int(rng()%60)-30,yl=Int(rng()%60)-30,yr=Int(rng()%60)-30;if(xl>xr)std::swap(xl,xr);if(yl>yr)std::swap(yl,yr);Int sum=0,pre=0,less=0,all=0;for(auto [pos,val]:ref){auto [x,y]=pos;all+=val;if(x<xr&&y<yr)pre+=val;if(xl<=x&&x<xr&&y<yr)less+=val;if(xl<=x&&x<xr&&yl<=y&&y<yr)sum+=val;}assert(s.get(xl,xr,yl,yr)==sum&&ss.get(xl,xr,yl,yr)==sum&&zero.get(xl,xr,yl,yr)==sum&&s.get_all()==all&&ss.get_all()==all&&zero.get_all()==all);assert(t.get(xl,xr,yl,yr)==sum&&z.get(xl,xr,yl,yr)==sum);assert(t.prefix(xr,yr)==pre&&t.getLess(xl,xr,yr)==less&&t.get_all()==all&&z.get_all()==all);}
    }
    for(Int h:{0,1,3,7})for(Int w:{0,1,5,8}){std::vector<std::vector<Int>> a(h,std::vector<Int>(w));auto t=initSegmentTree2D(a,[](Int x,Int y){return x+y;},Int(0));for(int q=0;q<300;++q){if(h&&w){Int x=rng()%h,y=rng()%w;a[x][y]=rng()%100;t.update(x,y,a[x][y]);}Int xl=rng()%(h+1),xr=rng()%(h+1),yl=h?rng()%(w+1):0,yr=h?rng()%(w+1):0;if(xl>xr)std::swap(xl,xr);if(yl>yr)std::swap(yl,yr);Int sum=0;for(Int i=xl;i<xr;++i)for(Int j=yl;j<yr;++j)sum+=a[i][j];assert(t.get(xl,xr,yl,yr)==sum);}}
    auto e=initCompressedFenwick2D(std::vector<std::pair<Int,Int>>{});assert(e.len()==0&&e.get_all()==0&&e.prefix(0,0)==0&&e.get(0,0)==0);
}
void dynamic_persistent_test(){
    std::mt19937_64 rng(19);Int n=Int(1)<<60;auto s=initDynamicSegmentTree(n,[](std::string a,const std::string& b){return a+b;},std::string{});std::map<Int,std::string> ref;
    for(int q=0;q<2000;++q){Int x=rng()%n;if(q%3==0&&!ref.empty())x=ref.begin()->first;std::string v=std::to_string(q)+",";s[x]=v;ref[x]=v;Int l=rng()%n,r=rng()%n;if(l>r)std::swap(l,r);std::string want,all;for(auto& [key,value]:ref){if(l<=key&&key<r)want+=value;all+=value;}assert(s.get(l,r)==want&&s.get_all()==all&&s.node_count()==Int(ref.size()));assert(std::string(s[x])==v);}
    s[from_end(1)]="end";assert(std::string(s[from_end(1)])=="end");auto copy=s;copy[0]="copy";assert(std::string(s[0]).empty());assert(copy.node_count()==s.node_count()+1);

}
void dynamic_lazy_test(){
    using S=RangeSum<Int>;using F=RangeAffine<Int>;auto op=[](S a,S b){return S{a.sum+b.sum,a.len+b.len};};auto map=[](F f,S s){return S{f.a*s.sum+f.b*s.len,s.len};};auto comp=[](F f,F g){return F{f.a*g.a,f.a*g.b+f.b};};auto init=[](Int l,Int r){return S{0,r-l};};std::mt19937 rng(20);
    for(Int n:{0,1,2,5,64,129,500}){auto s=initDynamicLazySegmentTree(n,op,S{},map,comp,F{1,0},init);std::vector<Int> a(n);std::set<Int> boundaries{0,n};assert(s.node_count()==Int(n>0));for(int q=0;q<2000;++q){Int l=rng()%(n+1),r=rng()%(n+1);if(l>r)std::swap(l,r);F f{Int(rng()%3)-1,Int(rng()%11)-5};if(q%4){Int before=s.node_count();s.apply(l,r,f);assert(s.node_count()<=before+2);if(l<r){boundaries.insert(l);boundaries.insert(r);}for(Int i=l;i<r;++i)a[i]=f.a*a[i]+f.b;}else if(n){Int p=rng()%n,x=Int(rng()%20)-10;s[p]=S{x,1};a[p]=x;boundaries.insert(p);boundaries.insert(p+1);}assert(s.node_count()==Int(boundaries.size())-1);l=rng()%(n+1);r=rng()%(n+1);if(l>r)std::swap(l,r);Int before=s.node_count(),sum=std::accumulate(a.begin()+l,a.begin()+r,Int(0));auto result=s.get(l,r);assert(result.sum==sum&&result.len==r-l);assert(s.node_count()==before);assert(s.get_all().sum==std::accumulate(a.begin(),a.end(),Int(0)));}auto copy=s;if(n){copy[from_end(1)]=S{99,1};assert(S(copy[from_end(1)]).sum==99&&S(s[from_end(1)]).sum==a.back());}}
    {Int n=Int(1)<<50;auto s=initDynamicLazySegmentTree(n,op,S{},map,comp,F{1,0},init);s.apply(123,n-456,F{1,3});assert(s.node_count()==3);assert(s.get(125,n-500).sum==(n-625)*3);assert(s.node_count()==3);s.apply(123,n-456,F{0,7});assert(s.node_count()==3&&s.get_all().sum==(n-579)*7);}
    {auto op=[](std::string a,std::string b){return a+b;};auto map=[](int f,std::string s){if(f)std::fill(s.begin(),s.end(),char(f));return s;};auto comp=[](int f,int g){return f?f:g;};auto init=[](Int l,Int r){std::string s;for(Int i=l;i<r;++i)s+=char('a'+i);return s;};auto t=initDynamicLazySegmentTree(20,op,std::string{},map,comp,0,init);auto a=init(0,20);t.apply(3,17,'z');std::fill(a.begin()+3,a.begin()+17,'z');t[8]=std::string("x");a[8]='x';for(Int l=0;l<=20;++l)for(Int r=l;r<=20;++r)assert(t.get(l,r)==a.substr(l,r-l));}
}
