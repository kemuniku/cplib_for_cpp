#include CPLIB_SEGTREE_HEADER
#include <random>
using namespace cplib;
int main(){
    std::mt19937 rng(12);auto add=[](Int a,Int b){return a+b;};
    for(Int n:{0,1,2,3,17,63,64,65,500}){
        std::vector<Int> a(n);for(auto& x:a)x=rng()%20;
        auto s=initSegmentTree(a,add,Int(0));auto st=cplib::initSegmentTree(a,add,Int(0));auto sv=cplib::initSegmentTree(a,add,Int(0));
        auto check=[&](auto& tree){assert(tree.len()==n);assert(tree.get_all()==std::accumulate(a.begin(),a.end(),Int(0)));for(int q=0;q<30;++q){Int l=rng()%(n+1),r=rng()%(n+1);if(l>r)std::swap(l,r);assert(tree.get(l,r)==std::accumulate(a.begin()+l,a.begin()+r,Int(0)));assert(tree[closed_slice(l,r-1)]==tree.get(l,r));Int cap=rng()%80;Int rr=l,sum=0;while(rr<n&&sum+a[rr]<=cap)sum+=a[rr++];assert(tree.max_right(l,[&](Int v){return v<=cap;})==rr);Int ll=r;sum=0;while(ll>0&&sum+a[ll-1]<=cap)sum+=a[--ll];assert(tree.min_left(r,[&](Int v){return v<=cap;})==ll);}};
        for(int q=0;q<100;++q){if(n){Int p=rng()%n,val=rng()%30;s[p]=val;st[p]=val;sv[p]=val;a[p]=val;}check(s);check(st);check(sv);}

#ifdef CPLIB_SEGTREE_VAR
        if(n){sv[from_end(1)]+=2;a[n-1]+=2;check(sv);auto copy=sv;auto moved=std::move(copy);moved[0]+=3;assert(moved.get_all()==sv.get_all()+3);auto assigned=cplib::initSegmentTree(Int(0),add,Int(0));assigned=sv;assigned[0]*=2;assert(assigned.get_all()==sv.get_all()+a[0]);assert(sv.get_all()==std::accumulate(a.begin(),a.end(),Int(0)));}
#endif
    }
    auto cat=[](std::string a,const std::string& b){return a+b;};std::vector<std::string> a={"ab","c","de","f","g"};auto s=initSegmentTree(a,cat,std::string{});auto st=cplib::initSegmentTree(a,cat,std::string{});auto sv=cplib::initSegmentTree(a,cat,std::string{});
    for(Int l=0;l<=5;++l)for(Int r=l;r<=5;++r){std::string want;for(Int i=l;i<r;++i)want+=a[i];assert(s.get(l,r)==want&&st.get(l,r)==want&&sv.get(l,r)==want);}
    for(Int l=0;l<5;++l){std::string want=s.get(l,4);assert(s.max_right(l,[&](std::string x){return want.starts_with(x);})==4);}
    for(Int r=1;r<=5;++r){std::string want=s.get(1,r);assert(s.min_left(r,[&](std::string x){return want.ends_with(x);})==1);}

#ifdef CPLIB_SEGTREE_VAR
    auto v=cplib::initSegmentTree(std::vector<Int>{5},[](Int x,Int y){return x+y;},Int(0));v[0]<<=2;v[0]>>=1;v[0]|=1;v[0]&=7;v[0]^=2;v[0]%=5;v[0].pow_assign(3);assert(v.get_all()==1);v[0]=-5;v[0].floor_div_assign(2);assert(v.get_all()==-3);
#endif
}
