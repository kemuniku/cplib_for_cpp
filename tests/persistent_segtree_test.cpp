#include <cplib/collections/persistent_segtree.hpp>
#include <random>
using namespace cplib;
int main(){std::mt19937_64 rng(19);
    for(Int size:{0,1,3,16,17,100}){std::vector<Int> a(size);std::iota(a.begin(),a.end(),0);auto s=cplib::initSegmentTree(a,[](Int x,Int y){return x+y;},Int(0));std::vector<decltype(s)> versions{s};std::vector<std::vector<Int>> refs{a};for(int q=0;size&&q<500;++q){Int v=rng()%versions.size(),i=rng()%size,x=rng()%100;auto b=refs[v];b[i]=x;versions.push_back(versions[v].update(i,x));refs.push_back(b);v=rng()%versions.size();Int l=rng()%(size+1),r=rng()%(size+1);if(l>r)std::swap(l,r);assert(versions[v].query(l,r)==std::accumulate(refs[v].begin()+l,refs[v].begin()+r,Int(0)));}assert(s.query(0,size)==size*(size-1)/2);}
}
