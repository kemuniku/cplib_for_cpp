#ifndef CPLIB_WARSHALL_HEADER
#define CPLIB_WARSHALL_HEADER <cplib/graph/warshall_floyd.hpp>
#endif
#include CPLIB_WARSHALL_HEADER
#include <random>
#include <iostream>
using namespace cplib;
template<class T> auto reference(const std::vector<std::vector<T>>& a,T inf){
    Int n=a.size();std::vector<std::vector<T>> result(n,std::vector<T>(n,inf));
    for(Int source=0;source<n;++source){
        auto& d=result[source];d[source]=0;std::vector<bool> bad(n);
        for(Int step=0;step<n;++step){
            auto next=d;
            for(Int u=0;u<n;++u)for(Int v=0;v<n;++v)if(d[u]!=inf && a[u][v]!=inf && d[u]+a[u][v]<next[v]){next[v]=d[u]+a[u][v];if(step==n-1)bad[v]=true;}
            d=next;
        }
        for(Int step=0;step<n;++step)for(Int u=0;u<n;++u)if(bad[u])for(Int v=0;v<n;++v)if(a[u][v]!=inf)bad[v]=true;
        for(Int v=0;v<n;++v)if(bad[v])d[v]=-inf;
    }
    return result;
}
template<class T> void check_type(){
    T inf=detail::distance_inf<T>();std::mt19937 rng(71);
    for(int trial=0;trial<160;++trial){
        Int n=rng()%20;std::vector<std::vector<T>> a(n,std::vector<T>(n,inf));
        for(Int i=0;i<n;++i)for(Int j=0;j<n;++j)if(rng()%4==0)a[i][j]=T(Int(rng()%30)-10);
        auto want=reference(a,inf);assert(cplib::warshall_floyd(a)==want);
    }
    for(Int n:{Int(31),Int(217),Int(257),Int(451)})for(bool dense:{false,true}){
        std::vector<std::vector<T>> a(n,std::vector<T>(n,inf));
        for(Int i=0;i<n;++i)for(Int j=0;j<n;++j)if(dense || rng()%20==0)a[i][j]=T(rng()%100);
        auto want=a;detail::warshall_inplace<detail::ScalarWarshall,true>(want,T(0),inf);assert(cplib::warshall_floyd_nonnegative(a)==want);
    }
}
int main(){
    auto g=initWeightedUnDirectedGraph(4);assert(g.add_edge(0,1,4)==0);assert(g.add_edge(1,2,3)==1);g.add_edge(0,1,8);
    auto s=initWeightedUnDirectedStaticGraph(4);for(auto e:g.edge_info)s.add_edge(e.src,e.dst,e.cost);s.build();
    for(Int i=0;i<4;++i){std::vector<std::tuple<Int,Int,Int>> a,b;for(auto e:g.to_and_cost_and_id(i))a.push_back(e);for(auto e:s.to_and_cost_and_id(i))b.push_back(e);assert(a==b);}
    assert(warshall_floyd(g)==warshall_floyd(s));assert(warshall_floyd(g)[0][2]==7);
    auto mixed=initUnWeightedDirectedStaticGraph(3);add_edge_static_impl(mixed,0,1,true);add_edge_static_impl(mixed,1,2,false);build(mixed);
    assert(warshall_floyd(mixed)[2][0]==INF64);assert(warshall_floyd(mixed)[2][1]==1);
    auto labels=initWeightedUnDirectedTableGraph(std::vector<std::string>{"a","b","c"});add_edge(labels,"a","c",Int(8));
    for(auto [v,c]:labels["c"]){assert(v=="a" && c==8);}
    for(auto [v,c,id]:labels.to_and_cost_and_id("a")){assert(v=="c" && c==8 && id==0);}
    check_type<Int>();check_type<std::int32_t>();check_type<double>();check_type<float>();
    std::cout<<"graph tests passed; AVX2="<<bool(__builtin_cpu_supports("avx2"))<<", AVX512F="<<bool(__builtin_cpu_supports("avx512f"))<<'\n';
}
