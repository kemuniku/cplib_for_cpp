#pragma once
#include <cplib/graph/topologicalsort.hpp>
namespace cplib {
// トポロジカル順序の数を求める。O(V2^V+E)時間・O(2^V+V)領域。
template<DirectedGraph G> Int count_topologicalsort(const G& g){
    Int n=g.len;assert(n<63);if(!isDAG(g))return 0;std::vector<UInt> pred(n);
    for(Int u=0;u<n;++u)for(auto [v,c]:g.to_and_cost(u))pred[v]|=UInt(1)<<u;
    UInt size=UInt(1)<<n;std::vector<Int> dp(size);dp[0]=1;
    for(UInt mask=0;mask<size;++mask)if(dp[mask])for(Int v=0;v<n;++v){UInt bit=UInt(1)<<v;if(!(mask&bit) && (mask&pred[v])==pred[v])dp[mask|bit]+=dp[mask];}
    return dp.back();
}
}
