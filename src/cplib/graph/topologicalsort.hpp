#pragma once
#include <cplib/graph/graph.hpp>
namespace cplib {
// スタックでトポロジカル順序を求める。O(V+E)。閉路部分は返さない。
template<DirectedGraph G> std::vector<Int> topologicalsort(const G& g){
    std::vector<Int> gin(g.len),stack,result;
    for(Int i=0;i<g.len;++i)for(auto [j,c]:g.to_and_cost(i))++gin[j];
    for(Int i=0;i<g.len;++i)if(!gin[i])stack.push_back(i);
    while(!stack.empty()){Int i=stack.back();stack.pop_back();result.push_back(i);for(auto [j,c]:g.to_and_cost(i))if(--gin[j]==0)stack.push_back(j);}
    return result;
}
template<DirectedGraph G> bool isDAG(const G& g){return Int(topologicalsort(g).size())==g.len;}
}
