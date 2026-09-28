#include <cplib/graph/planar_graph.hpp>
#include <fstream>
#include <iostream>
using namespace cplib;
int main(int argc,char** argv){assert(argc==2);std::ifstream in(argv[1]);Int n,m,expected;while(in>>n>>m>>expected){auto g=initUnWeightedUnDirectedGraph(n);auto s=initWeightedUnDirectedStaticGraph<Int>(n);for(Int i=0;i<m;++i){Int u,v;in>>u>>v;g.add_edge(u,v);s.add_edge(u,v,i);}assert(is_planar_graph(g)==bool(expected));assert(is_planar_graph(s)==bool(expected));s.build();assert(is_planar_graph(s)==bool(expected));}auto path=initUnWeightedUnDirectedGraph(200000);for(Int i=1;i<path.len;++i)path.add_edge(i-1,i);assert(is_planar_graph(path));}
