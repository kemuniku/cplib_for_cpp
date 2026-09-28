#pragma once
#include <cplib/graph/graph.hpp>
namespace cplib {
// 通行可能な上下左右の隣接マスを結ぶ。壁の頂点も残す。O(HW)。静的版は呼出側でbuildする。
template<bool ReturnStatic=false,class Grid,class T> auto grid_to_graph_impl(const Grid& a,const T& ok){Int h=a.size(),w=h?Int(a[0].size()):0;BasicGraph<void,ReturnStatic,false> out(h*w);for(Int i=0;i<h;++i)for(Int j=0;j<w;++j)if(a[i][j]==ok){if(i+1<h&&a[i+1][j]==ok)out.add_edge(i*w+j,(i+1)*w+j);if(j+1<w&&a[i][j+1]==ok)out.add_edge(i*w+j,i*w+j+1);}return out;}
template<bool ReturnStatic=false,class Grid,class T=char> auto grid_to_graph(const Grid& a,const T& ok=T('.')){return grid_to_graph_impl<ReturnStatic>(a,ok);}
}
