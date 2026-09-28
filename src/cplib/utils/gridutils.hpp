#pragma once
#include <cplib/common.hpp>
#include <cplib/iterator_views.hpp>
namespace cplib {
inline constexpr std::array<std::pair<Int,Int>,4> dij4{{{0,1},{-1,0},{0,-1},{1,0}}};
inline constexpr std::array<std::pair<Int,Int>,8> dij8{{{0,1},{-1,1},{-1,0},{-1,-1},{0,-1},{1,-1},{1,0},{1,1}}};
inline constexpr std::string_view dij4_str="RULD";
// 範囲内の隣接マスを列挙する。O(方向数)、追加領域O(1)。
inline auto griditer(Int i,Int j,Int h,Int w,std::span<const std::pair<Int,Int>> dij=dij4){
    return detail::filter_range(detail::transform_range(dij,[=](auto d){return std::pair{i+d.first,j+d.second};}),[=](auto p){return 0<=p.first&&p.first<h&&0<=p.second&&p.second<w;});
}
// 一致する最初のマスを返す。なければ(-1,-1)、O(全マス数)。
template<class Grid,class T> std::pair<Int,Int> gridfind(const Grid& grid,const T& value){for(Int i=0;i<Int(grid.size());++i)for(Int j=0;j<Int(grid[i].size());++j)if(grid[i][j]==value)return {i,j};return {-1,-1};}
template<class Grid,class T> auto gridfinds(const Grid& grid,const T& value){std::vector<std::pair<Int,Int>> r;for(Int i=0;i<Int(grid.size());++i)for(Int j=0;j<Int(grid[i].size());++j)if(grid[i][j]==value)r.emplace_back(i,j);return r;}
template<class Grid> Int height(const Grid& grid){return grid.size();}
template<class Grid> Int width(const Grid& grid){return grid.empty()?0:grid[0].size();}
template<class Grid> Int getid(const Grid& grid,Int i,Int j){return i*width(grid)+j;}
template<class Grid> std::pair<Int,Int> to_pos(const Grid& grid,Int id){return {id/width(grid),id%width(grid)};}
}
