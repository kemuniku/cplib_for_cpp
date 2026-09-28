#pragma once
#include <cplib/common.hpp>
namespace cplib {
// 最大要素の最初の位置を返す。空なら(-1,-1)、O(全要素数)。
template<class Grid> std::pair<Int,Int> maxIndex(const Grid& v){std::pair<Int,Int> r{-1,-1};for(Int i=0;i<Int(v.size());++i)for(Int j=0;j<Int(v[i].size());++j)if(r.first==-1 || v[i][j]>v[r.first][r.second])r={i,j};return r;}
// 最小要素の最初の位置を返す。空なら(-1,-1)、O(全要素数)。
template<class Grid> std::pair<Int,Int> minIndex(const Grid& v){std::pair<Int,Int> r{-1,-1};for(Int i=0;i<Int(v.size());++i)for(Int j=0;j<Int(v[i].size());++j)if(r.first==-1 || v[i][j]<v[r.first][r.second])r={i,j};return r;}
}
