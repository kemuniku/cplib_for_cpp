#pragma once
#include <cplib/collections/segtree.hpp>
namespace cplib {
// 座標圧縮とセグメント木で転倒数を数える。O(N log N)。
inline Int inversion_number(std::span<const Int> a){std::vector<Int> coords(a.begin(),a.end());std::sort(coords.begin(),coords.end());coords.erase(std::unique(coords.begin(),coords.end()),coords.end());auto seg=initSegmentTree(Int(coords.size()),[](Int l,Int r){return l+r;},Int(0));Int result=0;for(Int value:a){Int p=std::lower_bound(coords.begin(),coords.end(),value)-coords.begin();if(p+1<Int(coords.size()))result+=seg.get(p+1,coords.size());seg[p]=Int(seg[p])+1;}return result;}
}
