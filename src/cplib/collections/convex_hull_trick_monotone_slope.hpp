#pragma once
#include <cplib/collections/private/convex_hull_trick_impl.hpp>
namespace cplib {
class ConvexHullTrickMonotoneSlope {
    CHTMonotoneHull hull_;
public:
    explicit ConvexHullTrickMonotoneSlope(bool increasing=false):hull_(increasing){}
    void add_line(Int a,Int b){hull_.chtAddLine(a,b);}
    // 任意座標での最小値を二分探索。O(log N)。
    Int get_min(Int x)const{assert(!hull_.lines.empty());std::size_t l=0,r=hull_.lines.size()-1;while(l<r){auto m=l+(r-l)/2;if(chtValue(hull_.lines[m],x)>=chtValue(hull_.lines[m+1],x))l=m+1;else r=m;}return chtAnswer(chtValue(hull_.lines[l],x));}
};
inline ConvexHullTrickMonotoneSlope initConvexHullTrickMonotoneSlope(bool slopeIncreasing=false){return ConvexHullTrickMonotoneSlope(slopeIncreasing);}
}
