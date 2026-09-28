#pragma once
#include <cplib/collections/private/convex_hull_trick_impl.hpp>
namespace cplib {
class ConvexHullTrickMonotone {
    CHTMonotoneHull hull_;bool xIncreasing_=true,hasX_=false;Int lastX_=0;
public:
    explicit ConvexHullTrickMonotone(bool slopeIncreasing=false,bool xIncreasing=true):hull_(slopeIncreasing),xIncreasing_(xIncreasing){}
    void add_line(Int a,Int b){hull_.chtAddLine(a,b);}
    // 単調な座標を走査し、不要な直線を取り除く。償却 O(1)。
    Int get_min(Int x){auto& lines=hull_.lines;assert(!lines.empty());if(hasX_)assert(xIncreasing_?lastX_<=x:x<=lastX_);hasX_=true;lastX_=x;if(xIncreasing_){while(lines.size()>=2&&chtValue(lines[0],x)>=chtValue(lines[1],x))lines.pop_front();return chtAnswer(chtValue(lines.front(),x));}else{while(lines.size()>=2&&chtValue(lines.back(),x)>=chtValue(lines[lines.size()-2],x))lines.pop_back();return chtAnswer(chtValue(lines.back(),x));}}
};
inline ConvexHullTrickMonotone initConvexHullTrickMonotone(bool slopeIncreasing=false,bool xIncreasing=true){return ConvexHullTrickMonotone(slopeIncreasing,xIncreasing);}
}
