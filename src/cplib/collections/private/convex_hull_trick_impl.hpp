#pragma once
#include <cplib/common.hpp>
#include <deque>
namespace cplib {
struct CHTLine {Int a=0,b=0;__int128 start=0;friend bool operator<(const CHTLine& l,const CHTLine& r){return l.a>r.a;}friend bool operator<=(const CHTLine& l,const CHTLine& r){return l.a>=r.a;}};
inline __int128 chtValue(const CHTLine& line,Int x){return __int128(line.a)*x+line.b;}
inline Int chtAnswer(__int128 v){assert(v>=std::numeric_limits<Int>::min()&&v<=std::numeric_limits<Int>::max());return Int(v);}
// 傾きの小さい直線が最小になる最初の整数座標。
inline __int128 chtStart(const CHTLine& l,const CHTLine& r){__int128 n=__int128(r.b)-l.b,d=__int128(l.a)-r.a;assert(d>0);return n/d+(n%d>0);}
inline bool chtRedundant(const CHTLine& l,const CHTLine& m,const CHTLine& r){return chtStart(l,m)>=chtStart(m,r);}
struct CHTMonotoneHull {
    std::deque<CHTLine> lines;
    bool slopeIncreasing=false,hasSlope=false;Int lastSlope=0;
    explicit CHTMonotoneHull(bool increasing=false):slopeIncreasing(increasing){}
    // 傾きが単調な直線追加。償却 O(1)。
    void chtAddLine(Int a,Int b){if(hasSlope)assert(slopeIncreasing?lastSlope<=a:a<=lastSlope);hasSlope=true;lastSlope=a;CHTLine line{a,b};
        if(slopeIncreasing){if(!lines.empty()&&lines.front().a==a){if(lines.front().b<=b)return;lines.pop_front();}while(lines.size()>=2&&chtRedundant(line,lines[0],lines[1]))lines.pop_front();lines.push_front(line);}
        else{if(!lines.empty()&&lines.back().a==a){if(lines.back().b<=b)return;lines.pop_back();}while(lines.size()>=2&&chtRedundant(lines[lines.size()-2],lines.back(),line))lines.pop_back();lines.push_back(line);}
    }
};
inline CHTMonotoneHull initCHTMonotoneHull(bool increasing=false){return CHTMonotoneHull(increasing);}
inline void chtAddLine(CHTMonotoneHull& hull,Int a,Int b){hull.chtAddLine(a,b);}
}
