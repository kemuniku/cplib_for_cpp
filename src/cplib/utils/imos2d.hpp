#pragma once
#include <cplib/common.hpp>
namespace cplib {
class Imos2D {
    std::vector<std::vector<Int>> B;Int H,W;
public:
    // 二次元差分配列を構築する。O(HW)。
    Imos2D(Int h,Int w):B(h+1,std::vector<Int>(w+1)),H(h),W(w){}
    // 半開矩形に加算する。O(1)。
    void rectangle_add(Int il,Int ir,Int jl,Int jr,Int x){B[il][jl]+=x;B[il][jr]-=x;B[ir][jl]-=x;B[ir][jr]+=x;}
    // 加算結果を復元する。O(HW)、何度でも呼べる。
    auto build() const{
        std::vector<std::vector<Int>> result(H,std::vector<Int>(W));
        for(Int i=0;i<H;++i)for(Int j=0;j<W;++j){result[i][j]=B[i][j];if(i)result[i][j]+=result[i-1][j];if(j)result[i][j]+=result[i][j-1];if(i&&j)result[i][j]-=result[i-1][j-1];}
        return result;
    }
};
inline Imos2D initImos2D(Int h,Int w){return {h,w};}
inline void rectangle_add(Imos2D& s,Int il,Int ir,Int jl,Int jr,Int x){s.rectangle_add(il,ir,jl,jr,x);}
inline auto build(const Imos2D& s){return s.build();}
}
