#pragma once
#include <cplib/matrix/bit_matrix_ops.hpp>
#include <bit>
#include <memory>
#include <string_view>
namespace cplib {
template<Int H,Int W> class StaticMatrixMod2 {
 static_assert(H>=0&&W>=0);template<Int,Int> friend class StaticMatrixMod2;
 static constexpr Int Stride=(W+63)/64;
 std::array<std::array<UInt,Stride>,H> rows_{};
public:
 using value_type=bool;
 StaticMatrixMod2()=default;
 template<std::integral T> explicit StaticMatrixMod2(const std::array<std::array<T,W>,H>& rows){for(Int i=0;i<H;++i)for(Int j=0;j<W;++j)if(rows[i][j]&1)rows_[i][j>>6]|=UInt(1)<<(j&63);}
 constexpr Int h()const{return H;}constexpr Int w()const{return W;}
 bool operator()(Int i,Int j)const{assert(0<=i&&i<H&&0<=j&&j<W);return (rows_[i][j>>6]>>(j&63))&1;}
 detail::MatrixBitReference operator()(Int i,Int j){assert(0<=i&&i<H&&0<=j&&j<W);return {&rows_[i][j>>6],UInt(1)<<(j&63)};}
 struct Row{StaticMatrixMod2* matrix;Int row;auto operator[](Int col){return (*matrix)(row,col);}};struct ConstRow{const StaticMatrixMod2* matrix;Int row;bool operator[](Int col)const{return (*matrix)(row,col);}};
 Row operator[](Int row){return {this,row};}ConstRow operator[](Int row)const{return {this,row};}
 bool operator==(const StaticMatrixMod2&)const=default;
 std::string str()const{std::string out;for(Int i=0;i<H;++i){if(i)out+='\n';for(Int j=0;j<W;++j){if(j)out+=' ';out+=(*this)(i,j)?'1':'0';}}return out;}
 void setRowBits(Int i,std::string_view values){assert(0<=i&&i<H&&Int(values.size())<=W);rows_[i].fill(0);for(Int j=0;j<Int(values.size());++j)if(values[j]=='1')rows_[i][j>>6]|=UInt(1)<<(j&63);}
 std::string rowBits(Int i,Int width=W)const{assert(0<=i&&i<H&&0<=width&&width<=W);std::string out(width,'0');for(Int j=0;j<width;++j)out[j]=char('0'+((rows_[i][j>>6]>>(j&63))&1));return out;}
 static StaticMatrixMod2 identity() requires(H==W){StaticMatrixMod2 out;for(Int i=0;i<H;++i)out(i,i)=true;return out;}
 auto transposed()const{StaticMatrixMod2<W,H> out;for(Int i=0;i<H;++i)for(Int j=0;j<W;++j)if((*this)(i,j))out(j,i)=true;return out;}
 // 固定長版は分岐とFour Russians法の表の寸法もコンパイル時に決める。
 template<Int C> auto operator*(const StaticMatrixMod2<W,C>& b)const{
  StaticMatrixMod2<H,C> out;
  if constexpr(H<40||W<64){auto bt=b.transposed();for(Int i=0;i<H;++i)for(Int j=0;j<C;++j){int parity=0;for(Int k=0;k<Stride;++k)parity^=std::popcount(rows_[i][k]&bt.rows_[j][k])&1;if(parity)out(i,j)=true;}}
  else{constexpr Int count=(C+63)/64;std::array<std::array<UInt,count>,256> table{};for(Int start=0;start<W;start+=8){Int bits=std::min<Int>(8,W-start);for(unsigned mask=1;mask<(1u<<bits);++mask){Int previous=mask&(mask-1),bit=std::countr_zero(mask);for(Int k=0;k<count;++k)table[mask][k]=table[previous][k]^b.rows_[start+bit][k];}UInt mask=(UInt(1)<<bits)-1;for(Int i=0;i<H;++i){Int index=(rows_[i][start>>6]>>(start&63))&mask;for(Int k=0;k<count;++k)out.rows_[i][k]^=table[index][k];}}}
  return out;
 }
 StaticMatrixMod2& operator*=(const StaticMatrixMod2& b) requires(H==W){*this=(*this)*b;return *this;}
 StaticMatrixMod2 pow(Int exponent)const requires(H==W){assert(exponent>=0);auto out=identity(),base=*this;for(;exponent;exponent>>=1){if(exponent&1)out*=base;if(exponent>1)base*=base;}return out;}
 // 元の固定長版と同じくrankの作業行列はヒープに置く。
 Int rank()const{if constexpr(H==0||W==0)return 0;auto storage=std::make_unique<StaticMatrixMod2>(*this);auto& b=*storage;Int result=0;for(Int col=0;col<W;++col){Int pivot=result;while(pivot<H&&!b(pivot,col))++pivot;if(pivot==H)continue;std::swap(b.rows_[result],b.rows_[pivot]);for(Int i=result+1;i<H;++i)if(b(i,col))for(Int k=0;k<Stride;++k)b.rows_[i][k]^=b.rows_[result][k];if(++result==H)break;}return result;}
 bool determinant()const requires(H==W){return rank()==H;}
 // 左右の固定長行列を別々に掃き出す。動的版の拡大配列とは別実装。
 std::optional<StaticMatrixMod2> inverse()const requires(H==W){auto left=*this,right=identity();for(Int col=0;col<H;++col){Int pivot=col;while(pivot<H&&!left(pivot,col))++pivot;if(pivot==H)return std::nullopt;std::swap(left.rows_[col],left.rows_[pivot]);std::swap(right.rows_[col],right.rows_[pivot]);for(Int i=0;i<H;++i)if(i!=col&&left(i,col)){for(Int k=0;k<Stride;++k)left.rows_[i][k]^=left.rows_[col][k];for(Int k=0;k<Stride;++k)right.rows_[i][k]^=right.rows_[col][k];}}return right;}
 template<class B> auto solveLinearSystem(const B& b,Int height=H,Int width=W)const{assert(0<=height&&height<=H&&0<=width&&width<=W&&Int(b.size())==height);auto rows=initBitLinearSystem(height,width);Int stride=(width>>6)+1,full=width>>6,tail=width&63;for(Int i=0;i<height;++i){for(Int k=0;k<full;++k)rows[i*stride+k]=rows_[i][k];if(tail)rows[i*stride+full]=rows_[i][full]&((UInt(1)<<tail)-1);if(b[i])rows[i*stride+full]|=UInt(1)<<tail;}return solveBitLinearSystem(rows,height,width);}
 bool hafnian()const{assert(H==W);return fieldHafnian(matrixRows(*this,H,W));}
 StaticMatrixMod2 adjugate()const{assert(H==W);auto rows=*fieldAdjugateInverse(matrixRows(*this,H,W),true);StaticMatrixMod2 out;for(Int i=0;i<H;++i)for(Int j=0;j<W;++j)out(i,j)=bool(rows[i][j]);return out;}
};
template<Int H,Int W> auto initStaticMatrixMod2(){return StaticMatrixMod2<H,W>{};}
template<std::size_t H,std::size_t W,std::integral T> auto initStaticMatrixMod2(const std::array<std::array<T,W>,H>& rows){return StaticMatrixMod2<H,W>(rows);}
template<std::size_t H,std::size_t W,std::integral T> auto toStaticMatrixMod2(const std::array<std::array<T,W>,H>& rows){return initStaticMatrixMod2(rows);}
template<Int N> auto identityStaticMatrixMod2(){return StaticMatrixMod2<N,N>::identity();}
template<Int H,Int W> Int h(const StaticMatrixMod2<H,W>& a){return a.h();}template<Int H,Int W> Int w(const StaticMatrixMod2<H,W>& a){return a.w();}
template<Int H,Int W> auto to_string(const StaticMatrixMod2<H,W>& a){return a.str();}
template<Int H,Int W> void setRowBits(StaticMatrixMod2<H,W>& a,Int i,std::string_view values){a.setRowBits(i,values);}
template<Int H,Int W> auto rowBits(const StaticMatrixMod2<H,W>& a,Int i,Int width=W){return a.rowBits(i,width);}
template<Int H,Int W> auto transposed(const StaticMatrixMod2<H,W>& a){return a.transposed();}
template<Int N> auto pow(const StaticMatrixMod2<N,N>& a,Int exponent){return a.pow(exponent);}
template<Int H,Int W> Int rank(const StaticMatrixMod2<H,W>& a){return a.rank();}
template<Int N> bool determinant(const StaticMatrixMod2<N,N>& a){return a.determinant();}
template<Int N> auto inverse(const StaticMatrixMod2<N,N>& a){return a.inverse();}
template<Int H,Int W> bool hafnian(const StaticMatrixMod2<H,W>& a){return a.hafnian();}
template<Int H,Int W> auto adjugate(const StaticMatrixMod2<H,W>& a){return a.adjugate();}
template<Int H,Int W,class B> auto solveLinearSystem(const StaticMatrixMod2<H,W>& a,const B& b,Int height=H,Int width=W){return a.solveLinearSystem(b,height,width);}
}
