#pragma once
#include <cplib/common.hpp>
#include <optional>
namespace cplib {
template<class T> struct LinearSystemSolution {std::vector<T> particular;std::vector<std::vector<T>> basis;};
namespace detail {
template<class T> T field_sub(T a,T b){if constexpr(std::is_same_v<T,bool>)return a!=b;else return a-b;}
template<class T> T field_mul(T a,T b){if constexpr(std::is_same_v<T,bool>)return a&&b;else return a*b;}
template<class T> bool field_equal(T a,T b){if constexpr(requires{a.val();a.umod();})return a.val()==b.val();else return a==b;}
template<class T> T field_inv(T a){if constexpr(std::is_same_v<T,bool>)return true;else return T(1)/a;}
template<class T> struct FieldEchelonResult {std::vector<Int> pivots;T determinant=T(1);};
// 体上の前進消去。零判定はMontgomery表現の場合も公開値で行う。
template<class T> auto field_echelon(std::vector<std::vector<T>>& a,Int columns){
 FieldEchelonResult<T> result;
 for(Int col=0;col<columns;++col){Int row=result.pivots.size();if(row==Int(a.size()))break;Int pivot=row;while(pivot<Int(a.size())&&field_equal<T>(a[pivot][col],T(0)))++pivot;if(pivot==Int(a.size()))continue;if(pivot!=row){std::swap(a[pivot],a[row]);result.determinant=field_sub(T(0),result.determinant);}T value=a[row][col];result.determinant=field_mul(result.determinant,value);T inverse=field_inv(value);
  for(Int j=col;j<Int(a[row].size());++j)a[row][j]=field_mul<T>(a[row][j],inverse);
  for(Int i=row+1;i<Int(a.size());++i){T factor=a[i][col];if(field_equal(factor,T(0)))continue;a[i][col]=T(0);for(Int j=col+1;j<Int(a[i].size());++j)a[i][j]=field_sub<T>(a[i][j],field_mul<T>(factor,a[row][j]));}result.pivots.push_back(col);
 }
 return result;
}
}
template<class M> auto matrixRows(const M& a,Int height,Int width){assert(0<=height&&height<=a.h()&&0<=width&&width<=a.w());using T=typename M::value_type;std::vector<std::vector<T>> rows(height,std::vector<T>(width));for(Int i=0;i<height;++i)for(Int j=0;j<width;++j)rows[i][j]=a(i,j);return rows;}
template<class T> Int fieldRank(std::vector<std::vector<T>> rows,Int width){return detail::field_echelon(rows,width).pivots.size();}
template<class T> T fieldDeterminant(std::vector<std::vector<T>> rows){auto e=detail::field_echelon(rows,rows.size());return e.pivots.size()==rows.size()?e.determinant:T(0);}
// 特殊解と核の基底。O(hw min(h,w)+w² min(h,w))。
template<class T,class B> std::optional<LinearSystemSolution<T>> fieldSolve(std::vector<std::vector<T>> rows,Int width,const B& b){
 using namespace detail;assert(b.size()==rows.size());for(std::size_t i=0;i<rows.size();++i)rows[i].push_back(b[i]);auto pivots=field_echelon(rows,width).pivots;for(Int i=pivots.size();i<Int(rows.size());++i)if(!field_equal<T>(rows[i][width],T(0)))return std::nullopt;
 LinearSystemSolution<T> answer;answer.particular.resize(width);
 for(Int r=Int(pivots.size())-1;r>=0;--r){Int col=pivots[r];T value=rows[r][width];for(Int j=col+1;j<width;++j)value=field_sub(value,field_mul<T>(rows[r][j],answer.particular[j]));answer.particular[col]=value;}
 std::vector<bool> is_pivot(width);for(Int col:pivots)is_pivot[col]=true;
 for(Int free=0;free<width;++free){if(is_pivot[free])continue;std::vector<T> vector(width);vector[free]=T(1);for(Int r=Int(pivots.size())-1;r>=0;--r){Int col=pivots[r];for(Int j=col+1;j<width;++j)vector[col]=field_sub<T>(vector[col],field_mul<T>(rows[r][j],vector[j]));}answer.basis.push_back(std::move(vector));}
 return answer;
}
// 変換行列から逆行列・余因子行列を求め、特異な場合もO(n³)で処理。
template<class T> std::optional<std::vector<std::vector<T>>> fieldAdjugateInverse(std::vector<std::vector<T>> a,bool adjugate){
 using namespace detail;Int n=a.size();for(Int i=0;i<n;++i){a[i].resize(2*n);a[i][n+i]=T(1);}auto e=field_echelon(a,n);Int rank=e.pivots.size();if(!adjugate&&rank!=n)return std::nullopt;std::vector<std::vector<T>> answer(n,std::vector<T>(n));if(rank<n-1)return answer;
 for(Int r=rank-1;r>=0;--r){Int col=e.pivots[r];for(Int i=0;i<r;++i){T factor=a[i][col];if(field_equal(factor,T(0)))continue;a[i][col]=T(0);for(Int j=col+1;j<2*n;++j)a[i][j]=field_sub<T>(a[i][j],field_mul<T>(factor,a[r][j]));}}
 if(rank==n){T scale=adjugate?e.determinant:T(1);for(Int i=0;i<n;++i)for(Int j=0;j<n;++j)answer[i][j]=field_mul<T>(scale,a[i][n+j]);}
 else{Int free=0;for(Int col:e.pivots)if(col==free)++free;T scale=e.determinant;if((n-1-free)&1)scale=field_sub(T(0),scale);for(Int j=0;j<n;++j){answer[free][j]=field_mul<T>(scale,a[n-1][n+j]);for(Int r=0;r<rank;++r)answer[e.pivots[r]][j]=field_sub(T(0),field_mul<T>(a[r][free],answer[free][j]));}}
 return answer;
}
// 多項式を用いる包除原理。O(n² 2^(n/2))時間、O(n⁴)領域。
template<class T> T fieldHafnian(const std::vector<std::vector<T>>& rows){
 Int n=rows.size();assert(n%2==0);for(Int i=0;i<n;++i)for(Int j=0;j<i;++j)assert(detail::field_equal<T>(rows[i][j],rows[j][i]));
 if constexpr(std::is_same_v<T,bool>){auto a=rows;for(Int i=0;i<n;++i)a[i][i]=false;return fieldDeterminant(std::move(a));}
 else{
  Int degree=n/2,stride=degree+1;
  auto add_product=[&](std::vector<T>& target,Int offset,const std::vector<T>& a,Int x,Int y){for(Int i=0;i<degree;++i){if(detail::field_equal(a[x+i],T(0)))continue;for(Int j=0;j<degree-i;++j)target[offset+i+j+1]+=a[x+i]*a[y+j];}};
  auto solve=[&](auto&& self,const std::vector<T>& a,Int size)->std::vector<T>{std::vector<T> result(stride);if(!size){result[0]=T(1);return result;}Int m=size-2,end=m*(m-1)/2*stride;std::vector<T> reduced(a.begin(),a.begin()+end);auto without=self(self,reduced,m);Int u=m*(m-1)/2*stride,v=m*(m+1)/2*stride;for(Int i=0;i<m;++i)for(Int j=0;j<i;++j){Int offset=(i*(i-1)/2+j)*stride;add_product(reduced,offset,a,u+i*stride,v+j*stride);add_product(reduced,offset,a,v+i*stride,u+j*stride);}auto with_pair=self(self,reduced,m);for(Int i=0;i<=degree;++i)result[i]=with_pair[i]-without[i];for(Int i=0;i<degree;++i)for(Int j=0;j<degree-i;++j)result[i+j+1]+=with_pair[i]*a[v+m*stride+j];return result;};
  std::vector<T> a(n*(n-1)/2*stride);for(Int i=0;i<n;++i)for(Int j=0;j<i;++j)a[(i*(i-1)/2+j)*stride]=rows[i][j];return solve(solve,a,n)[degree];
 }
}
}
