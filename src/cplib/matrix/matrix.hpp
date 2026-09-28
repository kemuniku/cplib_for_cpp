#pragma once
#include <cplib/matrix/field_matrix_ops.hpp>
#include <sstream>
namespace cplib {
// 動的版は元と同じ行ごとの配列。0行の場合も列数を保持する。
template<class T> class Matrix {
 std::vector<std::vector<T>> arr_;Int empty_width_=0;
public:
 using value_type=T;
 Matrix()=default;
 Matrix(Int height,Int width,T value=T(0)):arr_(height,std::vector<T>(width,value)),empty_width_(width){assert(height>=0&&width>=0);}
 explicit Matrix(const std::vector<std::vector<T>>& rows):arr_(rows){for(const auto& row:rows)assert(row.size()==(rows.empty()?0:rows[0].size()));}
 explicit Matrix(const std::vector<T>& values,bool vertical=false){if(vertical){empty_width_=1;for(T value:values)arr_.push_back({value});}else arr_.push_back(values);}
 Int h()const{return arr_.size();}Int w()const{return arr_.empty()?empty_width_:Int(arr_[0].size());}
 auto& operator[](Int row){return arr_[row];}const auto& operator[](Int row)const{return arr_[row];}
 decltype(auto) operator()(Int row,Int col){return arr_[row][col];}decltype(auto) operator()(Int row,Int col)const{return arr_[row][col];}
 bool operator==(const Matrix& b)const{return h()==b.h()&&w()==b.w()&&arr_==b.arr_;}
 Matrix operator-()const{auto out=*this;for(Int i=0;i<h();++i)for(Int j=0;j<w();++j)out(i,j)=-(*this)(i,j);return out;}
 Matrix& operator*=(const Matrix& b){assert(w()==b.h());Matrix out(h(),b.w(),T(0));for(Int i=0;i<h();++i)for(Int j=0;j<b.w();++j)for(Int k=0;k<w();++k)out(i,j)+=(*this)(i,k)*b(k,j);*this=std::move(out);return *this;}
 Matrix& operator*=(T value){for(Int i=0;i<h();++i)for(Int j=0;j<w();++j)(*this)(i,j)*=value;return *this;}
 friend Matrix operator*(Matrix a,const Matrix& b){a*=b;return a;}friend Matrix operator*(Matrix a,T b){a*=b;return a;}friend Matrix operator*(T a,Matrix b){b*=a;return b;}
 Matrix& operator+=(const Matrix& b){assert(h()==b.h()&&w()==b.w());for(Int i=0;i<h();++i)for(Int j=0;j<w();++j){T a=(*this)(i,j);T rhs=b(i,j);(*this)(i,j)=a + rhs;}return *this;}
 Matrix& operator+=(T b){for(Int i=0;i<h();++i)for(Int j=0;j<w();++j){T a=(*this)(i,j);(*this)(i,j)=a + b;}return *this;}
 friend Matrix operator+(Matrix a,const Matrix& b){a+=b;return a;}
 friend Matrix operator+(Matrix a,T b){a+=b;return a;}
 friend Matrix operator+(T a,Matrix matrix){for(Int i=0;i<matrix.h();++i)for(Int j=0;j<matrix.w();++j){T b=matrix(i,j);matrix(i,j)=a + b;}return matrix;}
 Matrix& operator-=(const Matrix& b){assert(h()==b.h()&&w()==b.w());for(Int i=0;i<h();++i)for(Int j=0;j<w();++j){T a=(*this)(i,j);T rhs=b(i,j);(*this)(i,j)=a - rhs;}return *this;}
 Matrix& operator-=(T b){for(Int i=0;i<h();++i)for(Int j=0;j<w();++j){T a=(*this)(i,j);(*this)(i,j)=a - b;}return *this;}
 friend Matrix operator-(Matrix a,const Matrix& b){a-=b;return a;}
 friend Matrix operator-(Matrix a,T b){a-=b;return a;}
 friend Matrix operator-(T a,Matrix matrix){for(Int i=0;i<matrix.h();++i)for(Int j=0;j<matrix.w();++j){T b=matrix(i,j);matrix(i,j)=a - b;}return matrix;}
 Matrix& operator&=(const Matrix& b) requires(std::integral<T>){assert(h()==b.h()&&w()==b.w());for(Int i=0;i<h();++i)for(Int j=0;j<w();++j){T a=(*this)(i,j);T rhs=b(i,j);(*this)(i,j)=a & rhs;}return *this;}
 Matrix& operator&=(T b) requires(std::integral<T>){for(Int i=0;i<h();++i)for(Int j=0;j<w();++j){T a=(*this)(i,j);(*this)(i,j)=a & b;}return *this;}
 friend Matrix operator&(Matrix a,const Matrix& b) requires(std::integral<T>){a&=b;return a;}
 friend Matrix operator&(Matrix a,T b) requires(std::integral<T>){a&=b;return a;}
 friend Matrix operator&(T a,Matrix matrix) requires(std::integral<T>){for(Int i=0;i<matrix.h();++i)for(Int j=0;j<matrix.w();++j){T b=matrix(i,j);matrix(i,j)=a & b;}return matrix;}
 Matrix& operator|=(const Matrix& b) requires(std::integral<T>){assert(h()==b.h()&&w()==b.w());for(Int i=0;i<h();++i)for(Int j=0;j<w();++j){T a=(*this)(i,j);T rhs=b(i,j);(*this)(i,j)=a | rhs;}return *this;}
 Matrix& operator|=(T b) requires(std::integral<T>){for(Int i=0;i<h();++i)for(Int j=0;j<w();++j){T a=(*this)(i,j);(*this)(i,j)=a | b;}return *this;}
 friend Matrix operator|(Matrix a,const Matrix& b) requires(std::integral<T>){a|=b;return a;}
 friend Matrix operator|(Matrix a,T b) requires(std::integral<T>){a|=b;return a;}
 friend Matrix operator|(T a,Matrix matrix) requires(std::integral<T>){for(Int i=0;i<matrix.h();++i)for(Int j=0;j<matrix.w();++j){T b=matrix(i,j);matrix(i,j)=a | b;}return matrix;}
 Matrix& operator^=(const Matrix& b) requires(std::integral<T>){assert(h()==b.h()&&w()==b.w());for(Int i=0;i<h();++i)for(Int j=0;j<w();++j){T a=(*this)(i,j);T rhs=b(i,j);(*this)(i,j)=a ^ rhs;}return *this;}
 Matrix& operator^=(T b) requires(std::integral<T>){for(Int i=0;i<h();++i)for(Int j=0;j<w();++j){T a=(*this)(i,j);(*this)(i,j)=a ^ b;}return *this;}
 friend Matrix operator^(Matrix a,const Matrix& b) requires(std::integral<T>){a^=b;return a;}
 friend Matrix operator^(Matrix a,T b) requires(std::integral<T>){a^=b;return a;}
 friend Matrix operator^(T a,Matrix matrix) requires(std::integral<T>){for(Int i=0;i<matrix.h();++i)for(Int j=0;j<matrix.w();++j){T b=matrix(i,j);matrix(i,j)=a ^ b;}return matrix;}
 Matrix& operator<<=(const Matrix& b) requires(std::integral<T>){assert(h()==b.h()&&w()==b.w());for(Int i=0;i<h();++i)for(Int j=0;j<w();++j){T a=(*this)(i,j);T rhs=b(i,j);(*this)(i,j)=T(std::make_unsigned_t<T>(a) << rhs);}return *this;}
 Matrix& operator<<=(T b) requires(std::integral<T>){for(Int i=0;i<h();++i)for(Int j=0;j<w();++j){T a=(*this)(i,j);(*this)(i,j)=T(std::make_unsigned_t<T>(a) << b);}return *this;}
 friend Matrix operator<<(Matrix a,const Matrix& b) requires(std::integral<T>){a<<=b;return a;}
 friend Matrix operator<<(Matrix a,T b) requires(std::integral<T>){a<<=b;return a;}
 friend Matrix operator<<(T a,Matrix matrix) requires(std::integral<T>){for(Int i=0;i<matrix.h();++i)for(Int j=0;j<matrix.w();++j){T b=matrix(i,j);matrix(i,j)=T(std::make_unsigned_t<T>(a) << b);}return matrix;}
 Matrix& operator>>=(const Matrix& b) requires(std::integral<T>){assert(h()==b.h()&&w()==b.w());for(Int i=0;i<h();++i)for(Int j=0;j<w();++j){T a=(*this)(i,j);T rhs=b(i,j);(*this)(i,j)=a >> rhs;}return *this;}
 Matrix& operator>>=(T b) requires(std::integral<T>){for(Int i=0;i<h();++i)for(Int j=0;j<w();++j){T a=(*this)(i,j);(*this)(i,j)=a >> b;}return *this;}
 friend Matrix operator>>(Matrix a,const Matrix& b) requires(std::integral<T>){a>>=b;return a;}
 friend Matrix operator>>(Matrix a,T b) requires(std::integral<T>){a>>=b;return a;}
 friend Matrix operator>>(T a,Matrix matrix) requires(std::integral<T>){for(Int i=0;i<matrix.h();++i)for(Int j=0;j<matrix.w();++j){T b=matrix(i,j);matrix(i,j)=a >> b;}return matrix;}
 Matrix& operator/=(const Matrix& b) requires(std::integral<T>){assert(h()==b.h()&&w()==b.w());for(Int i=0;i<h();++i)for(Int j=0;j<w();++j){T a=(*this)(i,j);T rhs=b(i,j);(*this)(i,j)=a / rhs;}return *this;}
 Matrix& operator/=(T b) requires(std::integral<T>){for(Int i=0;i<h();++i)for(Int j=0;j<w();++j){T a=(*this)(i,j);(*this)(i,j)=a / b;}return *this;}
 friend Matrix operator/(Matrix a,const Matrix& b) requires(std::integral<T>){a/=b;return a;}
 friend Matrix operator/(Matrix a,T b) requires(std::integral<T>){a/=b;return a;}
 friend Matrix operator/(T a,Matrix matrix) requires(std::integral<T>){for(Int i=0;i<matrix.h();++i)for(Int j=0;j<matrix.w();++j){T b=matrix(i,j);matrix(i,j)=a / b;}return matrix;}
 Matrix& operator%=(const Matrix& b) requires(std::integral<T>){assert(h()==b.h()&&w()==b.w());for(Int i=0;i<h();++i)for(Int j=0;j<w();++j){T a=(*this)(i,j);T rhs=b(i,j);(*this)(i,j)=a % rhs;}return *this;}
 Matrix& operator%=(T b) requires(std::integral<T>){for(Int i=0;i<h();++i)for(Int j=0;j<w();++j){T a=(*this)(i,j);(*this)(i,j)=a % b;}return *this;}
 friend Matrix operator%(Matrix a,const Matrix& b) requires(std::integral<T>){a%=b;return a;}
 friend Matrix operator%(Matrix a,T b) requires(std::integral<T>){a%=b;return a;}
 friend Matrix operator%(T a,Matrix matrix) requires(std::integral<T>){for(Int i=0;i<matrix.h();++i)for(Int j=0;j<matrix.w();++j){T b=matrix(i,j);matrix(i,j)=a % b;}return matrix;}
 static Matrix identity(Int n,T one=T(1),T zero=T(0)){Matrix out(n,n,zero);for(Int i=0;i<n;++i)out(i,i)=one;return out;}
 Matrix pow(Int exponent)const{Matrix result=identity(h()),base=*this;while(exponent>0){if(exponent&1)result*=base;base*=base;exponent>>=1;}return result;}
 T sum()const{T result=T(0);for(const auto& row:arr_)for(T value:row)result+=value;return result;}
 std::string str()const{std::ostringstream out;for(Int i=0;i<h();++i){if(i)out<<'\n';for(Int j=0;j<w();++j){if(j)out<<' ';if constexpr(std::is_same_v<T,bool>)out<<((*this)(i,j)?"true":"false");else out<<(*this)(i,j);}}return out.str();}
 std::size_t hash()const{std::size_t seed=arr_.size();for(const auto& row:arr_){seed^=row.size()+0x9e3779b9+(seed<<6)+(seed>>2);for(T value:row)seed^=std::hash<T>{}(value)+0x9e3779b9+(seed<<6)+(seed>>2);}return seed;}
 Int rank()const{return fieldRank(matrixRows(*this,h(),w()),w());}
 T determinant()const{assert(h()==w());return fieldDeterminant(matrixRows(*this,h(),w()));}
 T hafnian()const{assert(h()==w());return fieldHafnian(matrixRows(*this,h(),w()));}
 template<class B> auto solveLinearSystem(const B& b)const{return fieldSolve(matrixRows(*this,h(),w()),w(),b);}
 std::optional<Matrix> inverse()const{assert(h()==w());auto rows=fieldAdjugateInverse(matrixRows(*this,h(),w()),false);if(!rows)return std::nullopt;return Matrix(*rows);}
 Matrix adjugate()const{assert(h()==w());return Matrix(*fieldAdjugateInverse(matrixRows(*this,h(),w()),true));}
};
template<class T> auto initMatrix(const std::vector<std::vector<T>>& rows){return Matrix<T>(rows);}
template<class T> auto toMatrix(const std::vector<std::vector<T>>& rows){return Matrix<T>(rows);}
template<class T> auto initMatrix(const std::vector<T>& values,bool vertical=false){return Matrix<T>(values,vertical);}
template<class T> auto initMatrix(Int h,Int w,T value){return Matrix<T>(h,w,value);}
template<class T> auto identity_matrix(Int n,T one=T(1),T zero=T(0)){return Matrix<T>::identity(n,one,zero);}
template<class T> Int h(const Matrix<T>& m){return m.h();}template<class T> Int w(const Matrix<T>& m){return m.w();}
template<class T> auto to_string(const Matrix<T>& m){return m.str();}template<class T> auto hash(const Matrix<T>& m){return m.hash();}
template<class T> auto pow(const Matrix<T>& m,Int n){return m.pow(n);}template<class T> T sum(const Matrix<T>& m){return m.sum();}
template<class T> Int rank(const Matrix<T>& m){return m.rank();}template<class T> T determinant(const Matrix<T>& m){return m.determinant();}template<class T> T hafnian(const Matrix<T>& m){return m.hafnian();}
template<class T,class B> auto solveLinearSystem(const Matrix<T>& m,const B& b){return m.solveLinearSystem(b);}
template<class T> auto inverse(const Matrix<T>& m){return m.inverse();}template<class T> auto adjugate(const Matrix<T>& m){return m.adjugate();}
}
namespace std {template<class T> struct hash<cplib::Matrix<T>>{size_t operator()(const cplib::Matrix<T>& m)const{return m.hash();}};}
