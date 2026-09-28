#pragma once
#include <cplib/matrix/matrix_avx2_common.hpp>
#include <array>
namespace cplib {
namespace impl=::cplib::detail::matrix_avx2;
template<Int H,Int W,impl::Element T> class StaticMatrix {
 static_assert(H>=0&&W>=0&&H<=INT32_MAX&&W<=INT32_MAX);static_assert(H==0||W<=(std::numeric_limits<Int>::max()/4)/H);
 std::array<T,H*W> values_{};std::uint32_t modulus_=0;
 template<Int,Int,impl::Element> friend class StaticMatrix;
 void check()const{[[maybe_unused]]auto p=impl::modulus<T>();assert(!modulus_||modulus_==p);}
 template<Int E> struct Reduction {std::array<std::uint32_t,H*(W+E)> values{};std::array<int,std::min(H,W)> pivots{};Int width=0,rank=0;std::uint32_t determinant=0;};
 template<Int E> auto reduce(Int height,Int width,bool reduced,const T* rhs=nullptr,bool identity=false)const{
  static_assert(E>=0&&W<=INT32_MAX-E);check();assert(height>=0&&height<=H&&width>=0&&width<=W);auto r=std::make_unique<Reduction<E>>();Int extra=identity?height:E;r->width=width+extra;
  fieldPrepareKernel(impl::raw(values_.data()),impl::raw(rhs),r->values.data(),height,width,extra,impl::modulus<T>(),MontgomeryModint<T>,identity,W);
  r->rank=fieldEliminateKernel(r->values.data(),height,r->width,width,r->pivots.data(),r->determinant,impl::modulus<T>(),reduced);return r;
 }
public:
 using value_type=T;StaticMatrix()=default;explicit StaticMatrix(T value):modulus_(impl::modulus<T>()){values_.fill(value);}
 explicit StaticMatrix(const std::array<std::array<T,W>,H>& rows):modulus_(impl::modulus<T>()){for(Int i=0;i<H;++i)std::copy(rows[i].begin(),rows[i].end(),values_.begin()+i*W);}
 static constexpr Int h(){return H;}static constexpr Int w(){return W;}
 const T& operator()(Int i,Int j)const{assert(i>=0&&i<H&&j>=0&&j<W);return values_[i*W+j];}T& operator()(Int i,Int j){assert(i>=0&&i<H&&j>=0&&j<W);return values_[i*W+j];}
 std::array<T,W> operator[](Int i)const{assert(i>=0&&i<H);std::array<T,W> out{};std::copy_n(values_.begin()+i*W,W,out.begin());return out;}
 struct Row {StaticMatrix* a;Int i;operator std::array<T,W>()const{return std::as_const(*a)[i];}T& operator[](Int j)const{return (*a)(i,j);}Row& operator=(const std::array<T,W>& row){a->setRow(i,row);return *this;}Row& operator=(const Row& b){return *this=std::array<T,W>(b);}};
 Row operator[](Int i){assert(i>=0&&i<H);return {this,i};}void setRow(Int i,const std::array<T,W>& row){assert(i>=0&&i<H);std::copy(row.begin(),row.end(),values_.begin()+i*W);}
 bool operator==(const StaticMatrix& b)const{check();b.check();for(Int i=0;i<H*W;++i)if(values_[i].val()!=b.values_[i].val())return false;return true;}
 std::string str()const{check();std::string out;for(Int i=0;i<H;++i){if(i)out+='\n';if constexpr(W>0)out+=matrixJoinValues(impl::raw(values_.data()+i*W),W,impl::modulus<T>(),MontgomeryModint<T>," ");}return out;}
 template<Int K> StaticMatrix<H,K,T> operator*(const StaticMatrix<W,K,T>& b)const{
  check();b.check();StaticMatrix<H,K,T> out(T(0));
  if constexpr(H<=4&&W<=4&&K<=4){for(Int i=0;i<H;++i)for(Int j=0;j<K;++j)for(Int k=0;k<W;++k)out(i,j)+=(*this)(i,k)*b(k,j);}
  else if constexpr(H>0&&W>0&&K>0){if constexpr(MontgomeryModint<T>)matrixProductMontgomeryKernel(impl::raw(values_.data()),impl::raw(b.values_.data()),impl::raw(out.values_.data()),H,W,K,impl::modulus<T>());else matrixProductKernel(impl::raw(values_.data()),impl::raw(b.values_.data()),impl::raw(out.values_.data()),H,W,K,impl::modulus<T>());}
  return out;
 }
 StaticMatrix& operator*=(const StaticMatrix<W,W,T>& b){return *this=*this*b;}
 StaticMatrix& operator+=(const StaticMatrix& b){check();b.check();for(Int i=0;i<H*W;++i)values_[i]+=b.values_[i];return *this;}
 friend StaticMatrix operator+(StaticMatrix a,const StaticMatrix& b){return a+=b;}
 StaticMatrix& operator+=(T value){check();for(auto& x:values_)x+=value;return *this;}
 friend StaticMatrix operator+(StaticMatrix a,T value){return a+=value;}
 StaticMatrix& operator-=(const StaticMatrix& b){check();b.check();for(Int i=0;i<H*W;++i)values_[i]-=b.values_[i];return *this;}
 friend StaticMatrix operator-(StaticMatrix a,const StaticMatrix& b){return a-=b;}
 StaticMatrix& operator-=(T value){check();for(auto& x:values_)x-=value;return *this;}
 friend StaticMatrix operator-(StaticMatrix a,T value){return a-=value;}
 StaticMatrix& operator*=(T value){check();for(auto& x:values_)x*=value;return *this;}
 friend StaticMatrix operator*(StaticMatrix a,T value){return a*=value;}friend StaticMatrix operator*(T value,StaticMatrix a){return a*=value;}
 StaticMatrix operator-()const{return *this*T(-1);}
 static StaticMatrix identity(Int n=H){static_assert(H==W);assert(n==H);StaticMatrix out(T(0));for(Int i=0;i<H;++i)out(i,i)=1;return out;}
 StaticMatrix pow(Int exponent)const{static_assert(H==W);assert(exponent>=0);check();StaticMatrix out=identity(),base=*this;while(exponent){if(exponent&1)out*=base;exponent>>=1;if(exponent)base*=base;}return out;}
 T sum()const{check();T out=0;for(auto x:values_)out+=x;return out;}
 Int rank(Int height=H,Int width=W)const{return reduce<0>(height,width,false)->rank;}
 T determinant(Int n=H)const{auto r=reduce<0>(n,n,false);return r->rank!=n?T(0):T(fieldCanonicalKernel(r->determinant,impl::modulus<T>()));}
 T hafnian(Int n=H)const{check();assert(n>=0&&n<=std::min(H,W)&&n%2==0);for(Int i=0;i<n;++i)for(Int j=0;j<i;++j)assert((*this)(i,j).val()==(*this)(j,i).val());return T(fieldHafnianKernel(impl::raw(values_.data()),n,impl::modulus<T>(),MontgomeryModint<T>,W));}
 std::optional<LinearSystemSolution<T>> solveLinearSystem(std::span<const T> rhs,Int height=H,Int width=W)const{assert(Int(rhs.size())==height);auto r=reduce<1>(height,width,true,rhs.data());return impl::solution<T>(*r,height,width);}
 std::optional<StaticMatrix> inverse(Int n=H)const{auto r=reduce<std::min(H,W)>(n,n,true,nullptr,true);if(r->rank!=n)return std::nullopt;StaticMatrix out(T(0));fieldInverseAdjugateKernel(r->values.data(),impl::raw(out.values_.data()),n,r->rank,r->pivots.data(),r->determinant,impl::modulus<T>(),MontgomeryModint<T>,false,W);return out;}
 StaticMatrix adjugate(Int n=H)const{auto r=reduce<std::min(H,W)>(n,n,true,nullptr,true);StaticMatrix out(T(0));fieldInverseAdjugateKernel(r->values.data(),impl::raw(out.values_.data()),n,r->rank,r->pivots.data(),r->determinant,impl::modulus<T>(),MontgomeryModint<T>,true,W);return out;}
};
template<Int H,Int W,class T> auto initMatrix(T value=T(0)){return StaticMatrix<H,W,T>(value);}
template<std::size_t H,std::size_t W,class T> auto toMatrix(const std::array<std::array<T,W>,H>& rows){return StaticMatrix<H,W,T>(rows);}
template<Int H,Int W,class T> auto identity_matrix(Int n=H){return StaticMatrix<H,W,T>::identity(n);}
template<Int H,Int W,Int K,class T> auto matrixProduct(const StaticMatrix<H,W,T>& a,const StaticMatrix<W,K,T>& b){return a*b;}
template<Int H,Int W,class T,class... A> auto h(const StaticMatrix<H,W,T>& a,A&&... args){return a.h(std::forward<A>(args)...);}
template<Int H,Int W,class T,class... A> auto w(const StaticMatrix<H,W,T>& a,A&&... args){return a.w(std::forward<A>(args)...);}
template<Int H,Int W,class T,class... A> auto rank(const StaticMatrix<H,W,T>& a,A&&... args){return a.rank(std::forward<A>(args)...);}
template<Int H,Int W,class T,class... A> auto determinant(const StaticMatrix<H,W,T>& a,A&&... args){return a.determinant(std::forward<A>(args)...);}
template<Int H,Int W,class T,class... A> auto hafnian(const StaticMatrix<H,W,T>& a,A&&... args){return a.hafnian(std::forward<A>(args)...);}
template<Int H,Int W,class T,class... A> auto inverse(const StaticMatrix<H,W,T>& a,A&&... args){return a.inverse(std::forward<A>(args)...);}
template<Int H,Int W,class T,class... A> auto adjugate(const StaticMatrix<H,W,T>& a,A&&... args){return a.adjugate(std::forward<A>(args)...);}
template<Int H,Int W,class T,class... A> auto sum(const StaticMatrix<H,W,T>& a,A&&... args){return a.sum(std::forward<A>(args)...);}
template<Int H,Int W,class T,class... A> auto pow(const StaticMatrix<H,W,T>& a,A&&... args){return a.pow(std::forward<A>(args)...);}
template<Int H,Int W,class T,class... A> auto solveLinearSystem(const StaticMatrix<H,W,T>& a,A&&... args){return a.solveLinearSystem(std::forward<A>(args)...);}
}
