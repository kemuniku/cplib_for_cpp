#include <cplib/matrix/matrix_product_avx2.hpp>
#include <random>
using namespace cplib;
template<class A,class B> void equal(const A& a,const B& b){assert(a.h()==b.h()&&a.w()==b.w());for(Int i=0;i<a.h();++i)for(Int j=0;j<a.w();++j)assert(a(i,j)==b(i,j));}
template<class T> void products(){std::mt19937 rng(11291);for(auto shape:std::vector<std::array<Int,3>>{{0,0,5},{4,0,9},{0,2,0},{1,1,1},{3,7,9},{9,17,15},{39,63,65},{64,64,64},{65,65,65},{80,96,112},{129,131,137},{257,259,261}}){auto [h,w,k]=shape;Matrix<T> a(h,w),b(w,k);Matrix<T> x(h,w),y(w,k);std::vector<std::uint32_t> aa(h*w),bb(w*k);for(Int i=0;i<h;++i)for(Int j=0;j<w;++j){aa[i*w+j]=rng()%T::umod();a(i,j)=x(i,j)=T(aa[i*w+j]);}for(Int i=0;i<w;++i)for(Int j=0;j<k;++j){bb[i*k+j]=rng()%T::umod();b(i,j)=y(i,j)=T(bb[i*k+j]);}auto c=a*b;equal(x*y,c);equal(cplib::matrixProduct(a,b),c);equal(cplib::matrixProduct(a,b),c);auto old=cplib::matrixProduct(aa,bb,h,w,k,T::umod()),cur=cplib::matrixProduct(aa,bb,h,w,k,T::umod());assert(old==cur);for(Int i=0;i<h;++i)for(Int j=0;j<k;++j)assert(cur[i*k+j]==c(i,j).val());}}
int main(){products<StaticMontgomeryModint<998244353>>();products<StaticBarrettModint<1000000007>>();}
