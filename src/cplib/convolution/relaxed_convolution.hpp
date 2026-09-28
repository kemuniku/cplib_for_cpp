#pragma once
#include <cplib/convolution/convolution.hpp>
#include <cplib/convolution/ntt.hpp>
#include <optional>
namespace cplib {
template<class T> class RelaxedConvolution {
    Int coefficientCount,currentIndex=0;std::vector<T> left,right,product;std::vector<std::vector<T>> leftPrefixTransforms,rightPrefixTransforms;
    void addProduct(Int first,const std::vector<T>& values,Int offset,Int count){if(offset>=Int(values.size()))return;for(Int i=0;i<std::min(count,Int(values.size())-offset);++i)product[first+i]+=values[offset+i];}
    void ensurePrefixTransforms(Int level,Int blockSize){if(Int(leftPrefixTransforms.size())<=level){leftPrefixTransforms.resize(level+1);rightPrefixTransforms.resize(level+1);}if(!leftPrefixTransforms[level].empty())return;Int size=blockSize*2;leftPrefixTransforms[level]=std::vector<T>(left.begin(),left.begin()+size);rightPrefixTransforms[level]=std::vector<T>(right.begin(),right.begin()+size);ntt(leftPrefixTransforms[level]);ntt(rightPrefixTransforms[level]);}
public:
    explicit RelaxedConvolution(Int count=0):coefficientCount(count),left(count),right(count),product(count){assert(count>=0);}
    Int len()const{return currentIndex;}
    T pendingCoefficient()const{assert(currentIndex<coefficientCount);return product[currentIndex];}
    std::vector<T> coefficients()const{return {product.begin(),product.begin()+currentIndex};}
    // 二入力のオンライン積。接頭辞NTTを再利用する元実装。全N項O(N log²N)。
    T add(T leftValue,T rightValue){assert(currentIndex<coefficientCount);Int index=currentIndex;left[index]=leftValue;right[index]=rightValue;product[index]+=leftValue*right[0];if(index>0)product[index]+=rightValue*left[0];if(++currentIndex>=coefficientCount)return product[index];Int level=0;for(Int blockSize=1;blockSize<=currentIndex;blockSize*=2,++level){if(currentIndex%(blockSize*2)!=blockSize)continue;Int count=std::min(blockSize,coefficientCount-currentIndex),size=blockSize*2;bool available=size>=64&&(T::umod()-1u)%std::uint32_t(size)==0;
            if(blockSize<=16){for(Int offset=0;offset<count;++offset){Int degree=blockSize+offset;if(currentIndex==blockSize){for(Int i=0;i<blockSize;++i){Int j=degree-i;if(j>=0&&j<blockSize)product[currentIndex+offset]+=left[i]*right[j];}}else for(Int i=0;i<blockSize;++i){Int j=degree-i;product[currentIndex+offset]+=left[currentIndex-blockSize+i]*right[j]+right[currentIndex-blockSize+i]*left[j];}}}
            else if(available){std::vector<T> a(size),b(size);if(currentIndex==blockSize){std::copy_n(left.begin(),blockSize,a.begin());std::copy_n(right.begin(),blockSize,b.begin());ntt(a);ntt(b);for(Int i=0;i<size;++i)a[i]*=b[i];}else{ensurePrefixTransforms(level,blockSize);std::copy_n(left.begin()+currentIndex-blockSize,blockSize,a.begin());std::copy_n(right.begin()+currentIndex-blockSize,blockSize,b.begin());ntt(a);ntt(b);for(Int i=0;i<size;++i)a[i]=a[i]*rightPrefixTransforms[level][i]+b[i]*leftPrefixTransforms[level][i];}intt(a);addProduct(currentIndex,a,blockSize,count);}
            else if(currentIndex==blockSize){auto values=convolution(std::vector<T>(left.begin(),left.begin()+blockSize),std::vector<T>(right.begin(),right.begin()+blockSize));addProduct(currentIndex,values,blockSize,count);}
            else{auto a=convolutionCyclicPowerOfTwo(std::vector<T>(left.begin()+currentIndex-blockSize,left.begin()+currentIndex),std::vector<T>(right.begin(),right.begin()+size),size);auto b=convolutionCyclicPowerOfTwo(std::vector<T>(right.begin()+currentIndex-blockSize,right.begin()+currentIndex),std::vector<T>(left.begin(),left.begin()+size),size);for(Int i=0;i<count;++i)product[currentIndex+i]+=a[blockSize+i]+b[blockSize+i];}break;
        }return product[index];}
    T append(T a,T b){return add(a,b);}T get(T a,T b){return add(a,b);}
};
namespace detail {
// 素数法のTonelli–Shanks。返す根の選び方も元実装と同一。
template<class T> std::optional<T> relaxedModSqrt(T a){Int p=T::umod();if(a.val()==0)return T(0);if(p==2)return a;if(a.pow((p-1)/2).val()!=1)return {};if(p%4==3)return a.pow((p+1)/4);Int q=p-1,s=0;while(!(q&1)){q>>=1;++s;}T z=2;while(z.pow((p-1)/2).val()!=p-1)z+=T(1);T c=z.pow(q),x=a.pow((q+1)/2),t=a.pow(q);Int m=s;while(t.val()!=1){Int i=1;T squared=t*t;while(i<m&&squared.val()!=1){squared*=squared;++i;}if(i==m)return {};T b=c.pow(Int(1)<<(m-i-1));x*=b;c=b*b;t*=c;m=i;}return x;}
template<class Derived,class T> struct RelaxedSingleAliases {auto append(T x){return static_cast<Derived*>(this)->add(x);}auto get(T x){return static_cast<Derived*>(this)->add(x);}};
}
template<class T> class RelaxedInv:public detail::RelaxedSingleAliases<RelaxedInv<T>,T> {
    template<class> friend class RelaxedLog;
    Int coefficientCount,currentIndex=0;RelaxedConvolution<T> convolution;std::vector<T> values;T constantInverse{};
public:
    explicit RelaxedInv(Int count=0):coefficientCount(count),convolution(count){assert(count>=0);values.reserve(count);}
    Int len()const{return currentIndex;}std::vector<T> coefficients()const{return values;}
    T add(T coefficient){assert(currentIndex<coefficientCount);T result;if(currentIndex==0){assert(coefficient.val()!=0);result=constantInverse=coefficient.inv();}else result=-(convolution.pendingCoefficient()+coefficient*values[0])*constantInverse;convolution.add(coefficient,result);values.push_back(result);++currentIndex;return result;}
};
template<class T> class RelaxedExp:public detail::RelaxedSingleAliases<RelaxedExp<T>,T> {
    Int coefficientCount,currentIndex=0;RelaxedConvolution<T> convolution;std::vector<T> values;
public:
    explicit RelaxedExp(Int count=0):coefficientCount(count),convolution(std::max(count-1,Int(0))){assert(count>=0&&count<=T::umod());values.reserve(count);}
    Int len()const{return currentIndex;}std::vector<T> coefficients()const{return values;}
    T add(T coefficient){assert(currentIndex<coefficientCount);T result;if(currentIndex==0){assert(coefficient.val()==0);result=T(1);}else result=convolution.add(values[currentIndex-1],coefficient*T(currentIndex))/T(currentIndex);values.push_back(result);++currentIndex;return result;}
};
template<class T> class RelaxedLog:public detail::RelaxedSingleAliases<RelaxedLog<T>,T> {
    Int coefficientCount,currentIndex=0;RelaxedInv<T> inverse;RelaxedConvolution<T> convolution;std::vector<T> values;
public:
    explicit RelaxedLog(Int count=0):coefficientCount(count),inverse(count),convolution(std::max(count-1,Int(0))){assert(count>=0&&count<=T::umod());values.reserve(count);}
    Int len()const{return currentIndex;}std::vector<T> coefficients()const{return values;}
    T add(T coefficient){assert(currentIndex<coefficientCount);T result=0;if(currentIndex==0){assert(coefficient.val()==1);inverse.add(coefficient);}else{inverse.add(coefficient);result=convolution.add(coefficient*T(currentIndex),inverse.values[currentIndex-1])/T(currentIndex);}values.push_back(result);++currentIndex;return result;}
};
template<class T> class RelaxedSqrt:public detail::RelaxedSingleAliases<RelaxedSqrt<T>,T> {
    Int coefficientCount,currentIndex=0;RelaxedConvolution<T> convolution;std::vector<T> values;T inverseDoubleConstant{};bool failed=false;
public:
    explicit RelaxedSqrt(Int count=0):coefficientCount(count),convolution(count){assert(count>=0);values.reserve(count);}
    Int len()const{return currentIndex;}std::optional<std::vector<T>> coefficients()const{if(failed)return {};return values;}
    // 逐次版は非零定数項のみ対応。根が存在しなければ以後もnullopt。
    std::optional<T> add(T coefficient){assert(currentIndex<coefficientCount);if(failed){++currentIndex;return {};}T value;if(currentIndex==0){assert(coefficient.val()!=0);auto root=detail::relaxedModSqrt(coefficient);if(!root){failed=true;++currentIndex;return {};}value=*root;inverseDoubleConstant=(value*T(2)).inv();}else value=(coefficient-convolution.pendingCoefficient())*inverseDoubleConstant;convolution.add(value,value);values.push_back(value);++currentIndex;return value;}
};
template<class T> class RelaxedPow:public detail::RelaxedSingleAliases<RelaxedPow<T>,T> {
    Int coefficientCount,currentIndex=0,exponent,order=-1,shift,bodyCount=0;T constantInverse{},constantPower{};RelaxedLog<T> logarithm;RelaxedExp<T> exponential;std::vector<T> body,values;
public:
    RelaxedPow(Int count,Int exponent):coefficientCount(count),exponent(exponent),shift(count){assert(count>=0&&count<=T::umod()&&exponent>=0);values.reserve(count);}
    Int len()const{return currentIndex;}std::vector<T> coefficients()const{return values;}
    T add(T coefficient){assert(currentIndex<coefficientCount);Int degree=currentIndex;T result=0;if(exponent==0)result=T(degree==0);else{if(order<0&&coefficient.val()!=0){order=degree;if(degree<=(coefficientCount-1)/exponent){shift=degree*exponent;bodyCount=coefficientCount-shift;constantInverse=coefficient.inv();constantPower=coefficient.pow(exponent);logarithm=RelaxedLog<T>(bodyCount);exponential=RelaxedExp<T>(bodyCount);body.push_back(exponential.add(logarithm.add(T(1))*T(exponent))*constantPower);}}else if(order>=0&&Int(body.size())<bodyCount)body.push_back(exponential.add(logarithm.add(coefficient*constantInverse)*T(exponent))*constantPower);if(degree>=shift)result=body[degree-shift];}values.push_back(result);++currentIndex;return result;}
};
template<class T> auto initRelaxedConvolution(Int count){return RelaxedConvolution<T>(count);}
template<class T> Int len(const RelaxedConvolution<T>& s){return s.len();}
template<class T> auto coefficients(const RelaxedConvolution<T>& s){return s.coefficients();}
template<class T> auto add(RelaxedConvolution<T>& s,T a,T b){return s.add(a,b);}
template<class T> auto append(RelaxedConvolution<T>& s,T a,T b){return s.add(a,b);}
template<class T> auto get(RelaxedConvolution<T>& s,T a,T b){return s.add(a,b);}
template<class T> auto initRelaxedInv(Int count){return RelaxedInv<T>(count);}
template<class T> Int len(const RelaxedInv<T>& s){return s.len();}
template<class T> auto coefficients(const RelaxedInv<T>& s){return s.coefficients();}
template<class T> auto add(RelaxedInv<T>& s,T a){return s.add(a);}
template<class T> auto append(RelaxedInv<T>& s,T a){return s.add(a);}
template<class T> auto get(RelaxedInv<T>& s,T a){return s.add(a);}
template<class T> auto initRelaxedExp(Int count){return RelaxedExp<T>(count);}
template<class T> Int len(const RelaxedExp<T>& s){return s.len();}
template<class T> auto coefficients(const RelaxedExp<T>& s){return s.coefficients();}
template<class T> auto add(RelaxedExp<T>& s,T a){return s.add(a);}
template<class T> auto append(RelaxedExp<T>& s,T a){return s.add(a);}
template<class T> auto get(RelaxedExp<T>& s,T a){return s.add(a);}
template<class T> auto initRelaxedLog(Int count){return RelaxedLog<T>(count);}
template<class T> Int len(const RelaxedLog<T>& s){return s.len();}
template<class T> auto coefficients(const RelaxedLog<T>& s){return s.coefficients();}
template<class T> auto add(RelaxedLog<T>& s,T a){return s.add(a);}
template<class T> auto append(RelaxedLog<T>& s,T a){return s.add(a);}
template<class T> auto get(RelaxedLog<T>& s,T a){return s.add(a);}
template<class T> auto initRelaxedSqrt(Int count){return RelaxedSqrt<T>(count);}
template<class T> Int len(const RelaxedSqrt<T>& s){return s.len();}
template<class T> auto coefficients(const RelaxedSqrt<T>& s){return s.coefficients();}
template<class T> auto add(RelaxedSqrt<T>& s,T a){return s.add(a);}
template<class T> auto append(RelaxedSqrt<T>& s,T a){return s.add(a);}
template<class T> auto get(RelaxedSqrt<T>& s,T a){return s.add(a);}
template<class T> auto initRelaxedPow(Int count,Int exponent){return RelaxedPow<T>(count,exponent);}
template<class T> Int len(const RelaxedPow<T>& s){return s.len();}
template<class T> auto coefficients(const RelaxedPow<T>& s){return s.coefficients();}
template<class T> auto add(RelaxedPow<T>& s,T a){return s.add(a);}
template<class T> auto append(RelaxedPow<T>& s,T a){return s.add(a);}
template<class T> auto get(RelaxedPow<T>& s,T a){return s.add(a);}
template<class T> T pendingCoefficient(const RelaxedConvolution<T>& s){return s.pendingCoefficient();}
template<class T> std::vector<T> invRelaxed(const std::vector<T>& f,Int n){if(n<=0)return {};RelaxedInv<T> state(n);std::vector<T> out(n);for(Int i=0;i<n;++i)out[i]=state.add(i<Int(f.size())?f[i]:T(0));return out;}
template<class T> auto invRelaxed(const std::vector<T>& f){return invRelaxed(f,Int(f.size()));}
template<class T> std::vector<T> expRelaxed(const std::vector<T>& f,Int n){if(n<=0)return {};RelaxedExp<T> state(n);std::vector<T> out(n);for(Int i=0;i<n;++i)out[i]=state.add(i<Int(f.size())?f[i]:T(0));return out;}
template<class T> auto expRelaxed(const std::vector<T>& f){return expRelaxed(f,Int(f.size()));}
template<class T> std::vector<T> logRelaxed(const std::vector<T>& f,Int n){if(n<=0)return {};RelaxedLog<T> state(n);std::vector<T> out(n);for(Int i=0;i<n;++i)out[i]=state.add(i<Int(f.size())?f[i]:T(0));return out;}
template<class T> auto logRelaxed(const std::vector<T>& f){return logRelaxed(f,Int(f.size()));}
template<class T> std::optional<std::vector<T>> sqrtRelaxed(const std::vector<T>& f,Int n){if(n<=0)return std::vector<T>{};Int order=0;while(order<std::min(Int(f.size()),n)&&f[order].val()==0)++order;if(order==std::min(Int(f.size()),n))return std::vector<T>(n);if(order&1)return {};Int shift=order/2,size=n-shift;RelaxedSqrt<T> state(size);std::vector<T> out(n);for(Int i=0;i<size;++i){auto root=state.add(order+i<Int(f.size())?f[order+i]:T(0));if(!root)return {};out[shift+i]=*root;}return out;}
template<class T> auto sqrtRelaxed(const std::vector<T>& f){return sqrtRelaxed(f,Int(f.size()));}
template<class T> std::vector<T> powRelaxed(const std::vector<T>& f,Int k,Int n){assert(k>=0);if(n<=0)return {};RelaxedPow<T> state(n,k);std::vector<T> out(n);for(Int i=0;i<n;++i)out[i]=state.add(i<Int(f.size())?f[i]:T(0));return out;}
template<class T> auto powRelaxed(const std::vector<T>& f,Int k){return powRelaxed(f,k,Int(f.size()));}
}
