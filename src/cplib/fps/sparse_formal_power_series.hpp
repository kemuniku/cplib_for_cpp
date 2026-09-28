#pragma once
#include <cplib/fps/bostan_mori.hpp>
#include <cplib/fps/product_tree.hpp>
#include <cplib/math/isqrt.hpp>
#include <memory>
#include <functional>
namespace cplib {
template<class T> struct SparseTerm {Int degree;T coefficient;bool operator==(const SparseTerm&)const=default;};
template<class T> using SparseFPS=std::vector<SparseTerm<T>>;
// 重複次数を統合して零項を除く。K項に対してO(K log K)。
template<class T> SparseFPS<T> initSparseFPS(SparseFPS<T> terms){for(auto term:terms)assert(term.degree>=0);std::sort(terms.begin(),terms.end(),[](auto a,auto b){return a.degree<b.degree;});SparseFPS<T> out;for(auto term:terms){if(term.coefficient.val()==0)continue;if(!out.empty()&&out.back().degree==term.degree){out.back().coefficient+=term.coefficient;if(out.back().coefficient.val()==0)out.pop_back();}else out.push_back(term);}return out;}
template<class T> SparseFPS<T> initSparseFPS(std::span<const SparseTerm<T>> terms){return initSparseFPS<T>(SparseFPS<T>(terms.begin(),terms.end()));}
template<class T> SparseFPS<T> initSparseFPS(std::initializer_list<SparseTerm<T>> terms){return initSparseFPS<T>(SparseFPS<T>(terms));}
// C++用の疎多項式リテラル。式木を各演算O(1)で構築し、最後にO(K)で列挙する。
// sfps<M>([](auto x){ return 1 + 2*x - (x^3); })。C++の^は優先順位が低いので括弧が必要。
template<Modint T> class SparseLiteral {
    struct Node {Int degree;T coefficient;std::shared_ptr<const Node> left,right;bool subtract=false;};
    std::shared_ptr<const Node> node;
    explicit SparseLiteral(std::shared_ptr<const Node> node):node(std::move(node)){}
public:
    SparseLiteral(T c):node(std::make_shared<Node>(Node{0,c,{},{},false})){}
    template<std::integral I> SparseLiteral(I c):SparseLiteral(T(c)){}
    static SparseLiteral variable(){return SparseLiteral(std::make_shared<Node>(Node{1,T(1),{},{},false}));}
    friend SparseLiteral operator+(SparseLiteral a,SparseLiteral b){return SparseLiteral(std::make_shared<Node>(Node{0,T(0),a.node,b.node,false}));}
    friend SparseLiteral operator-(SparseLiteral a,SparseLiteral b){return SparseLiteral(std::make_shared<Node>(Node{0,T(0),a.node,b.node,true}));}
    friend SparseLiteral operator-(SparseLiteral a){return SparseLiteral(T(0))-a;}
    friend SparseLiteral operator*(T c,SparseLiteral a){assert(!a.node->left);return SparseLiteral(std::make_shared<Node>(Node{a.node->degree,c*a.node->coefficient,{},{},false}));}
    friend SparseLiteral operator*(SparseLiteral a,T c){return c*a;}
    friend SparseLiteral operator^(SparseLiteral a,Int k){assert(k>=0&&!a.node->left&&a.node->degree==1&&a.node->coefficient==T(1));return SparseLiteral(std::make_shared<Node>(Node{k,T(1),{},{},false}));}
    SparseFPS<T> terms()const{SparseFPS<T> out;auto visit=[&](auto&& self,const std::shared_ptr<const Node>& p,bool neg)->void{if(!p->left){out.push_back({p->degree,neg?-p->coefficient:p->coefficient});return;}self(self,p->left,neg);self(self,p->right,neg!=p->subtract);};visit(visit,node,false);return out;}
};
template<Modint T> auto sfps(const SparseLiteral<T>& expression){return initSparseFPS<T>(expression.terms());}
template<Modint T,class Expression> requires std::invocable<Expression,SparseLiteral<T>> auto sfps(Expression expression){return sfps<T>(SparseLiteral<T>(expression(SparseLiteral<T>::variable())));}
template<Modint T,class Expression> auto SFPS(Expression expression){return sfps<T>(std::move(expression));}
template<class T,std::integral I> SparseFPS<T> operator+(SparseFPS<T> f,I c){T value=T(c);if(value.val()==0)return f;if(f.empty()||f[0].degree>0)f.insert(f.begin(),{0,value});else{f[0].coefficient+=value;if(f[0].coefficient.val()==0)f.erase(f.begin());}return f;}
template<class T,std::integral I> auto operator+(I c,const SparseFPS<T>& f){return f+c;}
template<class T,std::integral I> auto& operator+=(SparseFPS<T>& f,I c){return f=f+c;}
template<class T> bool isZero(const SparseFPS<T>& f){return f.empty();}
template<class T> Int degree(const SparseFPS<T>& f){return f.empty()?-1:f.back().degree;}
template<class T> T coefficient(const SparseFPS<T>& f,Int degree){if(degree<0)return T(0);auto it=std::lower_bound(f.begin(),f.end(),degree,[](const auto& term,Int d){return term.degree<d;});return it!=f.end()&&it->degree==degree?it->coefficient:T(0);}
template<class T> T constantTerm(const SparseFPS<T>& f){return coefficient(f,0);}
template<class T> std::vector<T> toDense(const SparseFPS<T>& f,Int n){if(n<=0)return {};std::vector<T> out(n);for(auto term:f){if(term.degree>=n)break;out[term.degree]=term.coefficient;}return out;}
template<class T> std::vector<T> mulPrefix(const std::vector<T>& f,const SparseFPS<T>& g,Int n){if(n<=0)return {};std::vector<T> out(n);for(auto term:g){if(term.degree>=n)break;for(Int i=0;i<std::min(Int(f.size()),n-term.degree);++i)out[i+term.degree]+=f[i]*term.coefficient;}return out;}
template<class T> std::vector<T> operator*(const std::vector<T>& f,const SparseFPS<T>& g){if(f.empty()||g.empty())return {};return mulPrefix(f,g,Int(f.size())+degree(g));}
template<class T> auto operator*(const SparseFPS<T>& f,const std::vector<T>& g){return g*f;}
template<class T> auto& operator*=(std::vector<T>& f,const SparseFPS<T>& g){return f=mulPrefix(f,g,f.size());}
namespace detail {
template<class T> SparseFPS<T> sparseTruncated(const SparseFPS<T>& f,Int n){SparseFPS<T> out;for(auto term:f){if(term.degree>=n)break;out.push_back(term);}return out;}
template<class T> void sparseDivideInPlace(std::vector<T>& f,const SparseFPS<T>& g,Int inputLen){T constant=constantTerm(g);assert(constant.val()!=0);T cinv=constant.inv();SparseFPS<T> terms;for(auto term:g){if(term.degree>=Int(f.size()))break;if(term.degree>0)terms.push_back({term.degree,term.coefficient*cinv});}if(constant.val()!=1)for(Int i=0;i<std::min(inputLen,Int(f.size()));++i)f[i]*=cinv;if(terms.empty())return;for(Int i=terms[0].degree;i<Int(f.size());++i){T value=f[i];for(auto term:terms){if(term.degree>i)break;value-=term.coefficient*f[i-term.degree];}f[i]=value;}}
// 非零定数項の冪。重みを次数ごとに差分更新する元実装。O(NK)。
template<Modint T> std::vector<T> sparsePowUnit(const SparseFPS<T>& f,T exponent,T constantRoot,Int n){if(n<=0)return {};T constant=constantTerm(f);assert(constant.val()!=0&&n<=T::umod());std::vector<T> out(n);out[0]=constantRoot;T cinv=constant.inv();struct Term {Int degree;T coefficient,weight;};std::vector<Term> terms;for(auto term:f){if(term.degree>=n)break;if(term.degree==0)continue;T c=term.coefficient*cinv;terms.push_back({term.degree,c,exponent*T(term.degree)*c});}if(terms.empty())return out;auto inverses=fpsIndexInverses<T>(n);for(Int degree=terms[0].degree;degree<n;++degree){T value=0;for(auto& term:terms){if(term.degree>degree)break;value+=term.weight*out[degree-term.degree];term.weight-=term.coefficient;}out[degree]=value*inverses[degree];}return out;}
}
// 密な分子を疎な分母で除算する。非零項Kに対してO(NK)。
template<class T> std::vector<T> divPrefix(const std::vector<T>& f,const SparseFPS<T>& g,Int n){if(n<=0)return {};auto out=prefix(f,n);detail::sparseDivideInPlace(out,g,f.size());return out;}
template<class T> auto operator/(const std::vector<T>& f,const SparseFPS<T>& g){return divPrefix(f,g,f.size());}
template<class T> auto& operator/=(std::vector<T>& f,const SparseFPS<T>& g){return f=divPrefix(f,g,f.size());}
template<class T> auto inv(const SparseFPS<T>& f,Int n){return divPrefix(std::vector<T>{T(1)},f,n);}
template<class T> std::vector<T> exp(const SparseFPS<T>& f,Int n){if(n<=0)return {};assert(n<=T::umod()&&constantTerm(f).val()==0);std::vector<T> out(n);out[0]=1;SparseFPS<T> terms;for(auto term:f){if(term.degree>=n)break;if(term.degree>0)terms.push_back({term.degree,term.coefficient*T(term.degree)});}if(terms.empty())return out;auto inverses=detail::fpsIndexInverses<T>(n);for(Int degree=terms[0].degree;degree<n;++degree){T value=0;for(auto term:terms){if(term.degree>degree)break;value+=term.coefficient*out[degree-term.degree];}out[degree]=value*inverses[degree];}return out;}
template<class T> std::vector<T> log(const SparseFPS<T>& f,Int n){if(n<=0)return {};assert(n<=T::umod()&&constantTerm(f).val()==1);std::vector<T> out(n);if(f.size()<=1||f[1].degree>=n)return out;for(auto term:f){if(term.degree>=n)break;out[term.degree]=term.coefficient*T(term.degree);}detail::sparseDivideInPlace(out,f,n);auto inverses=detail::fpsIndexInverses<T>(n);for(Int d=1;d<n;++d)out[d]*=inverses[d];return out;}
template<class T> std::vector<T> pow(const SparseFPS<T>& f,Int k,Int n){assert(k>=0);if(n<=0)return {};assert(n<=T::umod());if(k==0){std::vector<T> out(n);out[0]=1;return out;}if(k==1)return toDense(f,n);if(f.empty()||f[0].degree>(n-1)/k)return std::vector<T>(n);Int order=f[0].degree,shift=order*k;T leading=f[0].coefficient;SparseFPS<T> unit;for(auto term:f){if(term.degree-order>=n-shift)break;unit.push_back({term.degree-order,term.coefficient});}auto body=detail::sparsePowUnit(unit,T(k%T::umod()),leading.pow(k),n-shift);if(shift==0)return body;std::vector<T> out(n);std::copy(body.begin(),body.end(),out.begin()+shift);return out;}
template<class T> std::optional<std::vector<T>> sqrt(const SparseFPS<T>& f,Int n){if(n<=0)return std::vector<T>{};assert(n<=T::umod());if(f.empty()||f[0].degree>=n)return std::vector<T>(n);Int order=f[0].degree;if(order&1)return {};Int shift=order/2;auto leadingRoot=sqrt(std::vector<T>{f[0].coefficient},1);if(!leadingRoot)return {};SparseFPS<T> unit;for(auto term:f){if(term.degree-order>=n-shift)break;unit.push_back({term.degree-order,term.coefficient});}auto body=detail::sparsePowUnit(unit,T(1)/T(2),(*leadingRoot)[0],n-shift);if(shift==0)return body;std::vector<T> out(n);std::copy(body.begin(),body.end(),out.begin()+shift);return out;}
namespace detail {
template<class T> struct LinearPolynomial {T constant{},linear{};};
template<Modint T> auto multiplyPolynomialMatrices(const std::vector<std::vector<T>>& a,const std::vector<std::vector<T>>& b,Int d){std::vector<std::vector<T>> out(d*d);for(Int row=0;row<d;++row)for(Int middle=0;middle<d;++middle){const auto& left=a[row*d+middle];if(left.empty())continue;for(Int col=0;col<d;++col){const auto& right=b[middle*d+col];if(!right.empty())out[row*d+col]+=convolution(left,right);}}return out;}
template<Modint T> std::vector<T> shiftedLinearPolynomial(LinearPolynomial<T> p,Int shift){T constant=p.constant+p.linear*T(shift);if(p.linear.val()!=0)return {constant,p.linear};if(constant.val()!=0)return {constant};return {};}
template<Modint T> std::vector<std::vector<T>> blockMatrixProduct(const std::vector<LinearPolynomial<T>>& transition,Int dimension,Int left,Int right){if(right-left==1){std::vector<std::vector<T>> out(dimension*dimension);for(Int i=0;i<Int(out.size());++i)out[i]=shiftedLinearPolynomial(transition[i],left);return out;}Int middle=(left+right)>>1;auto lower=blockMatrixProduct(transition,dimension,left,middle),upper=blockMatrixProduct(transition,dimension,middle,right);return multiplyPolynomialMatrices(upper,lower,dimension);}
template<Modint T> std::vector<T> blockScalarProduct(LinearPolynomial<T> p,Int left,Int right){if(right-left==1)return shiftedLinearPolynomial(p,left);Int middle=(left+right)>>1;return convolution(blockScalarProduct(p,left,middle),blockScalarProduct(p,middle,right));}
template<Modint T> T evaluateLinearPolynomial(LinearPolynomial<T> p,Int x){return p.constant+p.linear*T(x);}
// 一次式の有理行列漸化式。平方分割したブロック行列を多点評価する。
template<Modint T> T nthTermPolynomialRecurrence(const std::vector<T>& initial,const std::vector<LinearPolynomial<T>>& transition,LinearPolynomial<T> denominator,Int target){Int d=initial.size();assert(d>0&&Int(transition.size())==d*d);if(target<d)return initial[target];Int count=target-d+1,blockSize=isqrt(count);if(blockSize*blockSize<count)++blockSize;Int blockCount=count/blockSize;auto state=initial;if(blockCount>0){auto matrix=blockMatrixProduct(transition,d,0,blockSize);auto denom=blockScalarProduct(denominator,0,blockSize);std::vector<T> points(blockCount);for(Int i=0;i<blockCount;++i)points[i]=T(i*blockSize);auto tree=initPolynomialProductTree(points);std::vector<std::vector<T>> values(d*d);for(Int i=0;i<d*d;++i)values[i]=matrix[i].empty()?std::vector<T>(blockCount):tree.evaluate(matrix[i]);auto dv=tree.evaluate(denom);for(Int block=0;block<blockCount;++block){assert(dv[block].val()!=0);std::vector<T> next(d);for(Int row=0;row<d;++row){for(Int col=0;col<d;++col)next[row]+=values[row*d+col][block]*state[col];next[row]/=dv[block];}state=std::move(next);}}for(Int step=blockCount*blockSize;step<count;++step){T denom=evaluateLinearPolynomial(denominator,step);assert(denom.val()!=0);std::vector<T> next(d);for(Int row=0;row<d;++row){for(Int col=0;col<d;++col)next[row]+=evaluateLinearPolynomial(transition[row*d+col],step)*state[col];next[row]/=denom;}state=std::move(next);}return state.back();}
template<class T> struct PolynomialRecurrence {std::vector<LinearPolynomial<T>> matrix;LinearPolynomial<T> denominator;};
template<Modint T> PolynomialRecurrence<T> expTransition(const SparseFPS<T>& f,Int d){PolynomialRecurrence<T> out{std::vector<LinearPolynomial<T>>(d*d),{T(d),T(1)}};for(Int row=0;row<d-1;++row)out.matrix[row*d+row+1]=out.denominator;for(auto term:f)if(term.degree>0&&term.degree<=d)out.matrix[(d-1)*d+d-term.degree].constant=T(term.degree)*term.coefficient;return out;}
template<Modint T> PolynomialRecurrence<T> powTransition(const SparseFPS<T>& f,T exponent,Int d){T constant=constantTerm(f);PolynomialRecurrence<T> out{std::vector<LinearPolynomial<T>>(d*d),{T(d)*constant,constant}};for(Int row=0;row<d-1;++row)out.matrix[row*d+row+1]=out.denominator;for(auto term:f)if(term.degree>0&&term.degree<=d)out.matrix[(d-1)*d+d-term.degree]={( (exponent+T(1))*T(term.degree)-T(d))*term.coefficient,-term.coefficient};return out;}
}
// 単項取得。逆元と対数はO(M(d)log k)、指数と冪はO(d³M(sqrt(k))log k)。
template<class T> T invCoefficient(const SparseFPS<T>& f,Int degree){assert(degree>=0&&constantTerm(f).val()!=0);auto relevant=detail::sparseTruncated(f,degree+1);return bostanMori(std::vector<T>{T(1)},toDense(relevant,cplib::degree(relevant)+1),degree);}
template<class T> T logCoefficient(const SparseFPS<T>& f,Int degree){assert(degree>=0&&constantTerm(f).val()==1);if(degree==0)return T(0);assert(degree<T::umod());auto relevant=detail::sparseTruncated(f,degree+1);auto denominator=toDense(relevant,cplib::degree(relevant)+1);std::vector<T> numerator(denominator.size()-1);for(auto term:relevant)if(term.degree>0)numerator[term.degree-1]=T(term.degree)*term.coefficient;return bostanMori(numerator,denominator,degree-1)/T(degree);}
template<class T> T expCoefficient(const SparseFPS<T>& f,Int degree){assert(degree>=0&&degree<T::umod()&&constantTerm(f).val()==0);if(degree==0)return T(1);auto relevant=detail::sparseTruncated(f,degree+1);Int d=cplib::degree(relevant);if(d<=0)return T(0);auto initial=exp(relevant,std::min(d,degree+1));if(degree<d)return initial[degree];auto rec=detail::expTransition(relevant,d);return detail::nthTermPolynomialRecurrence(initial,rec.matrix,rec.denominator,degree);}
template<class T> T powCoefficient(const SparseFPS<T>& f,Int k,Int degree){assert(k>=0&&degree>=0);if(k==0)return T(degree==0);if(f.empty()||f[0].degree>degree/k)return T(0);Int order=f[0].degree,shiftedDegree=degree-order*k;T leading=f[0].coefficient,scale=leading.pow(k);if(shiftedDegree==0)return scale;assert(shiftedDegree<T::umod());SparseFPS<T> terms;for(auto term:f){Int d=term.degree-order;if(d>shiftedDegree)break;terms.push_back({d,term.coefficient/leading});}auto unit=initSparseFPS<T>(std::move(terms));Int d=cplib::degree(unit);if(d<=0)return T(0);T exponent=T(k%T::umod());auto initial=detail::sparsePowUnit(unit,exponent,T(1),std::min(d,shiftedDegree+1));if(shiftedDegree<d)return initial[shiftedDegree]*scale;auto rec=detail::powTransition(unit,exponent,d);return detail::nthTermPolynomialRecurrence(initial,rec.matrix,rec.denominator,shiftedDegree)*scale;}
}
