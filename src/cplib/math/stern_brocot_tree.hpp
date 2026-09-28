#pragma once
#include <cplib/math/fractions.hpp>
#include <cplib/graph/graph.hpp>
#include <optional>
namespace cplib {
template<class T> struct SBTNode {
    T p,q,r,s,depth;
    T den()const{return q+s;}T num()const{return p+r;}
    operator Fraction<T>()const{return Fraction<T>(num(),den(),false);}
    bool operator==(const SBTNode&)const=default;
};
template<class T> T num(SBTNode<T> x){return x.num();}
template<class T> T den(SBTNode<T> x){return x.den();}
template<class T=Int> auto sbt_root(){return SBTNode<T>{T(0),T(1),T(1),T(0),T(0)};}
template<class T=Int> auto sbt_inf(){return SBTNode<T>{T(1),T(0),T(0),T(0),T(-1)};}
template<class T=Int> auto sbt_zero(){return SBTNode<T>{T(0),T(1),T(0),T(0),T(-1)};}
// 連分数展開と圧縮経路。Euclidの除算回数に比例する時間・領域。
template<class T> std::vector<T> continued_fraction_expansion(T a,T b){assert(a>=1&&b>=1);std::vector<T> result;while(true){result.push_back(a/b);if(a%b==0)break;a-=result.back()*b;std::swap(a,b);}return result;}
template<class T> std::vector<std::pair<char,T>> encode_path(T a,T b){auto c=continued_fraction_expansion(a,b);--c.back();std::vector<std::pair<char,T>> result;for(Int i=c[0]==0?1:0;i<Int(c.size());++i)result.emplace_back(i%2?'L':'R',c[i]);return result;}
template<class T> auto encode_path(SBTNode<T> x){return encode_path(x.num(),x.den());}
template<class T> SBTNode<T> move_left(SBTNode<T> x,T d){return {x.p,x.q,d*x.p+x.r,d*x.q+x.s,x.depth+d};}
template<class T> SBTNode<T> move_right(SBTNode<T> x,T d){return {x.p+d*x.r,x.q+d*x.s,x.r,x.s,x.depth+d};}
template<class T> auto to_fraction_tuple(SBTNode<T> x){return std::pair<T,T>{x.num(),x.den()};}
template<class T> auto to_SBTNode(T a,T b){auto now=sbt_root<T>();for(auto [c,d]:encode_path(a,b))now=c=='L'?move_left(now,d):move_right(now,d);return now;}
namespace detail {
template<class T> auto sbt_endpoint(T a,T b){return b==0?sbt_inf<T>():a==0?sbt_zero<T>():to_SBTNode(a,b);}
template<class T> T sbt_max_move_count(T base,T step,T n){if(base>n)return T(-1);if(step==0)return n;return (n-base)/step;}
template<class T> bool sbt_inner_bounded(SBTNode<T> x,T n){return x.num()<=n&&x.den()<=n;}
template<class T> T sbt_inner_left(SBTNode<T> x,T n){return std::min(sbt_max_move_count(x.num(),x.p,n),sbt_max_move_count(x.den(),x.q,n));}
template<class T> T sbt_inner_right(SBTNode<T> x,T n){return std::min(sbt_max_move_count(x.num(),x.r,n),sbt_max_move_count(x.den(),x.s,n));}
template<class T> T sbt_endpoint_left(SBTNode<T> x,T n){return std::min(sbt_max_move_count(x.r,x.p,n),sbt_max_move_count(x.s,x.q,n));}
template<class T> T sbt_endpoint_right(SBTNode<T> x,T n){return std::min(sbt_max_move_count(x.p,x.r,n),sbt_max_move_count(x.q,x.s,n));}
}
// 分子・分母がm以下の、直後・直前の有理数。存在しない場合は∞または0。
template<class T> SBTNode<T> min_greater_with_den_at_most(SBTNode<T> x,T m){assert(m>=1);if(detail::sbt_inner_bounded(x,m)){auto now=move_right(x,T(1));if(!detail::sbt_inner_bounded(now,m))return detail::sbt_endpoint(x.r,x.s);return move_left(now,detail::sbt_inner_left(now,m));}auto now=sbt_root<T>();for(auto [c,d]:encode_path(x)){if(c=='L'){T lim=std::min(d,detail::sbt_inner_left(now,m));now=move_left(now,lim);if(lim<d)return now;}else{T lim=std::min(d,detail::sbt_inner_right(now,m));now=move_right(now,lim);if(lim<d)return detail::sbt_endpoint(now.r,now.s);}}auto next=move_right(now,T(1));if(!detail::sbt_inner_bounded(next,m))return detail::sbt_endpoint(now.r,now.s);return move_left(next,detail::sbt_inner_left(next,m));}
template<class T> SBTNode<T> max_less_with_den_at_most(SBTNode<T> x,T m){assert(m>=1);if(detail::sbt_inner_bounded(x,m)){auto now=move_left(x,T(1));if(!detail::sbt_inner_bounded(now,m))return detail::sbt_endpoint(x.p,x.q);return move_right(now,detail::sbt_inner_right(now,m));}auto now=sbt_root<T>();for(auto [c,d]:encode_path(x)){if(c=='L'){T lim=std::min(d,detail::sbt_inner_left(now,m));now=move_left(now,lim);if(lim<d)return detail::sbt_endpoint(now.p,now.q);}else{T lim=std::min(d,detail::sbt_inner_right(now,m));now=move_right(now,lim);if(lim<d)return now;}}auto next=move_left(now,T(1));if(!detail::sbt_inner_bounded(next,m))return detail::sbt_endpoint(now.p,now.q);return move_right(next,detail::sbt_inner_right(next,m));}
template<class T> auto decode_path(const std::vector<std::pair<char,T>>& path){auto now=sbt_root<T>();for(auto [c,d]:path)now=c=='L'?move_left(now,d):move_right(now,d);return now;}
template<class T> auto LCA(T a,T b,T c,T d){auto x=continued_fraction_expansion(a,b),y=continued_fraction_expansion(c,d);--x.back();--y.back();auto now=sbt_root<T>();for(Int i=0;i<Int(std::min(x.size(),y.size()));++i){now=i%2?move_left(now,std::min(x[i],y[i])):move_right(now,std::min(x[i],y[i]));if(x[i]!=y[i])return now;}return now;}
template<class T> auto LCA(SBTNode<T> a,SBTNode<T> b){return LCA(a.num(),a.den(),b.num(),b.den());}
template<class T> std::optional<SBTNode<T>> ancestor(T a,T b,T k){auto c=continued_fraction_expansion(a,b);--c.back();auto now=sbt_root<T>();T count=0;for(Int i=0;i<Int(c.size());++i){T d=std::min(T(k-count),c[i]);now=i%2?move_left(now,d):move_right(now,d);count+=c[i];if(count>=k)return now;}return std::nullopt;}
template<class T> auto ancestor(SBTNode<T> x,T k){return ancestor(x.num(),x.den(),k);}
template<class T> auto get_range(SBTNode<T> x){return std::tuple{x.p,x.q,x.r,x.s};}
template<class T> auto get_range(T a,T b){return get_range(to_SBTNode(a,b));}
template<class T> auto get_range_fraction(SBTNode<T> x){return std::pair{Fraction<T>(x.p,x.q,false),Fraction<T>(x.r,x.s,false)};}
// 単調述語の境界を、分子・分母n以下の二つの端点で挟む。圧縮経路上で二分探索。
template<class Pred,class T> auto get_bounds(Pred is_ok,T n){assert(n>=1);auto now=sbt_root<T>();bool zero=is_ok(sbt_zero<T>()),infinity=is_ok(sbt_inf<T>());assert(zero!=infinity);bool current=is_ok(now),left=zero!=current;while(detail::sbt_inner_bounded(now,n)){T limit=left?detail::sbt_endpoint_left(now,n):detail::sbt_endpoint_right(now,n);if(limit<=0)break;auto move=[&](T d){return left?move_left(now,d):move_right(now,d);};T l=0,r=1;while(is_ok(move(r))==current){if(r==limit)return move(r);l=r;r+=std::min(r,T(limit-r));}while(r-l>1){T mid=l+(r-l)/2;if(is_ok(move(mid))==current)l=mid;else r=mid;}now=move(r);current=!current;left=!left;}return now;}
}
namespace std {
template<class T> struct hash<cplib::SBTNode<T>>{size_t operator()(const cplib::SBTNode<T>& x)const{size_t h=0;for(const T& v:{x.p,x.q,x.r,x.s,x.depth})h^=hash<T>{}(v)+0x9e3779b97f4a7c15ULL+(h<<6)+(h>>2);return h;}};
}
namespace cplib {
namespace detail {
inline int sbt_cmp_ein(SBTNode<Int> a,SBTNode<Int> b){if(a==b)return 0;auto [la,ra]=get_range_fraction(a);auto [lb,rb]=get_range_fraction(b);Fraction<Int> x=a,y=b;if(lb<x&&x<rb)return 1;if(la<y&&y<ra)return -1;return x<y?-1:x>y?1:0;}
}
inline auto initAuxiliaryWeightedTree(std::span<const SBTNode<Int>> input){assert(!input.empty());std::vector<SBTNode<Int>> v(input.begin(),input.end());auto cmp=[](auto a,auto b){return detail::sbt_cmp_ein(a,b)<0;};std::sort(v.begin(),v.end(),cmp);Int n=v.size();for(Int i=0;i<n-1;++i)v.push_back(LCA(v[i],v[i+1]));std::sort(v.begin(),v.end(),cmp);v.erase(std::unique(v.begin(),v.end()),v.end());auto result=initWeightedUnDirectedTableGraph<Int>(v);std::vector<SBTNode<Int>> stack{v[0]};for(Int i=1;i<Int(v.size());++i){while(!stack.empty()&&LCA(stack.back(),v[i])!=stack.back())stack.pop_back();if(!stack.empty())result.add_edge(stack.back(),v[i],v[i].depth-stack.back().depth);stack.push_back(v[i]);}return result;}
}
