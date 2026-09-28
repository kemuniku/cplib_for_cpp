#pragma once
#include <cplib/common.hpp>
#include <cmath>
#include <functional>
#include <sstream>
#include <type_traits>
namespace cplib {
inline double GEOMETRY_EPS=1e-10;
// 元の相対誤差判定は左辺xをそのまま用いる（abs(x)へは置換しない）。
template<class T,class S> bool geometry_eq(const T& x,const S& y){if constexpr(std::is_floating_point_v<T>||std::is_floating_point_v<S>){auto diff=std::abs(x-y);return diff<GEOMETRY_EPS*x||diff<GEOMETRY_EPS;}else return x==y;}
template<class T,class S> bool geometry_neq(const T& x,const S& y){return !geometry_eq(x,y);}
template<class T,class S> bool geometry_ge(const T& x,const S& y){return geometry_eq(x,y)||x>y;}
template<class T,class S> bool geometry_le(const T& x,const S& y){return geometry_eq(x,y)||x<y;}
template<class T,class S> bool geometry_gt(const T& x,const S& y){return !geometry_le(x,y);}
template<class T,class S> bool geometry_lt(const T& x,const S& y){return !geometry_ge(x,y);}
template<class T> struct Point {
    T x{},y{};
    Point operator+()const{return *this;}Point operator-()const{return {-x,-y};}
    Point& operator+=(Point q){x+=q.x;y+=q.y;return *this;}Point& operator-=(Point q){x-=q.x;y-=q.y;return *this;}
    template<class S> Point& operator*=(S c){if constexpr(std::is_floating_point_v<T>){x*=double(c);y*=double(c);}else{x*=c;y*=c;}return *this;}
    template<class S> Point& operator/=(S c){if constexpr(std::is_floating_point_v<T>){x/=double(c);y/=double(c);}else{x/=c;y/=c;}return *this;}
    friend Point operator+(Point p,Point q){return p+=q;}friend Point operator-(Point p,Point q){return p-=q;}
    template<class S> requires(!std::is_same_v<S,Point>) friend Point operator*(Point p,S c){return p*=c;}
    template<class S> requires(!std::is_same_v<S,Point>) friend Point operator*(S c,Point p){return p*=c;}
    template<class S> friend Point operator/(Point p,S c){return p/=c;}
    friend T operator*(Point p,Point q){return p.x*q.x+p.y*q.y;}
    friend bool operator<(Point p,Point q){return geometry_eq(p.x,q.x)?geometry_lt(p.y,q.y):geometry_lt(p.x,q.x);}
    friend bool operator>(Point p,Point q){return geometry_eq(p.x,q.x)?geometry_gt(p.y,q.y):geometry_gt(p.x,q.x);}
    friend bool operator==(Point p,Point q){return geometry_eq(p.x,q.x)&&geometry_eq(p.y,q.y);}
    friend bool operator<=(Point p,Point q){return !(p>q);}friend bool operator>=(Point p,Point q){return !(p<q);}
    friend std::ostream& operator<<(std::ostream& out,Point p){return out<<'('<<p.x<<", "<<p.y<<')';}
};
template<class T> auto initPoint(T x,T y){return Point<T>{x,y};}
template<class T> auto initPoint(std::pair<T,T> p){return Point<T>{p.first,p.second};}
template<class T> T dot(Point<T> p,Point<T> q){return p*q;}
template<class T> T cross(Point<T> p,Point<T> q){return p.x*q.y-p.y*q.x;}
template<class T> T norm(Point<T> p){return dot(p,p);}
template<class T> Int cmp(Point<T> p,Point<T> q){return p<q?-1:p==q?0:1;}
template<class T> std::string to_string(Point<T> p){std::ostringstream out;out<<p;return out.str();}
template<class T> struct Line {Point<T> s,t;friend std::ostream& operator<<(std::ostream& out,const Line& l){return out<<'('<<l.s<<", "<<l.t<<')';}};
template<class T> auto initLine(Point<T> s,Point<T> t){assert(s!=t);return Line<T>{s,t};}
template<class T> requires(!std::is_integral_v<T>) Line<T> initLine(T a,T b,T c){assert(geometry_neq(a,T(0))||geometry_neq(b,T(0)));if(geometry_eq(b,T(0)))return {{-c/a,T(0)},{-c/a,T(1)}};return {{T(0),-c/b},{T(1),(-a-c)/b}};}
template<class T> Point<T> vector(const Line<T>& l){return l.t-l.s;}
template<class T> std::string to_string(const Line<T>& l){std::ostringstream out;out<<l;return out.str();}
template<class T> struct Segment {Point<T> s,t;operator Line<T>()const{return initLine(s,t);}};
template<class T> auto initSegment(Point<T> s,Point<T> t){assert(s!=t);return Segment<T>{s,t};}
template<class T> auto toLine(const Segment<T>& s){return initLine(s.s,s.t);}
template<class T> auto vector(const Segment<T>& s){return s.t-s.s;}
}
namespace std {template<class T> struct hash<cplib::Point<T>> {size_t operator()(cplib::Point<T> p)const{size_t a=hash<T>{}(p.x),b=hash<T>{}(p.y);return a^(b+0x9e3779b97f4a7c15ULL+(a<<6)+(a>>2));}};}
