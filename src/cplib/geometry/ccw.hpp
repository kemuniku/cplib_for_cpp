#pragma once
#include <cplib/geometry/base.hpp>
namespace cplib {
inline constexpr Int COUNTER_CLOCKWISE=2,CLOCKWISE=-2,ON_SEGMENT=0,ONLINE_BACK=1,ONLINE_FRONT=-1;
template<class T> Int ccw(const Line<T>& l,Point<T> p,bool strict=false){T c=cross(vector(l),p-l.s);if(geometry_lt(c,0))return CLOCKWISE;if(geometry_gt(c,0))return COUNTER_CLOCKWISE;T d=dot(vector(l),p-l.s);if(strict){if(geometry_le(d,0))return ONLINE_BACK;if(geometry_ge(d,norm(vector(l))))return ONLINE_FRONT;}else{if(geometry_lt(d,0))return ONLINE_BACK;if(geometry_gt(d,norm(vector(l))))return ONLINE_FRONT;}return ON_SEGMENT;}
template<class T> Int ccw(const Segment<T>& l,Point<T> p,bool strict=false){return ccw(toLine(l),p,strict);}
template<class T> Int ccw(Point<T> a,Point<T> b,Point<T> p,bool strict=false){return ccw(initLine(a,b),p,strict);}
template<class T> bool online(const Line<T>& l,Point<T> p){Int c=ccw(l,p);return -1<=c&&c<=1;}
template<class T> bool online(const Segment<T>& l,Point<T> p){return online(toLine(l),p);}
template<class T> bool online(Point<T> a,Point<T> b,Point<T> p){return online(initLine(a,b),p);}
}
