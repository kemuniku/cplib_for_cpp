#pragma once
#include <cplib/geometry/ccw.hpp>
#include <cplib/geometry/angle.hpp>
namespace cplib {
template<class T> bool intersect(const Segment<T>& a,const Segment<T>& b,bool strict=false){if(strict){if(ccw(a,b.s,true)==ON_SEGMENT)return online(a,b.t);if(ccw(a,b.t,true)==ON_SEGMENT)return online(a,b.s);if(ccw(b,b.s,true)==ON_SEGMENT)return online(b,b.t);if(ccw(b,b.t,true)==ON_SEGMENT)return online(b,b.s);return ccw(a,b.s)*ccw(a,b.t)<0&&ccw(b,a.s)*ccw(b,a.t)<0;}return ccw(a,b.s)*ccw(a,b.t)<=0&&ccw(b,a.s)*ccw(b,a.t)<=0;}
template<class T> bool intersect(const Line<T>& a,const Line<T>& b){return !is_parallel(a,b)||online(a,b.s);}
template<class T> requires(!std::is_integral_v<T>) Point<T> cross_point(const Line<T>& a,const Line<T>& b){assert(intersect(a,b));if(is_parallel(a,b))return a.s;T d1=cross(vector(a),vector(b)),d2=cross(vector(a),a.t-b.s);return b.s+vector(b)*(d2/d1);}
template<class T> requires(!std::is_integral_v<T>) auto cross_point(const Segment<T>& a,const Segment<T>& b){return cross_point(toLine(a),toLine(b));}
template<class T> bool intersect(const Line<T>& a,const Segment<T>& b){return intersect(a,toLine(b));}
template<class T> bool intersect(const Segment<T>& a,const Line<T>& b){return intersect(toLine(a),b);}
template<class T> requires(!std::is_integral_v<T>) auto cross_point(const Line<T>& a,const Segment<T>& b){return cross_point(a,toLine(b));}
template<class T> requires(!std::is_integral_v<T>) auto cross_point(const Segment<T>& a,const Line<T>& b){return cross_point(toLine(a),b);}

}
