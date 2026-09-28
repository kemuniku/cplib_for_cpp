#pragma once
#include <cplib/geometry/intersect.hpp>
namespace cplib {
template<class T> T norm(Point<T> a,Point<T> b){return norm(a-b);}
template<class T> T norm(Point<T> p,const Line<T>& l){T area=cross(vector(l),p-l.s);return area*area/norm(vector(l));}
template<class T> T norm(Point<T> p,const Segment<T>& s){if(geometry_lt(dot(vector(s),p-s.s),0))return norm(p-s.s);if(geometry_lt(dot(-vector(s),p-s.t),0))return norm(p-s.t);return norm(p,toLine(s));}
template<class T> T norm(const Segment<T>& a,const Segment<T>& b){if(intersect(a,b))return T(0);return std::min({norm(a.s,b),norm(a.t,b),norm(b.s,a),norm(b.t,a)});}
template<std::floating_point T> double distance(Point<T> a,Point<T> b){return std::sqrt(norm(a,b));}
template<std::floating_point T> double distance(Point<T> p,const Line<T>& l){return std::sqrt(norm(p,l));}
template<std::floating_point T> double distance(Point<T> p,const Segment<T>& s){return std::sqrt(norm(p,s));}
template<std::floating_point T> double distance(const Segment<T>& a,const Segment<T>& b){return std::sqrt(norm(a,b));}
template<class T> T manhattan(Point<T> p){using std::abs;return abs(p.x)+abs(p.y);}
template<class T> T manhattan(Point<T> a,Point<T> b){return manhattan(a-b);}
}
