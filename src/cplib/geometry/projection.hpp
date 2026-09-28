#pragma once
#include <cplib/geometry/base.hpp>
namespace cplib {
template<class T> Point<T> projection(const Line<T>& l,Point<T> p){T t=dot(p-l.s,vector(l))/norm(vector(l));return l.s+vector(l)*t;}
template<class T> Point<T> reflection(const Line<T>& l,Point<T> p){return p+(projection(l,p)-p)*2;}
template<class T> auto projection(const Segment<T>& l,Point<T> p){return projection(toLine(l),p);}
template<class T> auto reflection(const Segment<T>& l,Point<T> p){return reflection(toLine(l),p);}
}
