#pragma once
#include <cplib/common.hpp>
#include <sstream>
namespace cplib {
template<class T> class MaxHeapQueue {
 std::vector<T> data;
 void up(Int i){T x=std::move(data[i]);while(i&&data[(i-1)/2]<x){data[i]=std::move(data[(i-1)/2]);i=(i-1)/2;}data[i]=std::move(x);}
 void down(Int i){T x=std::move(data[i]);for(Int c=i*2+1;c<len();c=i*2+1){if(c+1<len()&&!(data[c+1]<data[c]))++c;if(!(x<data[c]))break;data[i]=std::move(data[c]);i=c;}data[i]=std::move(x);}
public:
 MaxHeapQueue()=default;
 explicit MaxHeapQueue(std::span<const T> v):data(v.begin(),v.end()){for(Int i=len()/2;i-->0;)down(i);}
 Int len()const{return data.size();}const T& operator[](Int i)const{return data.at(i);}
 auto begin()const{return data.begin();}auto end()const{return data.end();}
 void push(T x){data.push_back(std::move(x));up(len()-1);}
 T pop(){assert(len());T x=std::move(data[0]);T last=std::move(data.back());data.pop_back();if(len()){data[0]=std::move(last);down(0);}return x;}
 Int find(const T& x)const{auto it=std::find(data.begin(),data.end(),x);return it==data.end()?-1:it-data.begin();}
 bool contains(const T& x)const{return find(x)>=0;}
 void del(Int i){assert(0<=i&&i<len());std::swap(data[i],data.back());data.pop_back();if(i<len()){if(i&&data[(i-1)/2]<data[i])up(i);else down(i);}}
 T replace(T x){assert(len());std::swap(x,data[0]);down(0);return x;}
 T pushpop(T x){if(len()&&x<data[0]){std::swap(x,data[0]);down(0);}return x;}
 void clear(){data.clear();}
 std::string str()const{std::ostringstream s;s<<'[';bool first=true;for(auto& x:data){if(!first)s<<", ";first=false;s<<x;}return s.str()+']';}
};
template<class T> auto initMaxHeapQueue(){return MaxHeapQueue<T>();}
template<class R> auto toMaxHeapQueue(const R& r){using T=typename R::value_type;return MaxHeapQueue<T>(r);}
}
