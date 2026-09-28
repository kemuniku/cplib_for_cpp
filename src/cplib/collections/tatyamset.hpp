#pragma once
#include <cplib/utils/backwards_index.hpp>
#include <cmath>
#include <optional>
#include <stdexcept>
namespace cplib {
// 平方分割による多重集合。バケット比 8、分割比 12 は原実装のまま。
template<class T> class SortedMultiSet {
    Int size_=0;
    std::pair<Int,Int> position(const T& x)const{for(Int i=0;i<Int(arr.size());++i)if(!(arr[i].back()<x))return {i,std::lower_bound(arr[i].begin(),arr[i].end(),x)-arr[i].begin()};return {Int(arr.size())-1,std::lower_bound(arr.back().begin(),arr.back().end(),x)-arr.back().begin()};}
    T innerpop(Int b,Int i){T out=arr[b][i];arr[b].erase(arr[b].begin()+i);--size_;if(arr[b].empty())arr.erase(arr.begin()+b);return out;}
    std::pair<Int,Int> locate(Int i)const{if(i<0){for(Int b=Int(arr.size());b-->0;){i+=arr[b].size();if(i>=0)return {b,i};}}else for(Int b=0;b<Int(arr.size());++b){if(i<Int(arr[b].size()))return {b,i};i-=arr[b].size();}throw std::out_of_range("SortedMultiSet index");}
public:
    std::vector<std::vector<T>> arr;
    explicit SortedMultiSet(std::vector<T> v={}){if(!std::is_sorted(v.begin(),v.end()))std::sort(v.begin(),v.end());size_=v.size();Int buckets=Int(std::ceil(std::sqrt(double(size_)/8)));for(Int i=0;i<buckets;++i)arr.emplace_back(v.begin()+size_*i/buckets,v.begin()+size_*(i+1)/buckets);}
    Int len()const{return size_;}
    bool contains(const T& x)const{if(!size_)return false;auto [b,i]=position(x);return i<Int(arr[b].size())&&arr[b][i]==x;}
    void incl(const T& x){if(!size_){arr={{x}};size_=1;return;}auto [b,i]=position(x);arr[b].insert(arr[b].begin()+i,x);++size_;if(arr[b].size()>arr.size()*12){auto mid=arr[b].size()/2;std::vector<T> right(arr[b].begin()+mid,arr[b].end());arr[b].resize(mid);arr.insert(arr.begin()+b+1,std::move(right));}}
    bool excl(const T& x){if(!size_)return false;auto [b,i]=position(x);if(i==Int(arr[b].size())||!(arr[b][i]==x))return false;innerpop(b,i);return true;}
    std::optional<T> lt(const T& x)const{for(Int b=Int(arr.size());b-->0;)if(arr[b].front()<x)return *std::prev(std::lower_bound(arr[b].begin(),arr[b].end(),x));return {};}
    std::optional<T> le(const T& x)const{for(Int b=Int(arr.size());b-->0;)if(!(x<arr[b].front()))return *std::prev(std::upper_bound(arr[b].begin(),arr[b].end(),x));return {};}
    std::optional<T> gt(const T& x)const{for(const auto& bucket:arr)if(x<bucket.back())return *std::upper_bound(bucket.begin(),bucket.end(),x);return {};}
    std::optional<T> ge(const T& x)const{for(const auto& bucket:arr)if(!(bucket.back()<x))return *std::lower_bound(bucket.begin(),bucket.end(),x);return {};}
    T operator[](Int i)const{auto [b,j]=locate(i);return arr[b][j];}T operator[](BackwardsIndex i)const{return (*this)[size_-i.value];}
    T pop(Int i=-1){auto [b,j]=locate(i);return innerpop(b,j);}
    Int index(const T& x)const{Int out=0;for(const auto& b:arr){if(!(b.back()<x))return out+std::lower_bound(b.begin(),b.end(),x)-b.begin();out+=b.size();}return out;}
    Int index_right(const T& x)const{Int out=0;for(const auto& b:arr){if(x<b.back())return out+std::upper_bound(b.begin(),b.end(),x)-b.begin();out+=b.size();}return out;}
    Int count(const T& x)const{return index_right(x)-index(x);}Int lowerBound(const T& x)const{return index(x);}Int upperBound(const T& x)const{return index_right(x);}
    struct Iterator {const SortedMultiSet* owner;std::size_t b=0,i=0;const T& operator*()const{return owner->arr[b][i];}Iterator& operator++(){if(++i==owner->arr[b].size()){++b;i=0;}return *this;}bool operator==(const Iterator&)const=default;};
    Iterator begin()const{return {this,0,0};}Iterator end()const{return {this,arr.size(),0};}const SortedMultiSet& items()const{return *this;}
};
template<class T> auto initSortedMultiset(const std::vector<T>& v={}){return SortedMultiSet<T>(v);}
}
