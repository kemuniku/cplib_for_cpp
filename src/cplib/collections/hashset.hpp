#pragma once
#include <cplib/common.hpp>
#include <functional>
namespace cplib {
namespace detail {enum class OpenHashState {empty,active,inactive};inline std::size_t open_hash_capacity(std::size_t n){return n?std::size_t(1)<<(std::bit_width(n)+1):4;}}
// 線形探索と削除済み印による集合。再挿入時の fill 加算も元実装と同じ。
template<class T,class Hash=std::hash<T>> class HashSet {
    std::size_t length_=0,fill_=0,mask_=3;
    std::size_t find(const T& x)const{auto p=Hash{}(x)&mask_;while(values[p].state!=detail::OpenHashState::empty&&!(values[p].value==x))p=(p+1)&mask_;return p;}
    void add_item(const T& x){auto p=find(x);if(values[p].state==detail::OpenHashState::active)return;++length_;++fill_;values[p]={detail::OpenHashState::active,x};}
    void resize(){auto old=std::move(values);values.resize(detail::open_hash_capacity(length_));mask_=values.size()-1;length_=fill_=0;for(const auto& n:old)if(n.state==detail::OpenHashState::active)add_item(n.value);}
public:
    struct Node {detail::OpenHashState state=detail::OpenHashState::empty;T value{};};
    std::vector<Node> values=std::vector<Node>(4);
    Int len()const{return length_;}
    void incl(const T& x){add_item(x);if(detail::open_hash_capacity(fill_)>values.size())resize();}
    bool contains(const T& x)const{return values[find(x)].state==detail::OpenHashState::active;}
    void excl(const T& x){auto p=find(x);if(values[p].state!=detail::OpenHashState::active)return;--length_;values[p].state=detail::OpenHashState::inactive;}
    struct Iterator {const HashSet* owner;std::size_t pos;void skip(){while(pos<owner->values.size()&&owner->values[pos].state!=detail::OpenHashState::active)++pos;}const T& operator*()const{return owner->values[pos].value;}Iterator& operator++(){++pos;skip();return *this;}bool operator==(const Iterator&)const=default;};
    Iterator begin()const{Iterator it{this,0};it.skip();return it;}Iterator end()const{return {this,values.size()};}const HashSet& items()const{return *this;}
};
template<class T> HashSet<T> initHashSet(){return {};}
template<class T> HashSet<T> toHashSet(const std::vector<T>& v){HashSet<T> out;for(const T& x:v)out.incl(x);return out;}
}
