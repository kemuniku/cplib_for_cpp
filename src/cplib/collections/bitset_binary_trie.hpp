#pragma once
#include <cplib/utils/backwards_index.hpp>
#include <memory>
#include <stdexcept>
namespace cplib {
// 固定した長さのビット列の多重集合。添字0からfalse < trueの辞書順。
// キーは葉だけに保持し、削除したノードを再利用する。各操作O(1+N)、lenはO(1)。
template<class T> class BitSetBinaryTrie {
 struct Node{std::array<Int,2> children{-1,-1};Int count=0;std::shared_ptr<const T> key;};
 std::vector<Node> nodes_;std::vector<Int> free_;Int height_=0;
 void check(const T& x)const{if(nodes_.empty())throw std::invalid_argument("BitSetBinaryTrie is not initialized");if(x.len()!=height_)throw std::invalid_argument("BitSet length must match trie height");}
 Int new_node(){if(!free_.empty()){Int node=free_.back();free_.pop_back();return node;}nodes_.emplace_back();return nodes_.size()-1;}
 T kth(Int k,const T* mask)const{if(k<0||k>=len())throw std::out_of_range("BitSetBinaryTrie index out of bounds");Int node=0;for(Int i=0;i<height_;++i){int first=mask?bool((*mask)[i]):0;Int child=nodes_[node].children[first],count=child<0?0:nodes_[child].count;if(k<count)node=child;else{k-=count;node=nodes_[node].children[1-first];}}return *nodes_[node].key;}
 Int bound(const T& x,const T* mask,bool inclusive)const{check(x);Int node=0,result=0;for(Int i=0;i<height_;++i){int first=mask?bool((*mask)[i]):0,bit=bool(x[i]);if(bit){Int child=nodes_[node].children[first];if(child>=0)result+=nodes_[child].count;}node=nodes_[node].children[first^bit];if(node<0)return result;}return result+(inclusive?nodes_[node].count:0);}
public:
 BitSetBinaryTrie()=default;
 explicit BitSetBinaryTrie(Int h):height_(h){if(h<0)throw std::invalid_argument("BitSet length must be non-negative");nodes_.emplace_back();}
 explicit BitSetBinaryTrie(const T& sample):BitSetBinaryTrie(sample.len()){}
 Int len()const{return nodes_.empty()?0:nodes_[0].count;}
 Int count(const T& x)const{check(x);Int node=0;for(Int i=0;i<height_;++i){node=nodes_[node].children[bool(x[i])];if(node<0)return 0;}return nodes_[node].count;}
 bool contains(const T& x)const{return count(x)>0;}
 void incl(const T& x,Int multiplicity=1){check(x);if(multiplicity<0)throw std::invalid_argument("Multiplicity must be non-negative");if(multiplicity>std::numeric_limits<Int>::max()-len())throw std::overflow_error("BitSetBinaryTrie size overflow");if(!multiplicity)return;Int node=0;nodes_[node].count+=multiplicity;for(Int i=0;i<height_;++i){int bit=bool(x[i]);if(nodes_[node].children[bit]<0){Int child=new_node();nodes_[node].children[bit]=child;}node=nodes_[node].children[bit];nodes_[node].count+=multiplicity;}if(!nodes_[node].key)nodes_[node].key=std::make_shared<T>(x);}
 void excl(const T& x,Int multiplicity=1){
  check(x);if(multiplicity<0)throw std::invalid_argument("Multiplicity must be non-negative");if(!multiplicity)return;std::vector<Int> path{0};
  for(Int i=0;i<height_;++i){Int child=nodes_[path.back()].children[bool(x[i])];if(child<0)throw std::invalid_argument("Not enough copies of key");path.push_back(child);}if(nodes_[path.back()].count<multiplicity)throw std::invalid_argument("Not enough copies of key");
  for(Int node:path)nodes_[node].count-=multiplicity;
  if(!nodes_[path.back()].count)nodes_[path.back()].key.reset();
  for(Int i=height_;i>0;--i){Int node=path[i];if(nodes_[node].count)break;nodes_[path[i-1]].children[bool(x[i-1])]=-1;nodes_[node]=Node{};free_.push_back(node);}
 }
 T get_kth(Int k)const{return kth(k,nullptr);}T get_kth(Int k,const T& mask)const{check(mask);return kth(k,&mask);}
 T operator[](Int k)const{return get_kth(k);}T operator[](BackwardsIndex k)const{return get_kth(len()-k.value);}
 Int lowerBound(const T& x)const{return bound(x,nullptr,false);}Int upperBound(const T& x)const{return bound(x,nullptr,true);}
 Int lowerBound(const T& x,const T& mask)const{check(mask);return bound(x,&mask,false);}Int upperBound(const T& x,const T& mask)const{check(mask);return bound(x,&mask,true);}
};
template<class T> auto initBitSetBinaryTrie(Int h){return BitSetBinaryTrie<T>(h);}
template<class T> auto initBitSetBinaryTrie(const T& sample){return BitSetBinaryTrie<T>(sample);}
template<class T> Int len(const BitSetBinaryTrie<T>& trie){return trie.len();}
template<class T> Int count(const BitSetBinaryTrie<T>& trie,const T& x){return trie.count(x);}
template<class T> bool contains(const BitSetBinaryTrie<T>& trie,const T& x){return trie.contains(x);}
template<class T> void incl(BitSetBinaryTrie<T>& trie,const T& x,Int count=1){trie.incl(x,count);}
template<class T> void excl(BitSetBinaryTrie<T>& trie,const T& x,Int count=1){trie.excl(x,count);}
template<class T> T get_kth(const BitSetBinaryTrie<T>& trie,Int k){return trie.get_kth(k);}
template<class T> T get_kth(const BitSetBinaryTrie<T>& trie,Int k,const T& mask){return trie.get_kth(k,mask);}
template<class T> Int lowerBound(const BitSetBinaryTrie<T>& trie,const T& x){return trie.lowerBound(x);}
template<class T> Int upperBound(const BitSetBinaryTrie<T>& trie,const T& x){return trie.upperBound(x);}
template<class T> Int lowerBound(const BitSetBinaryTrie<T>& trie,const T& x,const T& mask){return trie.lowerBound(x,mask);}
template<class T> Int upperBound(const BitSetBinaryTrie<T>& trie,const T& x,const T& mask){return trie.upperBound(x,mask);}
}
