#pragma once
#include <cplib/utils/backwards_index.hpp>
#include <cplib/utils/itertools.hpp>
#include <unordered_set>
#include <functional>
#include <cplib/collections/private/sort_order.hpp>
namespace cplib {

namespace detail {
struct RangeSortEmpty{};
template<class T,bool Products> class RangeSortEngine {
 struct Product{T ascending{},descending{};};
 struct Node:std::conditional_t<Products,Product,RangeSortEmpty>{Int left=0,right=0,count=0,key=0,bit=0;};
 std::vector<Node> nodes;std::vector<Int> freeNodes,roots,next,starts;std::vector<bool> reversed;std::vector<T> products;Int keyLimit,size=1;std::function<T(T,T)> merge;T e;
#ifndef NDEBUG
 std::unordered_set<Int> activeKeys;
#endif
 Int newNode(){if(!freeNodes.empty()){Int n=freeNodes.back();freeNodes.pop_back();return n;}nodes.emplace_back();return nodes.size()-1;}
 void releaseNode(Int n){nodes[n]={};freeNodes.push_back(n);}
 void pull(Int n){Int a=nodes[n].left,b=nodes[n].right;nodes[n].key=nodes[a].key;nodes[n].count=nodes[a].count+nodes[b].count;if constexpr(Products){nodes[n].ascending=merge(nodes[a].ascending,nodes[b].ascending);nodes[n].descending=merge(nodes[b].descending,nodes[a].descending);}}
 Int branch(Int a,Int b,Int bit){Int n=newNode();nodes[n].left=a;nodes[n].right=b;nodes[n].bit=bit;pull(n);return n;}
 Int singleton(Int key,T value){Int n=newNode();nodes[n].count=1;nodes[n].key=key;nodes[n].bit=-1;if constexpr(Products)nodes[n].ascending=nodes[n].descending=value;return n;}
 std::pair<Int,Int> split(Int node,Int k){if(!k)return {0,node};if(k==nodes[node].count)return {node,0};Int l=nodes[node].left,r=nodes[node].right,c=nodes[l].count;if(k==c){releaseNode(node);return {l,r};}if(k<c){auto [a,b]=split(l,k);nodes[node].left=b;pull(node);return {a,node};}auto [a,b]=split(r,k-c);nodes[node].right=a;pull(node);return {node,b};}
 Int meld(Int a,Int b){Int ab=nodes[a].bit,bb=nodes[b].bit;UInt diff=UInt(nodes[a].key)^UInt(nodes[b].key);Int bit=diff?Int(std::bit_width(diff))-1:-1;if(bit>std::max(ab,bb))return nodes[a].key<nodes[b].key?branch(a,b,bit):branch(b,a,bit);if(ab<bb)return meld(b,a);if(ab>bb){if((UInt(nodes[b].key)>>ab)&1){Int child=meld(nodes[a].right,b);nodes[a].right=child;}else{Int child=meld(nodes[a].left,b);nodes[a].left=child;}}else{assert(ab>=0);Int l=meld(nodes[a].left,nodes[b].left),r=meld(nodes[a].right,nodes[b].right);nodes[a].left=l;nodes[a].right=r;releaseNode(b);}pull(a);return a;}
 Int blockStart(Int index)const{if(roots[index])return index;for(Int n=index+size;n>1;n>>=1)if((n&1)&&starts[n-1]!=-1)return starts[n-1];return -1;}
 void setStart(Int i,Int v){Int n=i+size;starts[n]=v;while((n>>=1)>0){Int value=std::max(starts[n*2],starts[n*2+1]);if(starts[n]==value)break;starts[n]=value;}}
 void setProduct(Int i,T v){if constexpr(Products){Int n=i+size;products[n]=v;while((n>>=1)>0)products[n]=merge(products[n*2],products[n*2+1]);}}
 void refresh(Int i){if constexpr(Products){auto& n=nodes[roots[i]];setProduct(i,reversed[i]?n.descending:n.ascending);}}
 T product(Int l,Int r)const{T a=e,b=e;for(l+=size,r+=size;l<r;l>>=1,r>>=1){if(l&1)a=merge(a,products[l++]);if(r&1)b=merge(products[--r],b);}return merge(a,b);}
 void cut(Int i){if(i==len()||roots[i])return;Int start=blockStart(i),root=roots[start],k=i-start;auto [a,b]=split(root,reversed[start]?nodes[root].count-k:k);if(reversed[start])std::swap(a,b);roots[start]=a;roots[i]=b;reversed[i]=reversed[start];next[i]=next[start];next[start]=i;setStart(i,i);refresh(start);refresh(i);}
 Int locate(Int i)const{assert(0<=i&&i<len());Int start=blockStart(i),node=roots[start],k=i-start;if(reversed[start])k=nodes[node].count-1-k;while(nodes[node].count>1){Int l=nodes[node].left;if(k<nodes[l].count)node=l;else{k-=nodes[l].count;node=nodes[node].right;}}return node;}
 void updateValue(Int node,Int key,T value){if(nodes[node].count==1){nodes[node].ascending=nodes[node].descending=value;return;}updateValue((UInt(key)>>nodes[node].bit)&1?nodes[node].right:nodes[node].left,key,value);pull(node);}
 T nodeProduct(Int node,Int l,Int r,bool rev)const{if(!l&&r==nodes[node].count)return rev?nodes[node].descending:nodes[node].ascending;Int a=nodes[node].left,b=nodes[node].right,c=nodes[a].count;if(r<=c)return nodeProduct(a,l,r,rev);if(l>=c)return nodeProduct(b,l-c,r-c,rev);T x=nodeProduct(a,l,c,rev),y=nodeProduct(b,0,r-c,rev);return rev?merge(y,x):merge(x,y);}
 T blockProduct(Int start,Int l,Int r)const{return reversed[start]?nodeProduct(roots[start],next[start]-r,next[start]-l,true):nodeProduct(roots[start],l-start,r-start,false);}
public:
 template<class Op> RangeSortEngine(std::span<const Int> keys,std::span<const T> values,Int limit,Op op,T identity):roots(keys.size()),next(keys.size()),reversed(keys.size()),keyLimit(limit),merge(op),e(identity){assert(limit>=0);if constexpr(Products)assert(keys.size()==values.size());while(size<len())size*=2;starts.assign(size*2,-1);if constexpr(Products)products.assign(size*2,e);nodes.reserve(std::max<Int>(1,len()*2));nodes.emplace_back();for(Int i=0;i<len();++i){Int key=keys[i];assert(0<=key&&key<limit);
#ifndef NDEBUG
 assert(activeKeys.insert(key).second);
#endif
 roots[i]=singleton(key,Products?values[i]:e);next[i]=i+1;starts[size+i]=i;if constexpr(Products)products[size+i]=values[i];}for(Int i=size-1;i>0;--i){starts[i]=std::max(starts[i*2],starts[i*2+1]);if constexpr(Products)products[i]=merge(products[i*2],products[i*2+1]);}}
 Int len()const{return roots.size();}Int key(Int i)const{return nodes[locate(i)].key;}
 T value(Int i)const requires Products{return nodes[locate(i)].ascending;}
 void update(Int i,T value) requires Products{Int k=key(i),start=blockStart(i);updateValue(roots[start],k,value);refresh(start);}
 void update(Int i,Int k,T value){assert(0<=i&&i<len()&&0<=k&&k<keyLimit);
#ifndef NDEBUG
 Int old=key(i);assert(k==old||!activeKeys.contains(k));activeKeys.erase(old);activeKeys.insert(k);
#endif
 cut(i);cut(i+1);releaseNode(roots[i]);roots[i]=singleton(k,value);reversed[i]=false;refresh(i);}
 T get(Int l,Int r)const requires Products{assert(0<=l&&l<=r&&r<=len());if(l==r)return e;Int a=blockStart(l),b=blockStart(r-1);if(a==b)return blockProduct(a,l,r);T out=blockProduct(a,l,next[a]);if(next[a]<b)out=merge(out,product(next[a],b));return merge(out,blockProduct(b,b,r));}
 T get_all()const requires Products{return products[1];}
 void sort(Int l,Int r,SortOrder order=Ascending){assert(0<=l&&l<=r&&r<=len());if(r-l<=1)return;cut(l);cut(r);Int root=roots[l];for(Int start=next[l];start<r;start=next[start]){root=meld(root,roots[start]);roots[start]=0;setProduct(start,e);setStart(start,-1);}roots[l]=root;next[l]=r;reversed[l]=order==Descending;refresh(l);}
 Generator<Int> items()const{std::vector<Int> stack;for(Int start=0;start<len();start=next[start]){stack.push_back(roots[start]);while(!stack.empty()){Int node=stack.back();stack.pop_back();if(nodes[node].count==1)co_yield nodes[node].key;else if(reversed[start]){stack.push_back(nodes[node].left);stack.push_back(nodes[node].right);}else{stack.push_back(nodes[node].right);stack.push_back(nodes[node].left);}}}}
};
}
}
