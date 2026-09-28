#pragma once
#include <cplib/collections/persistent_lazysegtree.hpp>
namespace cplib {
template<class T=Int> struct SegmentTreeNode{T value;std::shared_ptr<SegmentTreeNode> left,right;};
namespace detail {struct PersistentNoAction{};}
template<class T=Int> class PSegmentTree {
 PersistentLazySegmentTree<T,detail::PersistentNoAction> tree;
 explicit PSegmentTree(decltype(tree) t):tree(std::move(t)){}
public:
 template<class V,class Op> PSegmentTree(const V& v,Op op,T e):tree(v,op,e,[](auto,T x){return x;},[](auto,auto){return detail::PersistentNoAction{};},detail::PersistentNoAction{}){}
 Int len()const{return tree.len();}PSegmentTree update(Int p,T v)const{return PSegmentTree(tree.update(p,std::move(v)));}
 PSegmentTree copy_range(const PSegmentTree& src,Int l,Int r)const{return PSegmentTree(tree.copy_range(src.tree,l,r));}
 template<class L,class R> PSegmentTree copy_range(const PSegmentTree& src,ClosedSlice<L,R> s)const{return PSegmentTree(tree.copy_range(src.tree,s));}
 T get(Int l,Int r)const{return tree.get(l,r);}T query(Int l,Int r)const{return get(l,r);}T get_all()const{return tree.get_all();}
 template<class I> T operator[](I i)const{return tree[i];}template<class L,class R> T get(ClosedSlice<L,R> s)const{return tree.get(s);}
 template<class I> void set(I i,T v){tree.set(i,std::move(v));}
 template<class Pred> Int max_right(Int l,Pred f)const{return tree.max_right(l,f);}template<class Pred> Int min_left(Int r,Pred f)const{return tree.min_left(r,f);}
 std::string str()const{return tree.str();}
};
template<class V,class Op,class T> PSegmentTree(const V&,Op,T)->PSegmentTree<T>;
template<class T> using PersistentSegmentTree=PSegmentTree<T>;
template<class V,class Op,class T> auto initPersistentSegmentTree(const V& v,Op op,T e){return PSegmentTree<T>(v,op,e);}
template<class V,class Op,class T> auto initSegmentTree(const V& v,Op op,T e){return initPersistentSegmentTree(v,op,e);}
template<class T> auto update(const PSegmentTree<T>& s,Int i,T v){return s.update(i,v);}template<class T> T query(const PSegmentTree<T>& s,Int l,Int r){return s.query(l,r);}
template<class... A> auto newPersistentSegWith(A&&... args){return initPersistentSegmentTree(std::forward<A>(args)...);}
template<class... A> auto newSegWith(A&&... args){return initPersistentSegmentTree(std::forward<A>(args)...);}
}
