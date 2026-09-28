#pragma once
#include <cplib/common.hpp>
#include <cplib/utils/backwards_index.hpp>
#include <memory>
namespace cplib {
// 区間をノードに持つ遅延AVL木。構築O(1)、各操作O(log(K+2))、Q回更新後O(Q+1)領域。
template<class S,class F> class DynamicLazySegmentTree {
    struct Node;using Ptr=std::shared_ptr<Node>;
    struct Node{Int a,b,lo,hi,height=1;S value,product;F tag,lazy;bool pending=false;Ptr left,right;};
    Ptr root;Int length=0,nodes=0;std::function<S(S,S)> merge;S identity;std::function<S(F,S)> mapping;std::function<F(F,F)> composition;F id;std::function<S(Int,Int)> initial;
    Ptr makeNode(Int a,Int b,F tag){S value=mapping(tag,initial(a,b));++nodes;return std::make_shared<Node>(Node{a,b,a,b,1,value,value,tag,id,false,nullptr,nullptr});}
    static Int height(const Ptr& p){return p?p->height:0;}
    static Ptr clone(const Ptr& p){if(!p)return nullptr;auto q=std::make_shared<Node>(*p);q->left=clone(p->left);q->right=clone(p->right);return q;}
    void pull(const Ptr& p){p->product=p->value;p->lo=p->a;p->hi=p->b;if(p->left){p->product=merge(p->left->product,p->product);p->lo=p->left->lo;}if(p->right){p->product=merge(p->product,p->right->product);p->hi=p->right->hi;}p->height=std::max(height(p->left),height(p->right))+1;}
    void allApply(const Ptr& p,F f){if(!p)return;p->value=mapping(f,p->value);p->product=mapping(f,p->product);p->tag=composition(f,p->tag);if(p->pending)p->lazy=composition(f,p->lazy);else{p->lazy=f;p->pending=true;}}
    void push(const Ptr& p){if(p->pending){allApply(p->left,p->lazy);allApply(p->right,p->lazy);p->lazy=id;p->pending=false;}}
    Ptr rotateLeft(Ptr p){push(p);auto result=p->right;push(result);p->right=result->left;result->left=p;pull(p);pull(result);return result;}
    Ptr rotateRight(Ptr p){push(p);auto result=p->left;push(result);p->left=result->right;result->right=p;pull(p);pull(result);return result;}
    Ptr balance(Ptr p){pull(p);if(height(p->left)>height(p->right)+1){if(height(p->left->left)<height(p->left->right))p->left=rotateLeft(p->left);return rotateRight(p);}if(height(p->right)>height(p->left)+1){if(height(p->right->right)<height(p->right->left))p->right=rotateRight(p->right);return rotateLeft(p);}return p;}
    Ptr insertFirst(Ptr p,Ptr added){if(!p)return added;push(p);p->left=insertFirst(p->left,added);return balance(p);}
    Ptr splitAt(Ptr p,Int index){if(!p||index==p->a||index==p->b)return p;push(p);if(index<p->a)p->left=splitAt(p->left,index);else if(p->b<index)p->right=splitAt(p->right,index);else{auto added=makeNode(index,p->b,p->tag);p->b=index;p->value=mapping(p->tag,initial(p->a,index));p->right=insertFirst(p->right,added);}return balance(p);}
    void applyNode(const Ptr& p,Int l,Int r,F f){if(!p||r<=p->lo||p->hi<=l)return;if(l<=p->lo&&p->hi<=r){allApply(p,f);return;}push(p);applyNode(p->left,l,r,f);if(l<=p->a&&p->b<=r){p->value=mapping(f,p->value);p->tag=composition(f,p->tag);}applyNode(p->right,l,r,f);pull(p);}
    S getNode(const Ptr& p,Int l,Int r){if(!p||r<=p->lo||p->hi<=l)return identity;if(l<=p->lo&&p->hi<=r)return p->product;push(p);S result=getNode(p->left,l,r);Int a=std::max(l,p->a),b=std::min(r,p->b);if(a<b){S value=a==p->a&&b==p->b?p->value:mapping(p->tag,initial(a,b));result=merge(result,value);}return merge(result,getNode(p->right,l,r));}
    void setNode(const Ptr& p,Int index,S value){push(p);if(index<p->a)setNode(p->left,index,value);else if(index>=p->b)setNode(p->right,index,value);else{p->value=value;p->tag=id;}pull(p);}
public:
    template<class Op,class Map,class Comp,class Initial> DynamicLazySegmentTree(Int n,Op op,S e,Map map,Comp comp,F id_,Initial init):length(n),merge(op),identity(e),mapping(map),composition(comp),id(id_),initial(init){assert(n>=0);if(n>0)root=makeNode(0,n,id);}
    DynamicLazySegmentTree(const DynamicLazySegmentTree& p):root(clone(p.root)),length(p.length),nodes(p.nodes),merge(p.merge),identity(p.identity),mapping(p.mapping),composition(p.composition),id(p.id),initial(p.initial){}
    DynamicLazySegmentTree(DynamicLazySegmentTree&&)=default;
    DynamicLazySegmentTree& operator=(DynamicLazySegmentTree&&)=default;
    DynamicLazySegmentTree& operator=(const DynamicLazySegmentTree& p){if(this!=&p){DynamicLazySegmentTree copy(p);*this=std::move(copy);}return *this;}
    // 更新は両端の境界を作る。追加ノードは高々2個。
    void apply(Int l,Int r,F f){assert(0<=l&&l<=r&&r<=length);if(l==r)return;root=splitAt(root,l);root=splitAt(root,r);applyNode(root,l,r,f);}
    S get(Int l,Int r){assert(0<=l&&l<=r&&r<=length);return l==r?identity:getNode(root,l,r);}
    void update(Int p,S value){assert(0<=p&&p<length);root=splitAt(root,p);root=splitAt(root,p+1);setNode(root,p,value);}
    template<class L,class R> S get(ClosedSlice<L,R> s){return get(resolve_index(length,s.a),resolve_index(length,s.b)+1);}
    template<class L,class R> void apply(ClosedSlice<L,R> s,F f){apply(resolve_index(length,s.a),resolve_index(length,s.b)+1,f);}
    template<class L,class R> S operator[](ClosedSlice<L,R> s){return get(s);}
    struct Reference{DynamicLazySegmentTree* tree;Int index;operator S()const{return tree->get(index,index+1);}Reference& operator=(S value){tree->update(index,value);return *this;}Reference& operator=(const Reference& r){return *this=S(r);}};
    Reference operator[](Int p){assert(0<=p&&p<length);return {this,p};}
    Reference operator[](BackwardsIndex p){return (*this)[length-p.value];}
    Int len()const{return length;}
    Int node_count()const{return nodes;}
    S get_all()const{return root?root->product:identity;}
};
template<class S,class F,class Op,class Map,class Comp,class Initial> auto initDynamicLazySegmentTree(Int n,Op op,S e,Map map,Comp comp,F id,Initial initial){return DynamicLazySegmentTree<S,F>(n,op,e,map,comp,id,initial);}
template<class... Args> auto newDynamicLazySegWith(Args&&... args){return initDynamicLazySegmentTree(std::forward<Args>(args)...);}
template<class S,class F> Int len(const DynamicLazySegmentTree<S,F>& t){return t.len();}
template<class S,class F> Int node_count(const DynamicLazySegmentTree<S,F>& t){return t.node_count();}
template<class S,class F> S get_all(const DynamicLazySegmentTree<S,F>& t){return t.get_all();}
template<class S,class F> S get(DynamicLazySegmentTree<S,F>& t,Int l,Int r){return t.get(l,r);}
template<class S,class F,class L,class R> S get(DynamicLazySegmentTree<S,F>& t,ClosedSlice<L,R> s){return t.get(s);}
template<class S,class F> void apply(DynamicLazySegmentTree<S,F>& t,Int l,Int r,F f){t.apply(l,r,f);}
template<class S,class F,class L,class R> void apply(DynamicLazySegmentTree<S,F>& t,ClosedSlice<L,R> s,F f){t.apply(s,f);}
template<class S,class F> void update(DynamicLazySegmentTree<S,F>& t,Int p,S value){t.update(p,value);}
}
