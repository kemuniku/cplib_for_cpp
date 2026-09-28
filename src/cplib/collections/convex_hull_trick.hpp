#pragma once
#include <cplib/collections/avltreenode.hpp>
#include <cplib/collections/private/convex_hull_trick_impl.hpp>
namespace cplib {
class ConvexHullTrick {
    AvlNodePtr<CHTLine> root_;
    static void set_start(AvlNodePtr<CHTLine> n,AvlNodePtr<CHTLine> previous){if(n)n->key.start=previous?chtStart(previous->key,n->key):-(__int128(1)<<65);}
    static auto clone(AvlNodePtr<CHTLine> n,AvlNodePtr<CHTLine> p={}) -> AvlNodePtr<CHTLine>{if(!n)return {};auto v=std::make_shared<AvlTreeNode<CHTLine>>(n->key);v->h=n->h;v->len=n->len;v->p=p;v->l=clone(n->l,v);v->r=clone(n->r,v);return v;}
public:
    ConvexHullTrick()=default;ConvexHullTrick(const ConvexHullTrick& v):root_(clone(v.root_)){}ConvexHullTrick& operator=(const ConvexHullTrick& v){if(this!=&v)root_=clone(v.root_);return *this;}ConvexHullTrick(ConvexHullTrick&&)=default;ConvexHullTrick& operator=(ConvexHullTrick&&)=default;
    // 任意の傾きを追加し、不要な直線を AVL 木から除く。償却 O(log N)。
    void add_line(Int a,Int b){CHTLine line{a,b};auto [left,right]=lower_bound_node(root_,line);if(right&&right->key.a==a){if(right->key.b<=b)return;auto n=next(right);root_=erase(root_,right,n);right=n;}if(left&&right&&chtRedundant(left->key,line,right->key))return;
        while(left){auto p=prev(left);if(!p||!chtRedundant(p->key,left->key,line))break;root_=erase(root_,left,right);left=p;}
        while(right){auto n=next(right);if(!n||!chtRedundant(line,right->key,n->key))break;root_=erase(root_,right,n);right=n;}
        auto n=std::make_shared<AvlTreeNode<CHTLine>>(line);set_start(n,left);set_start(right,n);root_=insert(root_,n);
    }
    // 最小となる区間の左端で探索。O(log N)。
    Int get_min(Int x)const{assert(root_);auto n=root_;AvlNodePtr<CHTLine> best;while(n){if(n->key.start<=__int128(x)){best=n;n=n->r;}else n=n->l;}return chtAnswer(chtValue(best->key,x));}
};
inline ConvexHullTrick initConvexHullTrick(){return {};}
}
