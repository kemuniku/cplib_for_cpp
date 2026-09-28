#pragma once
#include <cplib/common.hpp>
#include <memory>
namespace cplib {
// 旧版の番兵ノード付き AVL 木。現行版とは独立した高さ管理を保持。
template<class K> struct AvlTreeNode;
template<class K> std::shared_ptr<AvlTreeNode<K>> get_avltree_nilnode();
template<class K> struct AvlTreeNode {
    using Ptr=std::shared_ptr<AvlTreeNode>;
    Ptr l,r;std::weak_ptr<AvlTreeNode> p;Int h=1,len=1;K key;
    explicit AvlTreeNode(K value,bool sentinel=false):key(std::move(value)){if(!sentinel){l=r=get_avltree_nilnode<K>();p=l;}}
};
template<class K> using AvlNodePtr=std::shared_ptr<AvlTreeNode<K>>;
template<class K> AvlNodePtr<K> get_avltree_nilnode(){
    static AvlTreeNode<K> node(K{},true);
    static AvlNodePtr<K> ptr=[](){AvlNodePtr<K> out(&node,[](auto*){});node.h=node.len=0;node.l=node.r=out;node.p=out;return out;}();
    return ptr;
}
namespace detail {
template<class K> void avl_update(AvlNodePtr<K> n){n->h=std::max(n->l->h,n->r->h)+1;n->len=n->l->len+n->r->len+1;}
template<class K> void avl_children(AvlNodePtr<K> n,AvlNodePtr<K> l,AvlNodePtr<K> r){n->l=l;if(l!=get_avltree_nilnode<K>())l->p=n;n->r=r;if(r!=get_avltree_nilnode<K>())r->p=n;avl_update(n);}
template<class K> AvlNodePtr<K> avl_rebalance(AvlNodePtr<K> n){auto l=n->l,r=n->r;Int lh=l?l->h:0,rh=r?r->h:0;
    if(lh+1<rh){auto rl=r->l,rr=r->r;if((rl?rl->h:0)<=(rr?rr->h:0)){r->p=n->p;avl_children(n,l,rl);avl_children(r,n,rr);return r;}else{rl->p=n->p;auto rll=rl->l,rlr=rl->r;avl_children(n,l,rll);avl_children(r,rlr,rr);avl_children(rl,n,r);return rl;}}
    if(rh+1<lh){auto ll=l->l,lr=l->r;if((lr?lr->h:0)<=(ll?ll->h:0)){l->p=n->p;avl_children(n,lr,r);avl_children(l,ll,n);return l;}else{lr->p=n->p;auto lrl=lr->l,lrr=lr->r;avl_children(n,lrr,r);avl_children(l,ll,lrl);avl_children(lr,l,n);return lr;}}
    avl_update(n);return n;
}
template<class K> AvlNodePtr<K> avl_to_root(AvlNodePtr<K> n){while(n->p.lock()!=get_avltree_nilnode<K>()){auto p=n->p.lock();if(p->l==n)p->l=avl_rebalance(n);else p->r=avl_rebalance(n);n=p;}return avl_rebalance(n);}
template<class K> auto avl_search(AvlNodePtr<K> n,const K& key,bool strict){AvlNodePtr<K> l=get_avltree_nilnode<K>(),r=l;while(n!=get_avltree_nilnode<K>()){if(strict?key<n->key:!(n->key<key)){r=n;n=n->l;}else{l=n;n=n->r;}}return std::pair{l,r};}
}
template<class K> AvlNodePtr<K> rootOf(AvlNodePtr<K> n){while(n->p.lock()!=get_avltree_nilnode<K>())n=n->p.lock();return n;}
template<class K> auto lower_bound_node(AvlNodePtr<K> n,const K& key){return detail::avl_search(n,key,false);}
template<class K> auto upper_bound_node(AvlNodePtr<K> n,const K& key){return detail::avl_search(n,key,true);}
template<class K> AvlNodePtr<K> insert(AvlNodePtr<K> n,AvlNodePtr<K> x){if(n==get_avltree_nilnode<K>())return x;auto [l,r]=lower_bound_node(n,x->key);if(l!=get_avltree_nilnode<K>()&&l->r==get_avltree_nilnode<K>()){detail::avl_children(l,l->l,x);return detail::avl_to_root(l);}detail::avl_children(r,x,r->r);return detail::avl_to_root(r);}
template<class K> AvlNodePtr<K> erase(AvlNodePtr<K>,AvlNodePtr<K> x,AvlNodePtr<K> nxt){auto xp=x->p.lock();AvlNodePtr<K> result=get_avltree_nilnode<K>();
    if(x->r==get_avltree_nilnode<K>()){auto xl=x->l;if(xl!=get_avltree_nilnode<K>())xl->p=xp;if(xp!=get_avltree_nilnode<K>()){if(xp->l==x)xp->l=xl;else xp->r=xl;}result=xp!=get_avltree_nilnode<K>()?detail::avl_to_root(xp):xl;}
    else{auto nxtp=nxt->p.lock(),nxtr=nxt->r;if(xp!=get_avltree_nilnode<K>()){if(xp->l==x)xp->l=nxt;else xp->r=nxt;}nxt->p=xp;nxt->l=x->l;if(nxt->l!=get_avltree_nilnode<K>())nxt->l->p=nxt;if(x->r==nxt){detail::avl_update(nxt);result=detail::avl_to_root(nxt);}else{if(nxtp->l==nxt)nxtp->l=nxtr;else nxtp->r=nxtr;if(nxtr!=get_avltree_nilnode<K>())nxtr->p=nxtp;nxt->r=x->r;nxt->r->p=nxt;detail::avl_update(nxt);result=detail::avl_to_root(nxtp);}}
    x->l=x->r=get_avltree_nilnode<K>();x->p=x->l;detail::avl_update(x);return result;
}
template<class K> AvlNodePtr<K> next(AvlNodePtr<K> n){if(n->r!=get_avltree_nilnode<K>()){n=n->r;while(n->l!=get_avltree_nilnode<K>())n=n->l;return n;}while(n->p.lock()!=get_avltree_nilnode<K>()){auto p=n->p.lock();if(p->r!=n)return p;n=p;}return get_avltree_nilnode<K>();}
template<class K> AvlNodePtr<K> prev(AvlNodePtr<K> n){if(n->l!=get_avltree_nilnode<K>()){n=n->l;while(n->r!=get_avltree_nilnode<K>())n=n->r;return n;}while(n->p.lock()!=get_avltree_nilnode<K>()){auto p=n->p.lock();if(p->l!=n)return p;n=p;}return get_avltree_nilnode<K>();}
template<class K> AvlNodePtr<K> get(AvlNodePtr<K> n,Int idx){assert(idx>=0);if(idx>=n->len)return get_avltree_nilnode<K>();for(;;){Int left=n->l?n->l->len:0;if(left==idx)return n;if(left<idx){idx-=left+1;n=n->r;}else n=n->l;}}
template<class K> Int index(AvlNodePtr<K> n){if(n==get_avltree_nilnode<K>())return 0;Int result=n->l?n->l->len:0;while(n->p.lock()!=get_avltree_nilnode<K>()){auto p=n->p.lock();if(p->r==n)result+=1+(p->l?p->l->len:0);n=p;}return result;}
}
