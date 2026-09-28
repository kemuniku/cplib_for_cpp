#pragma once
#include <cplib/tree/heavylightdecomposition.hpp>
#include <memory>
namespace cplib {
enum StaticTopTreeNodeKind {sttLeaf,sttCompress,sttRake};
struct StaticTopTreeNode {StaticTopTreeNodeKind kind=sttLeaf;Int parent=-1,left=-1,right=-1,upper=-1,lower=-1,size=0;};
struct StaticTopTreeData {
 Int numVertices=0,root=0;std::vector<StaticTopTreeNode> nodes;
 // 累積頂点数の中央に最も近い境界で、順序を保って結合する。
 Int mergeRange(const std::vector<Int>& items,const std::vector<Int>& prefix,Int left,Int right,StaticTopTreeNodeKind kind){if(right-left==1)return items[left];Int weight=prefix[right]-prefix[left],target=prefix[left]+(weight+1)/2,lo=left+1,hi=right-1;while(lo<hi){Int mid=(lo+hi)/2;if(prefix[mid]<target)lo=mid+1;else hi=mid;}Int mid=lo;if(mid>left+1&&std::abs(2*(prefix[mid-1]-prefix[left])-weight)<=std::abs(2*(prefix[mid]-prefix[left])-weight))--mid;
  Int l=mergeRange(items,prefix,left,mid,kind),r=mergeRange(items,prefix,mid,right,kind),upper=nodes[l].upper,lower=kind==sttCompress?nodes[r].lower:nodes[l].lower;
  if(kind==sttCompress)assert(nodes[l].lower==nodes[r].upper);else assert(nodes[l].upper==nodes[r].upper);
  Int out=nodes.size();nodes.push_back({kind,-1,l,r,upper,lower,weight});nodes[l].parent=nodes[r].parent=out;return out;
 }
 // 重み付き平衡結合。O(K log K)。
 Int mergeBalanced(const std::vector<Int>& items,StaticTopTreeNodeKind kind){assert(!items.empty());if(items.size()==1)return items[0];std::vector<Int> prefix(items.size()+1);for(std::size_t i=0;i<items.size();++i)prefix[i+1]=prefix[i]+nodes[items[i]].size;return mergeRange(items,prefix,0,items.size(),kind);}
};
using StaticTopTree=std::shared_ptr<StaticTopTreeData>;
// 非空の木のHLDから高さO(log N)の構造を作る。O(N log N)時間、O(N)空間。
inline StaticTopTree initStaticTopTree(const HeavyLightDecomposition& hld){Int n=hld.numVertices();assert(n>0);auto tree=std::make_shared<StaticTopTreeData>();tree->numVertices=n;tree->nodes.reserve(2*n-1);for(Int v=0;v<n;++v)tree->nodes.push_back({sttLeaf,-1,-1,-1,hld.parentOf(v),v,1});std::vector<Int> pathRoot(n);for(Int i=n-1;i>=0;--i){Int head=hld.toVtx(i);if(hld.heavyRootOf(head)!=head)continue;std::vector<Int> path{head};Int v=head;while(true){Int heavy=hld.heavyChildOf(v);if(heavy==-1)break;std::vector<Int> branches{heavy};for(Int child:hld.children(v))if(child!=heavy)branches.push_back(pathRoot[child]);path.push_back(tree->mergeBalanced(branches,sttRake));v=heavy;}pathRoot[head]=tree->mergeBalanced(path,sttCompress);}tree->root=pathRoot[hld.toVtx(0)];assert(Int(tree->nodes.size())==2*n-1);return tree;}
inline StaticTopTree initStaticTopTreeFromParent(std::span<const Int> parent,Int root=0){return initStaticTopTree(initHldFromParent(parent,root));}
}
