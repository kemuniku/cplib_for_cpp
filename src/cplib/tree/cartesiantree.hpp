#pragma once
#include <cplib/common.hpp>
namespace cplib {
struct CartesianTreeNode{Int p=-1,l=-1,r=-1;};
// 最小値のCartesian treeを親・左右の子で返す。同値は左を優先。O(n)。
inline std::vector<CartesianTreeNode> cartesian_tree_tuple(std::span<const Int> a){
    std::vector<CartesianTreeNode> result(a.size());std::vector<Int> stack;
    for(Int i=0;i<Int(a.size());++i){Int s=-1;while(!stack.empty()&&a[stack.back()]>a[i]){s=stack.back();stack.pop_back();}if(s!=-1){if(result[i].l!=-1)std::swap(result[i].l,result[s].l);else result[i].l=s;result[s].p=i;}if(!stack.empty()){result[i].p=stack.back();result[stack.back()].r=i;}stack.push_back(i);}return result;
}
}
