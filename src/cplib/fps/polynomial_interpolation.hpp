#pragma once
#include <cplib/fps/product_tree.hpp>
namespace cplib {
// 積木と導関数の多点評価から補間する。各xは異なること。O(M(N)log N)。
template<Modint T> std::vector<T> polynomialInterpolation(const std::vector<T>& xs,const std::vector<T>& ys){assert(xs.size()==ys.size());if(xs.empty())return {};auto tree=initPolynomialProductTree(xs);auto denominators=multipointEvaluation(derivative(tree.nodes[1]),xs);std::vector<std::vector<T>> partial(tree.nodes.size());for(Int i=0;i<tree.leafCount;++i)if(i<Int(xs.size())){assert(denominators[i].val()!=0);partial[tree.leafCount+i]={ys[i]/denominators[i]};}for(Int i=tree.leafCount-1;i>=1;--i)partial[i]=partial[i*2]*tree.nodes[i*2+1]+partial[i*2+1]*tree.nodes[i*2];return prefix(partial[1],xs.size());}
}
