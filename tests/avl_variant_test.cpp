#include CPLIB_AVL_HEADER
#include <random>
#include <set>
using namespace cplib;
void avl_version_test(){
    std::mt19937 rng(19);auto m=initAvlSortedMultiSet<Int>();auto s=initAvlSortedSet<Int>();std::multiset<Int> expected;std::set<Int> unique;
    for(int q=0;q<30000;++q){Int x=Int(rng()%1000)-500;if(rng()%2){m.incl(x);s.incl(x);expected.insert(x);unique.insert(x);}else{auto it=expected.find(x);bool has=it!=expected.end();assert(m.excl(x)==has);if(has)expected.erase(it);assert(s.excl(x)==bool(unique.erase(x)));}assert(m.len()==Int(expected.size())&&s.len()==Int(unique.size()));assert(m.count(x)==Int(expected.count(x)));assert(m.lowerBound(x)==std::distance(expected.begin(),expected.lower_bound(x)));assert(m.upperBound(x)==std::distance(expected.begin(),expected.upper_bound(x)));if(q%100==0){std::vector<Int> v;for(Int y:m)v.push_back(y);assert(std::equal(v.begin(),v.end(),expected.begin(),expected.end()));for(Int i=0;i<m.len();++i){auto n=get(m.root,i);assert(n->key==v[i]&&index(n)==i&&rootOf(n)==m.root);assert(m[i]==v[i]);if(i+1<m.len())assert(next(n)->key==v[i+1]);if(i)assert(prev(n)->key==v[i-1]);}auto copy=m;copy.incl(10000);assert(!m.contains(10000));}}
    while(m.len()){assert(m.pop()==*expected.rbegin());expected.erase(std::prev(expected.end()));}

}

int main(){avl_version_test();}
