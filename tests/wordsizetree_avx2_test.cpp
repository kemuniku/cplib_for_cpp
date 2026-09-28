#include <cplib/collections/wordsizetree_avx2.hpp>
#include <random>
#include <set>
#include <memory>
using namespace cplib;
void wordsize_test(){
    std::mt19937 rng(923);constexpr Int U=WordsizeTreeAvx2Capacity;
    for(Int n:{Int(0),Int(63),Int(65),Int(256),Int(10000),U}){std::vector<bool> bits(n);std::string str(n,'0');std::set<Int> values;for(Int i=0;i<n;++i)if(rng()%1000==0){bits[i]=true;str[i]='1';values.insert(i);}auto avx=std::make_unique<WordsizeTreeAvx2>(std::string_view(str));auto avxBool=std::make_unique<WordsizeTreeAvx2>(bits);for(int q=0;q<10000;++q){Int x=rng()%U;if(q<6)x=std::array<Int,6>{0,63,64,255,256,U-1}[q];if(rng()%2){avx->incl(x);avxBool->incl(x);values.insert(x);}else{avx->excl(x);avxBool->excl(x);values.erase(x);}assert((*avx)[x]==values.contains(x)&&(*avxBool)[x]==values.contains(x));Int y=rng()%U;if(q%10==0)y=x;auto lb=values.lower_bound(y),ub=values.upper_bound(y);Int ge=lb==values.end()?-1:*lb,le=ub==values.begin()?-1:*std::prev(ub);assert(avx->ge(y)==ge&&avxBool->ge(y)==ge);assert(avx->le(y)==le&&avxBool->le(y)==le);}for(Int x:values){avx->excl(x);}assert(avx->ge(-1)==-1&&avx->le(U)==-1);assert(avx->le(-1)==-1&&avx->ge(U)==-1);}
}
int main(){wordsize_test();}
