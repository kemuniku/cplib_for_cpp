#pragma once
#include <cplib/str/run_length_encode.hpp>
#include <unordered_map>
namespace cplib {
class PrimeFactorTable {
    std::vector<Int> table;
public:
    // 最大素因数表を構築する。O(n log log n)時間・O(n)領域。
    explicit PrimeFactorTable(Int maxn):table(maxn+1) {
        for(Int i=2;i<=maxn;++i) if(!table[i]) for(Int j=i;j<=maxn;j+=i) table[j]=i;
    }
    // 素因数を昇順に返す。O(log x)。
    std::vector<Int> primefactor(Int x) const {
        assert(x>=1 && static_cast<std::size_t>(x)<table.size());
        std::vector<Int> result;
        while(x!=1){result.push_back(table[x]);x/=table[x];}
        std::reverse(result.begin(),result.end());return result;
    }
    // 素因数ごとの指数を返す。O(log x)。
    std::unordered_map<Int,Int> primefactor_table(Int x) const {
        std::unordered_map<Int,Int> result; for(auto p:primefactor(x)) ++result[p];return result;
    }
    // 素因数と指数の組を昇順で返す。O(log x)。
    auto primefactor_tuple(Int x) const {return run_length_encode(primefactor(x));}
};
// 最大素因数表を構築する。
inline PrimeFactorTable initPrimeFactorTable(Int maxn){return PrimeFactorTable(maxn);}
inline auto primefactor(const PrimeFactorTable& t,Int x){return t.primefactor(x);}
inline auto primefactor_table(const PrimeFactorTable& t,Int x){return t.primefactor_table(x);}
inline auto primefactor_tuple(const PrimeFactorTable& t,Int x){return t.primefactor_tuple(x);}
}
