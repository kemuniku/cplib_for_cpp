#pragma once
#include <cplib/collections/private/scalar_bitset.hpp>
namespace cplib {
using BitSet=ScalarBitSet<-1>;
inline BitSet initBitSet(Int n){return BitSet(n);}inline BitSet initBitSet(const std::vector<bool>& v,Int n){return BitSet(v,n);}inline BitSet initBitSet(const std::vector<bool>& v){return BitSet(v);}
inline BitSet initBitSetFromIndexes(const std::vector<Int>& indexes,Int n){BitSet out(n);for(Int i:indexes)out.set(i,true);return out;}
}
