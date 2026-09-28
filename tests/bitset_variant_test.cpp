#include CPLIB_BITSET_HEADER
#include <random>
using namespace cplib;

template<Int N> void check() {
    std::mt19937 rng(N);
    std::vector<bool> values(N);
    std::vector<Int> indices;
    for (Int i=0;i<N;++i) if ((values[i]=rng()%2)) indices.push_back(i);
#ifdef CPLIB_STATIC_BITSET
    cplib::BitSet<N> bits=cplib::initBitSet<N>(values);
    auto from_indices=cplib::initBitSetFromIndexes<N>(indices);
#else
    cplib::BitSet bits=cplib::initBitSet(values);
    auto from_indices=cplib::initBitSetFromIndexes(indices,N);
#endif
    assert(bits==from_indices);
    assert(bits.len()==N && bits.popcount()==Int(indices.size()));
    for (Int i=0;i<N;++i) assert(bool(bits[i])==values[i]);
}
int main(){check<0>();check<1>();check<63>();check<64>();check<65>();check<257>();check<1003>();}
