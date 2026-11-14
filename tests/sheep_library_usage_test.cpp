#include "sheep_library_macros.hpp"
#include <cplib/collections/bitset.hpp>
#include <cplib/collections/bitvector.hpp>
#include <cplib/collections/fenwick.hpp>
#include <cplib/collections/segtree.hpp>
#include <cplib/collections/QSWAG.hpp>
#include <cplib/geometry/convex_polygon.hpp>
#include <cplib/graph/dijkstra.hpp>
#include <cplib/str/hash_string.hpp>
#include <cplib/str/repeated_static_string.hpp>
#include <cplib/str/trie.hpp>
#include <cplib/str/lcs_bitset.hpp>
#include <cplib/matrix/matrix.hpp>
using namespace cplib;

#ifndef len
#error sheep len must remain a macro
#endif
#ifndef out
#error sheep out must remain a macro
#endif
#ifndef read
#error sheep read must remain a macro
#endif

int main() {
    vector<Int> values{3, 1, 2};
    sort(all(values));
    assert(values == vector<Int>({1, 2, 3}));
    auto fenwick = initFenwickTree(values);
    assert((cplib::len)(fenwick) == 3);
    assert(fenwick.size() == 3 && len(fenwick) == 3);
    assert(fenwick.get(0, 3) == 6);
    auto segment = initSegmentTree(values, std::plus<Int>{}, Int(0));
    assert(segment.size() == 3 && len(segment) == 3);
    segment.update(1, 5);
    assert(segment.get(0, 3) == 9);
    auto queue = initQSWAG(std::plus<Int>{}, Int(0));
    queue.push(2); queue.push(3);
    assert(queue.size() == 2 && len(queue) == 2 && queue.fold() == 5);

    auto bits = initBitSet(vector<bool>{true, true});
    assert(bits.size() == 2 && len(bits) == 2 && bits.allOnes());
    BitVector bit_vector(3);
    setValue(bit_vector, 1); bit_vector.build();
    assert(bit_vector[1]);
    set<int> ordered{4};
    hash<int> hasher;
    assert(ordered.count(4) == 1 && hasher(4) == std::hash<int>{}(4));
    Fraction<int> fraction(1, 2);
    assert(hashValue(fraction) == std::hash<Fraction<int>>{}(fraction));

    auto rolling = initRollingHash("abc");
    HashString hashed = rolling;
    assert(hashed.size() == 3 && len(hashed) == 3);
    auto text = toStaticString("ab");
    auto repeated = initRepeatedStaticString(text, 5);
    assert(repeated.size() == 5 && len(repeated) == 5 && repeated[4] == 'a');
    Trie<> trie;
    trie.incl("ab");
    assert(trie.size() == 1 && len(trie) == 1);
    assert(cplib::bitset_lcs::LCS(std::string(80, 'a'), std::string(75, 'a')) == 75);
    auto polygon = initPolygon<Int>({{0, 0}, {2, 0}, {0, 2}});
    assert(polygon.size() == 3 && len(polygon) == 3 && area(polygon) == 4);
    auto graph = initUnWeightedDirectedGraph(3);
    assert(graph.size() == 3 && len(graph) == 3);

    std::istringstream input_stream("7 8 abc");
    auto *saved_input = cin.rdbuf(input_stream.rdbuf());
    ii(a, b);
    si(word);
    cin.rdbuf(saved_input);
    assert(a == 7 && b == 8 && word == "abc");
    std::ostringstream output_stream;
    auto *saved_output = cout.rdbuf(output_stream.rdbuf());
    out(a, b);
    cout.rdbuf(saved_output);
    assert(output_stream.str() == "7 8\n");
}
