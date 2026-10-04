"""Compile independent size() contracts; variants intentionally have overlapping names."""
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
import os
import shlex
import subprocess
import tempfile

ROOT = Path(__file__).resolve().parents[1]
CASES = []


def case(header, body, extra=""):
    CASES.append((header, body, extra))


for header in ["avlset", "avlset_old", "raw_ptr_avlset", "tatyamset"]:
    case("collections/" + header, """
auto s = initAvlSortedMultiSet<Int>(); check(s, 0);
s.incl(4); check(s, 1); s.incl(4); check(s, 2);
s.excl(4); check(s, 1); s.excl(4); check(s, 0);
""" if header != "tatyamset" else """
auto s = initSortedMultiset<Int>(); check(s, 0);
s.incl(4); check(s, 1); s.incl(4); check(s, 2);
s.excl(4); check(s, 1); s.excl(4); check(s, 0);
""")

case("str/static_string", """
StaticString<char> empty; check(empty, 0);
auto s = toStaticString("abc", true); check(s, 3);
check(s.substr(1, 2), 1); check(s.substr(1, 1), 0);
check(s.reversed(), 3); check(toStaticString(std::vector<Int>{9}), 1);
assert(s.base->size == 3 && s.base->S.size() == 6);
""")
case("str/merged_static_string", """
MergedStaticString<char> s; check(s, 0);
auto text = toStaticString("abc"); s &= text.substr(0, 1); check(s, 1);
s &= text; check(s, 4); check(s[closed_slice(1, 2)], 2);
check(mergeStaticStrings(text, text), 6);
""")
case("str/fixedlength_merged_static_string", """
FixedLengthMergedStaticString<char, 0> e; check(e, 0);
auto s = toStaticString("a"); check(e & s, 1); check(s & s, 2);
""")
case("str/rolling_hash", """
for (const auto& text : {std::string{}, std::string("a"), std::string("abc")}) {
    auto s = initRollingHash(text); check(s, text.size()); s.build(); check(s, text.size());
}
""")
for header in ["hash_string", "can_reverse_hash_string"]:
    case("str/" + header, """
RollingHash e; check(e, 0); auto s = initRollingHash("abc"); check(s, 3);
check(s.substr(1, 1), 0); check(s.substr(1, 2), 1); check(*s.R, 3);
HashString h = s; static_assert(std::is_same_v<decltype(h.size), Int>);
assert(h.size == 3 && h.len() == 3);
""")
case("str/trie", """
Trie<> t; check(t, 0); t.incl("a"); check(t, 1);
t.incl("a"); check(t, 2); t.incl(""); check(t, 3);
t.excl("a"); check(t, 2); t.excl(""); check(t, 1);
""")
case("str/double_ended_palindromic_tree", """
DoubleEndedPalindromicTree t; check(t, 0);
t.push_back('a'); check(t, 1); t.push_front('b'); check(t, 2);
t.pop_back(); check(t, 1); t.pop_front(); check(t, 0);
""")
case("geometry/polygon", """
Polygon<Int> p; check(p, 0); p.v.push_back({1, 2}); check(p, 1);
p.v.push_back({2, 3}); check(p, 2); p.v.clear(); check(p, 0);
""")
case("utils/grid_searcher", """
GridSearcher s; check(s, 0); s.incl(1, 2); check(s, 1);
s.incl(1, 2); check(s, 2); s.excl(1, 2); check(s, 1);
s.excl(1, 2); check(s, 0);
""")
case("utils/auto_rollback", """
TemporaryRollbackLog log; Int x = 1; check(log, 0);
Temporary(log, [&](auto& scope) {
    scope.set(x, 2); check(scope, 1); scope.set(x, 3); check(scope, 1);
}); check(log, 0); assert(x == 1);
""")

for header in ["segtree", "segtree_static_op", "segtree_var"]:
    case("collections/" + header, """
for (Int n : {0, 1, 4}) {
    auto s = initSegmentTree(std::vector<Int>(n, 2), std::plus<Int>{}, Int(0));
    check(s, n); if (n) { s.update(0, 7); check(s, n); assert(s.get(0, 1) == 7); }
}
""")
case("graph/private/combined_segtree", """
for (Int n : {0, 1, 4}) {
    auto s = initSegmentTree(std::vector<Int>(n, 2), std::plus<Int>{}, Int(0));
    check(s, n); if (n) { s.update(0, 7); check(s, n); }
}
""")
for header in ["lazysegtree", "lazysegtree_static_op"]:
    case("collections/" + header, """
auto map = [](Int f, Int s) { return f ? f : s; };
auto comp = [](Int f, Int g) { return f ? f : g; };
for (Int n : {0, 1, 4}) {
    auto s = initLazySegmentTree(std::vector<Int>(n), std::plus<Int>{}, Int(0), map, comp, Int(0));
    check(s, n); if (n) { s.update(0, 7); check(s, n); }
}
""")
for header in ["dualsegtree", "dualsegtree_static_op"]:
    case("collections/" + header, """
for (Int n : {0, 1, 4}) {
    auto s = initDualSegmentTree(std::vector<Int>(n), std::plus<Int>{}, std::plus<Int>{}, Int(0));
    check(s, n); s.apply(0, n, 3); check(s, n);
}
""")
case("collections/segtree_beats_template", """
for (Int n : {0, 1, 4}) {
    auto s = initRangeChminChmaxRangeSumMaxMin(std::vector<Int>(n));
    check(s, n); if (n) { s.update(0, 7); check(s, n); }
}
""")
case("collections/fenwick", """
for (Int n : {0, 1, 4}) {
    FenwickTree<Int> s(n); check(s, n);
    if (n) { s.add(0, 7); check(s, n); assert(s.get(0, 1) == 7); }
}
""")
case("collections/fenwick_avx2", """
for (Int n : {0, 1, 65}) {
    FenwickTreeAvx2 s(n); check(s, n);
    if (n) { s.add(0, 7); check(s, n); assert(s.get(0, 1) == 7); }
}
""")
case("collections/root_rangesum", """
for (Int n : {0, 1, 4}) {
    RootRangeSum<Int> s{std::vector<Int>(n)}; check(s, n);
    if (n) { s.update(0, 7); check(s, n); assert(s.get(0, 1) == 7); }
}
""")
case("collections/range_linear_add_range_min", """
for (Int n : {0, 1, 4}) {
    RangeLinearAddRangeMin s{std::vector<Int>(n)}; check(s, n);
    s.add(0, n, 1, 2); check(s, n);
}
""")
case("collections/dynamic_segtree", """
for (Int n : std::array<Int, 3>{0, 1, Int(1) << 50}) {
    auto s = initDynamicSegmentTree(n, std::plus<Int>{}, Int(0)); check(s, n);
    if (n) { s.update(n - 1, 7); check(s, n); assert(s.node_count() < 100); }
}
""")
case("collections/dynamic_lazysegtree", """
for (Int n : std::array<Int, 3>{0, 1, Int(1) << 50}) {
    auto s = initDynamicLazySegmentTree(n, std::plus<Int>{}, Int(0),
        std::plus<Int>{}, std::plus<Int>{}, Int(0), [](Int, Int) { return Int(0); });
    check(s, n); s.apply(0, n, 2); check(s, n);
}
""")
case("collections/persistent_segtree", """
for (Int n : {0, 1, 4}) {
    auto s = initPersistentSegmentTree(std::vector<Int>(n), std::plus<Int>{}, Int(0));
    check(s, n); if (n) { auto t = s.update(0, 7); check(t, n); check(s, n); }
}
""")
case("collections/persistent_lazysegtree", """
for (Int n : {0, 1, 4}) {
    auto s = initPersistentLazySegmentTree(std::vector<Int>(n), std::plus<Int>{}, Int(0),
        std::plus<Int>{}, std::plus<Int>{}, Int(0));
    check(s, n); if (n) { auto t = s.update(0, 7); check(t, n); check(s, n); }
}
""")
case("collections/persistent_array", """
for (Int n : {0, 1, 4}) {
    auto s = initPersistentArray(std::vector<Int>(n)); check(s, n);
    if (n) { auto t = s.change_value(0, 7); check(t, n); check(s, n); assert(s[0] == 0 && t[0] == 7); }
}
""")

for header in ["bitset", "bitset_avx2", "bitset_avx512", "staticbitset", "staticbitset_avx2", "staticbitset_avx512"]:
    fixed = header.startswith("static")
    case("collections/" + header, """
BitSet<65> s; check(s, 65); s[64] = true; check(s, 65);
BitSet<0> e; check(e, 0); BitSet<1> one; check(one, 1);
""" if fixed else """
for (Int n : {0, 1, 65}) {
    BitSet s(n); check(s, n); if (n) { s[n - 1] = true; check(s, n); }
}
""")
case("collections/combined", """
static_assert(BitSet<0>::size() == 0 && BitSet<1>::size() == 1 && BitSet<65>::size() == 65);
BitSet<65> s; check(s, 65); s[64] = true; check(s, 65); BitSet<0> e; check(e, 0);
""")
case("collections/bitvector", """
for (Int n : {0, 1, 63, 64, 65}) {
    BitVector s(n); check(s, n); if (n) s.set(n - 1); s.build(); check(s, n);
    assert(s.rank(n) == Int(n > 0));
}
""")
case("collections/binary_trie", """
BinaryTrie s(8); check(s, 0); s.incl(4); check(s, 1);
s.incl(4, 2); check(s, 3); s.excl(4, 3); check(s, 0);
""")
case("collections/bitset_binary_trie", """
BitSetBinaryTrie<BitSet> s(3); BitSet key(3); check(s, 0);
s.incl(key); check(s, 1); s.incl(key, 2); check(s, 3); s.excl(key, 3); check(s, 0);
""", "#include <cplib/collections/bitset.hpp>\n")
case("collections/intset", """
IntSet s(10); check(s, 0); s.incl(3); check(s, 1);
s.incl(3); check(s, 1); s.excl(3); check(s, 0);
""")
case("collections/hashset", """
HashSet<Int> s; check(s, 0); s.incl(3); check(s, 1);
s.incl(3); check(s, 1); s.excl(3); check(s, 0);
""")
case("collections/hashtable", """
HashTable<Int, Int> s; check(s, 0); s[3] = 7; check(s, 1);
s[3] = 8; check(s, 1); s.del(3); check(s, 0);
""")
case("collections/defaultdict", """
auto s = initDefaultDict<Int, Int>(0); check(s, 0);
s[3] = 7; check(s, 1); s[3] = 8; check(s, 1); s.del(3); check(s, 0);
""")
for header, typename in [("max_heapqueue", "MaxHeapQueue<Int>"), ("deletable_heapqueue", "Deletable_HeapQueue<Int>"), ("radix_heap", "RadixHeap<Int, Int>")]:
    case("collections/" + header, ("""
RadixHeap<Int, Int> s; check(s, 0); s.push(2, 3); check(s, 1);
s.push(2, 4); check(s, 2); s.pop(); check(s, 1); s.pop(); check(s, 0);
""" if header == "radix_heap" else f"""
{typename} s; check(s, 0); s.push(2); check(s, 1);
s.push(3); check(s, 2); s.pop(); check(s, 1); s.pop(); check(s, 0);
"""))
for header, make, push, pop in [("QSWAG", "initSWAG", "push", "pop"), ("SWAG", "initSWAG", "addLast", "popFirst")]:
    case("collections/" + header, f"""
auto s = {make}(std::plus<Int>{{}}, Int(0)); check(s, 0);
s.{push}(2); check(s, 1); s.{push}(3); check(s, 2);
assert(s.fold() == 5); s.{pop}(); check(s, 1); s.{pop}(); check(s, 0);
""")

for header, typename in [("unionfind", "UnionFind"), ("group_unionfind", "UnionFind"), ("weightedunionfind", "WeightedUnionFind<>"), ("rollback_unionfind", "RollbackUnionFind"), ("ppunionfind", "PartialPersistentUnionFind")]:
    op = "s.unite(0, 1, 3)" if header == "weightedunionfind" else "s.unite(0, 1)"
    more = "assert(s.size(0) == 2);" if header == "ppunionfind" else ""
    case("collections/" + header, f"""
{typename} empty(0), one(1); check(empty, 0); check(one, 1);
{typename} s(3); check(s, 3); {op}; check(s, 3); {more}
""")
for header in ["wordsizetree", "wordsizetree_avx2"]:
    typename = "WordsizeTree" if header == "wordsizetree" else "WordsizeTreeAvx2"
    case("collections/" + header, f"""
auto s = std::make_unique<{typename}>(); check(*s, 0);
s->incl(63); check(*s, 1); s->incl(63); check(*s, 1);
s->incl(64); check(*s, 2); s->excl(65); check(*s, 2);
s->excl(63); check(*s, 1); s->excl(63); check(*s, 1); s->excl(64); check(*s, 0);
auto built = std::make_unique<{typename}>(std::vector<bool>{{true, false, true}}); check(*built, 2);
auto copy = std::make_unique<{typename}>(*built); copy->excl(0); check(*copy, 1); check(*built, 2);
{typename} single(std::vector<bool>{{true}}); check(single, 1);
s->incl((Int(1) << 24) - 1); check(*s, 1); s->excl((Int(1) << 24) - 1); check(*s, 0);
std::set<Int> expected;
for (Int i = 0; i < 1000; ++i) {{
    Int key = (i * 137 + i / 7) % 251;
    if (i % 3) {{ s->incl(key); expected.insert(key); }}
    else {{ s->excl(key); expected.erase(key); }}
    check(*s, expected.size());
}}
""", "#include <memory>\n#include <set>\n")
case("collections/waveletmatrix", """
for (Int n : {0, 1, 4}) { WaveletMatrix s(std::vector<Int>(n, 2)); check(s, n); }
""")
case("collections/waveletmatrix_fenwick", """
for (Int n : {0, 1, 4}) {
    auto s = initWaveletMatrixFenwick(std::vector<std::pair<Int, Int>>(n, {2, 3}));
    check(s, n); if (n) { s.add(0, 4); check(s, n); }
}
""")
for header in ["graph/graph", "graph/private/combined_graph"]:
    case(header, """
for (Int n : {0, 1, 4}) {
    auto s = initUnWeightedDirectedGraph(n); check(s, n);
    if (n) { s.add_edge(0, 0); check(s, n); assert(s.len == n); }
}
""")
case("graph/range_edge_graph", """
auto s = initWeightedRangeGraph<Int>(3); check(s, s.graph().len);
auto before = s.size(); s.add_point_to_point_edge(0, 1, 2); check(s, before);
""")
case("tree/link_cut_tree", """
for (Int n : {0, 1, 4}) {
    std::vector<Int> values(n);
    detail::LinkCutTreeState<Int> state(values, std::plus<Int>{}, Int(0), {}); check(state, n);
    auto s = initLinkCutTree(values, std::plus<Int>{}, Int(0)); check(s, n);
    if (n) { s.update(0, 7); check(s, n); }
}
""")

case("collections/staticrangecount", """
for (Int n : {0, 1, 4}) { auto s = initStaticRangeCount(std::vector<Int>(n, 2)); check(s, n); assert(s.count(0, n, 2) == n); }
""")
case("collections/compressed_segtree", """
for (Int n : {0, 1, 4}) {
    std::vector<Int> keys(n); std::iota(keys.begin(), keys.end(), 0);
    auto s = initCompressedSegmentTree(keys, std::plus<Int>{}, Int(0)); check(s, n);
    if (n) { s.update(0, 7); check(s, n); }
}
auto duplicate = initCompressedSegmentTree(std::vector<Int>{1, 1, 2}, std::plus<Int>{}, Int(0)); check(duplicate, 2);
""")
case("collections/compressed_lazysegtree", """
for (Int n : {0, 1, 4}) {
    std::vector<Int> keys(n); std::iota(keys.begin(), keys.end(), 0);
    auto s = initCompressedLazySegmentTree(keys, std::plus<Int>{}, Int(0),
        std::plus<Int>{}, std::plus<Int>{}, Int(0)); check(s, n);
    if (n) { s.update(0, 7); check(s, n); }
    auto intervals = initCompressedLazySegmentTree(keys, std::plus<Int>{}, Int(0),
        std::plus<Int>{}, std::plus<Int>{}, Int(0), [](Int l, Int r) { return r - l; });
    check(intervals, std::max<Int>(0, n - 1));
}
""")
for header in ["compressed_fenwick2d", "compressed_segtree2d"]:
    make = "initCompressedFenwick2D(points)" if header == "compressed_fenwick2d" else "initCompressedSegmentTree2D(points, std::plus<Int>{}, Int(0))"
    update = "s.add(0, 0, 7)" if header == "compressed_fenwick2d" else "s.update(0, 0, 7)"
    case("collections/" + header, f"""
for (Int n : {{0, 1, 4}}) {{
    std::vector<std::tuple<Int, Int, Int>> points;
    for (Int i = 0; i < n; ++i) points.emplace_back(i, i, 2);
    auto s = {make}; check(s, n); if (n) {{ {update}; check(s, n); }}
}}
""")
case("collections/segtree_beats", """
struct S { Int sum = 0; bool fail = false; };
for (Int n : {0, 1, 4}) {
    auto op = [](S a, S b) { return S{a.sum + b.sum, false}; };
    auto map = [](Int, S s) { return s; };
    SegmentTreeBeats<S, Int> s(n, op, S{}, map, std::plus<Int>{}, Int(0));
    check(s, n); if (n) { s.update(0, S{7, false}); check(s, n); }
}
""")
case("collections/range_sort_array", """
for (Int n : {0, 1, 4}) {
    std::vector<Int> keys(n); std::iota(keys.begin(), keys.end(), 0);
    RangeSortArray s(keys, 10); check(s, n); s.sort(0, n, Descending); check(s, n);
    if (n) { s.update(0, 9); check(s, n); }
}
""")
case("collections/range_sort_segtree", """
for (Int n : {0, 1, 4}) {
    std::vector<Int> keys(n); std::iota(keys.begin(), keys.end(), 0);
    auto s = initRangeSortSegmentTree(keys, keys, 10, std::plus<Int>{}, Int(0));
    check(s, n); s.sort(0, n, Descending); check(s, n);
    if (n) { s.update(0, 9, 7); check(s, n); }
}
""")
for header in ["range_reverse_array", "range_reverse_array_monoid", "range_reverse_dualsegtree", "range_reverse_lazysegtree"]:
    make = {
        "range_reverse_array": "initRangeReverseArray(v)",
        "range_reverse_array_monoid": "initRangeReverseArrayMonoid(v, std::plus<Int>{}, Int(0))",
        "range_reverse_dualsegtree": "initRangeReverseDualSegmentTree(v, std::plus<Int>{}, std::plus<Int>{}, Int(0))",
        "range_reverse_lazysegtree": "initRangeReverseLazySegmentTree(v, std::plus<Int>{}, Int(0), std::plus<Int>{}, std::plus<Int>{}, Int(0))",
    }[header]
    case("collections/" + header, f"""
std::vector<Int> v; auto s = {make}; check(s, 0);
s.insert(0, 7); check(s, 1); s.insert(1, 8); check(s, 2);
s.reverse(0, 2); check(s, 2); s.erase(0); check(s, 1); s.erase(0); check(s, 0);
""" if header != "range_reverse_array" else """
for (Int n : {0, 1, 4}) {
    auto s = initRangeReverseArray(std::vector<Int>(n)); check(s, n);
    s.reverse(0, n); check(s, n); if (n) { s.update(0, 7); check(s, n); }
}
""")
for header in ["retroactive_priority_queue", "compressed_retroactive_priority_queue", "dynamic_retroactive_priority_queue"]:
    make = {
        "retroactive_priority_queue": "RetroactivePriorityQueue<Int>(3)",
        "compressed_retroactive_priority_queue": "CompressedRetroactivePriorityQueue<Int, Int>(std::vector<Int>{0, 1, 2})",
        "dynamic_retroactive_priority_queue": "DynamicRetroactivePriorityQueue<Int, Int>()",
    }[header]
    case("collections/" + header, f"""
auto s = {make}; check(s, 0); s.setPush(0, 7); check(s, 1);
s.setPush(1, 8); check(s, 2); s.setPop(2); check(s, 1);
s.erase(2); check(s, 2); s.erase(1); check(s, 1); s.erase(0); check(s, 0);
""")
case("collections/indexed_retroactive_priority_queue", """
IndexedRetroactivePriorityQueue<Int> s; check(s, 0);
s.insertPush(0, 7); check(s, 1); s.insertPop(1); check(s, 0);
s.erase(1); check(s, 1); s.erase(0); check(s, 0);
""")
case("convolution/relaxed_convolution", """
using M = StaticMontgomeryModint<998244353>;
RelaxedConvolution<M> s(4); check(s, 0); s.append(1, 1); check(s, 1); s.append(0, 0); check(s, 2);
RelaxedInv<M> inv(4); RelaxedExp<M> exp(4); RelaxedLog<M> log(4); RelaxedSqrt<M> sqrt(4); RelaxedPow<M> pow(4, 2);
check(inv, 0); check(exp, 0); check(log, 0); check(sqrt, 0); check(pow, 0);
for (Int i = 0; i < 4; ++i) {
    inv.append(i == 0); exp.append(0); log.append(i == 0); sqrt.append(i == 0); pow.append(i == 0);
    check(inv, i + 1); check(exp, i + 1); check(log, i + 1); check(sqrt, i + 1); check(pow, i + 1);
}
""")
case("convolution/semi_relaxed_convolution", """
using M = StaticMontgomeryModint<998244353>;
auto s = initSemiRelaxedConvolution(std::vector<M>{1, 2}); check(s, 0);
s.append(1); check(s, 1); s.append(2); check(s, 2);
""")
case("str/static_string", """
auto s = toStaticString("abc"); assert(s.size() == 3 && len(s) == 3 && sz(s) == 3);
assert((s.len)() == 3 && (cplib::len)(s) == 3);
""", "#define LOCAL_TEST\n#include <cplib/str/static_string.hpp>\n#include <cplib/tmpl/sheep.hpp>\n")

PREAMBLE = """
#include <functional>
#include <iterator>
#include <ranges>
#include <type_traits>
{extra}
#include <cplib/{header}.hpp>
using namespace cplib;
template<class T> void check(const T& value, Int expected) {{
    static_assert(std::is_same_v<decltype(value.size()), Int>);
    assert(value.size() == expected && std::size(value) == expected);
    assert(std::ranges::size(value) == expected);
    if constexpr (requires {{ (value.len)(); }}) assert((value.len)() == expected);
}}
int main() {{ {body} }}
"""


def run_case(item):
    header, body, extra = item
    with tempfile.TemporaryDirectory(prefix="cplib_size_") as directory:
        path = Path(directory)
        source, binary = path / "test.cpp", path / "test"
        source.write_text(PREAMBLE.format(header=header, body=body, extra=extra))
        command = shlex.split(os.environ.get("CXX", "g++"))
        command += ["-std=c++20", "-O2", "-I", str(ROOT / "src"), str(source), "-o", str(binary)]
        result = subprocess.run(command, capture_output=True, text=True, timeout=90)
        if result.returncode:
            return header, result.stderr
        result = subprocess.run([str(binary)], capture_output=True, text=True, timeout=30)
        return header, result.stderr if result.returncode else ""


if __name__ == "__main__":
    with ThreadPoolExecutor(max_workers=3) as pool:
        failures = [(header, error) for header, error in pool.map(run_case, CASES) if error]
    for header, error in failures:
        print(header + ":\n" + error)
    print(f"{len(CASES) - len(failures)}/{len(CASES)} size API configurations passed")
    raise SystemExit(bool(failures))
