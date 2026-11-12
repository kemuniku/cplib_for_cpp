#include <cplib/collections/avl_seq.hpp>
#include <cplib/collections/range_reverse_array.hpp>
#include <cplib/collections/avl_monoid.hpp>
#include <cplib/collections/avl_dualsegtree.hpp>
#include <cplib/collections/avl_lazysegtree.hpp>
#include <cplib/collections/double_ended_priority_queue.hpp>
#include <cplib/collections/lazy_leftist_heap.hpp>
#include <cplib/collections/static_range_distinct.hpp>
#include <cplib/collections/defaultdict.hpp>
#include <cplib/collections/SWAG.hpp>
#include <cplib/str/suffix_automaton_array.hpp>
#include <cplib/str/suffix_automaton_table.hpp>
#include <cplib/str/palindromic_tree.hpp>
#include <cplib/str/palindromic_tree_table.hpp>
#include <cplib/str/trie_string.hpp>
#include <cplib/str/hamming_distance.hpp>
#include <random>
#include <set>
#include <map>

using namespace cplib;

struct Sum {
    Int value = 0, length = 0;
    bool operator==(const Sum &) const = default;
};
void testAvl(std::mt19937_64 &rng) {
    auto op = [](Sum a, Sum b) { return Sum{a.value + b.value, a.length + b.length}; };
    auto mapping = [](Int delta, Sum x) { return Sum{x.value + delta * x.length, x.length}; };
    auto composition = [](Int newer, Int older) { return newer + older; };
    auto tree = initAvlLazySegmentTree(std::vector<Sum>{}, op, Sum{}, mapping, composition, Int(0));
    std::vector<Sum> values;
    for (Int step = 0; step < 8000; ++step) {
        Int kind = rng() % 6;
        if (kind == 0 || values.empty()) {
            Int index = rng() % (values.size() + 1), value = rng() % 20;
            tree.insert(index, {value, 1});
            values.insert(values.begin() + index, {value, 1});
        } else if (kind == 1) {
            Int l = rng() % values.size(), r = l + rng() % (values.size() - l + 1);
            tree.erase(l, r);
            values.erase(values.begin() + l, values.begin() + r);
        } else if (kind == 2) {
            Int l = rng() % (values.size() + 1), r = rng() % (values.size() + 1);
            if (l > r)
                std::swap(l, r);
            Int delta = rng() % 5;
            tree.apply(l, r, delta);
            for (Int i = l; i < r; ++i)
                values[i].value += delta;
        } else if (kind == 3) {
            Int index = rng() % values.size(), value = rng() % 20;
            tree.update(index, {value, 1});
            values[index] = {value, 1};
        } else if (kind == 4) {
            Int k = rng() % (values.size() + 1);
            auto [left, right] = tree.split(k);
            assert(tree.size() == 0 && left.size() == k && right.size() == Int(values.size()) - k);
            tree = merge(left, right);
            assert(left.size() == 0 && right.size() == 0);
        } else {
            auto copy = tree;
            tree.apply(0, tree.size(), 1);
            assert(copy.toSeq() == values);
            for (auto &value : values)
                ++value.value;
        }
        assert(tree.size() == Int(values.size()) && tree.toSeq() == values);
        Int l = rng() % (values.size() + 1), r = rng() % (values.size() + 1);
        if (l > r)
            std::swap(l, r);
        Sum expected;
        for (Int i = l; i < r; ++i)
            expected = op(expected, values[i]);
        assert(tree.fold(l, r) == expected);
        Int limit = rng() % 150, right = l, sum = 0;
        while (right < Int(values.size()) && sum + values[right].value <= limit)
            sum += values[right++].value;
        assert(tree.max_right(l, [&](Sum x) { return x.value <= limit; }) == right);
        Int left = r;
        sum = 0;
        while (left > 0 && sum + values[left - 1].value <= limit)
            sum += values[--left].value;
        assert(tree.min_left(r, [&](Sum x) { return x.value <= limit; }) == left);
    }
    std::vector<Sum> iterated;
    for (const auto &value : tree)
        iterated.push_back(value);
    assert(iterated == values);
    struct Affine { Int a, b; };
    auto affine = initAvlLazySegmentTree(std::vector<Sum>(20, Sum{1, 1}), op, Sum{},
        [](Affine f, Sum x) { return Sum{f.a * x.value + f.b * x.length, x.length}; },
        [](Affine f, Affine g) { return Affine{f.a * g.a, f.a * g.b + f.b}; }, Affine{1, 0});
    std::vector<Int> reference(20, 1);
    for (Int q = 0; q < 1000; ++q) {
        Int l = rng() % 21, r = rng() % 21;
        if (l > r)
            std::swap(l, r);
        Affine action{Int(rng() % 2), Int(rng() % 10)};
        affine.apply(l, r, action);
        for (Int i = l; i < r; ++i)
            reference[i] = action.a * reference[i] + action.b;
        assert(affine.fold().value == std::accumulate(reference.begin(), reference.end(), Int(0)));
    }
    auto reversible = initRangeReverseArray(std::vector<Int>{1, 2, 4});
    assert(reversible.insort(2) == 2);
    reversible.reverse(0, reversible.size());
    auto [leftHalf, rightHalf] = reversible.split(2);
    assert(reversible.size() == 0);
    reversible = merge(leftHalf, rightHalf);
    assert((reversible.toSeq() == std::vector<Int>{4, 2, 2, 1}));
    auto seq = toAvlSeq(std::vector<Int>{1, 2, 4});
    assert(seq.insort(2) == 2);
    assert(seq[from_end(1)] == 4);
    auto strings = initAvlMonoid(std::vector<std::string>{"a", "b", "c"},
        [](auto a, auto b) { return a + b; }, std::string{});
    assert(strings.fold(1, 3) == "bc");
    strings.update(1, "xy");
    assert(strings.fold() == "axyc");
    auto dual = initAvlDualSegmentTree(std::vector<Int>{1, 2, 3},
        [](Int f, Int x) { return f + x; }, composition, Int(0));
    dual.apply(0, 3, 10);
    dual.insert(1, 0);
    assert((dual.toSeq() == std::vector<Int>{11, 0, 12, 13}));
}
void testQueuesAndDistinct(std::mt19937_64 &rng) {
    DoubleEndedPriorityQueue<Int> queue;
    std::multiset<Int> expected;
    for (Int step = 0; step < 20000; ++step) {
        Int kind = rng() % 3;
        if (kind == 0 || expected.empty()) {
            Int value = Int(rng() % 1001) - 500;
            queue.push(value);
            expected.insert(value);
        } else if (kind == 1) {
            assert(queue.popMin() == *expected.begin());
            expected.erase(expected.begin());
        } else {
            auto it = std::prev(expected.end());
            assert(queue.popMax() == *it);
            expected.erase(it);
        }
        assert(queue.size() == Int(expected.size()));
        if (!expected.empty())
            assert(queue.min() == *expected.begin() && queue.max() == *expected.rbegin());
    }
    std::vector<Int> input(1000);
    for (auto &x : input)
        x = Int(rng() % 21) - 10;
    auto built = toDoubleEndedPriorityQueue(input);
    std::sort(input.begin(), input.end());
    for (Int x : input)
        assert(built.popMin() == x);
    LazyLeftistHeapPool<Int, Int> pool;
    Int a = -1, b = -1;
    for (Int i = 0; i < 100; ++i)
        a = pool.meld(a, pool.singleton(i, i));
    pool.addAll(a, 20);
    for (Int i = 0; i < 100; ++i)
        b = pool.meld(b, pool.singleton(i, i + 100));
    Int root = pool.meld(a, b), previous = -1;
    while (root >= 0) {
        assert(pool.top(root).first >= previous);
        previous = pool.top(root).first;
        root = pool.pop(root);
    }
    auto distinct = initStaticRangeDistinct(input);
    auto comparison = initStaticRangeDistinct(input, [](Int a, Int b) { return a < b ? -1 : a > b ? 1 : 0; });
    for (Int q = 0; q < 3000; ++q) {
        Int l = rng() % 1001, r = rng() % 1001;
        if (l > r)
            std::swap(l, r);
        Int count = std::set<Int>(input.begin() + l, input.begin() + r).size();
        assert(distinct.count(l, r) == count && comparison.count(l, r) == count);
    }
    DefaultDict<Int, Int> dict(7);
    assert(dict.get(3) == 7 && dict.size() == 0);
    auto swag = newSwagWith([](Int a, Int b) { return a + b; }, Int(0));
    swag.addLast(5);
    assert(swag.fold() == 5);
}
void testStrings(std::mt19937_64 &rng) {
    for (Int trial = 0; trial < 100; ++trial) {
        auto fixed = initSuffixAutomatonArray<'a', 'c'>();
        auto table = initSuffixAutomatonTable<char>();
        PalindromicTree palindrome(3);
        TablePalindromicTree<char> palindromeTable;
        std::string s;
        for (Int step = 0; step < 30; ++step) {
            char c = 'a' + rng() % 3;
            s += c;
            fixed.extend(c);
            table.extend(c);
            auto node = palindrome.add(c - 'a');
            auto tableNode = palindromeTable.add(c);
            assert(node->size() == tableNode->size());
            std::set<std::string> substrings;
            std::map<std::vector<Int>, Int> palindromes;
            Int longest = 0;
            for (Int l = 0; l < Int(s.size()); ++l) {
                for (Int r = l + 1; r <= Int(s.size()); ++r) {
                    auto sub = s.substr(l, r - l), reverse = sub;
                    substrings.insert(sub);
                    std::reverse(reverse.begin(), reverse.end());
                    if (sub == reverse) {
                        std::vector<Int> key;
                        for (char x : sub)
                            key.push_back(x - 'a');
                        ++palindromes[key];
                        if (r == Int(s.size()))
                            longest = std::max(longest, r - l);
                    }
                }
            }
            assert(fixed.countDistinctSubstrings() == Int(substrings.size()));
            assert(table.countDistinctSubstrings() == Int(substrings.size()));
            assert(fixed.contains(s) && table.contains(s));
            assert(fixed.contains("") && table.contains(""));
            assert(!fixed.contains("d") && !table.contains("d"));
            assert(palindrome.size() == Int(s.size()) && node->size() == longest);
            palindrome.update_count();
            palindrome.update_count();
            palindromeTable.update_count();
            palindromeTable.update_count();
            assert(palindrome.nodes.size() == palindromes.size() + 2);
            for (std::size_t i = 2; i < palindrome.nodes.size(); ++i) {
                const auto &p = palindrome.nodes[i];
                assert(p->count() == palindromes.at(p->get_palindrome()));
                assert(p->count() == palindromeTable.nodes[i]->count());
                auto suffix = p->suffix_link;
                auto series = suffix;
                while (series->number >= 2 && p->diff() == series->diff())
                    series = series->suffix_link;
                assert(series == p->series_link);
            }
        }
        auto copy = palindrome;
        palindrome.add(0);
        assert(copy.size() == Int(s.size()));
    }
    assert(SuffixAutomatonTable<char>("ab").size() == 2);
    assert(SuffixAutomatonArray<>("ab").size() == 2);
    assert(TablePalindromicTree<char>("ab").size() == 2);
    auto alphabet = initPalindromicTree<'a', 'z'>("ababa");
    assert(alphabet.size() == 5 && alphabet.last_node->character == 'a');
    auto palindromeValues = alphabet.get_palindrome(alphabet.last_node);
    assert(std::string(palindromeValues.begin(), palindromeValues.end()) == "ababa");
    auto reservedAlphabet = initPalindromicTree<'a', 'z'>(10);
    reservedAlphabet.add('a');
    enum class Symbol { A = -2, B = -1, C = 0 };
    auto symbols = initPalindromicTree<Symbol::A, Symbol::C>();
    symbols.add(Symbol::A);
    symbols.add(Symbol::B);
    symbols.add(Symbol::A);
    assert((symbols.get_palindrome(symbols.last_node) == std::vector<Symbol>{Symbol::A, Symbol::B, Symbol::A}));
    PalindromicTree largeAlphabet(1000000);
    auto firstPalindrome = largeAlphabet.add(123456);
    assert(largeAlphabet.nodes[0]->link.size() == 1);
    assert(largeAlphabet.nodes[0]->link.hasKey(123456));
    assert(largeAlphabet.nodes[0]->link[123456] == firstPalindrome);
    assert(largeAlphabet.nodes[0]->link.getOrDefault(-1) == nullptr);
    auto copied = largeAlphabet;
    assert(copied.nodes[0]->link[123456] != firstPalindrome);
    auto arena = initTrieStringArray<'a', 'c'>();
    auto root = arena.root();
    std::vector<decltype(root)> nodes{root};
    std::vector<std::string> values{std::string{}};
    for (Int i = 0; i < 10000; ++i) {
        Int parent = rng() % nodes.size();
        char c = 'a' + rng() % 3;
        nodes.push_back(nodes[parent].add(c));
        values.push_back(values[parent] + c);
        assert(nodes.back().to_string() == values.back());
        Int other = rng() % nodes.size(), common = 0;
        while (common < Int(std::min(values.back().size(), values[other].size())) && values.back()[common] == values[other][common])
            ++common;
        assert(nodes.back().lcp(nodes[other]) == common);
        assert((nodes.back() < nodes[other]) == (values.back() < values[other]));
        assert(nodes.back().prefix(common) == nodes[other].prefix(common));
        assert(nodes.back().pop() == nodes[parent]);
    }
    auto surviving = []() {
        auto local = initTrieString();
        return local.root().add('x').add('y');
    }();
    assert(surviving.to_string() == "xy");
    auto integers = initTrieStringTable<Int>();
    auto a = integers.root().add(-3).add(1000000000000LL);
    assert((a.toSeq() == std::vector<Int>{-3, 1000000000000LL}));
    for (Int n = 0; n < 300; ++n) {
        std::string a(n, '\0'), b(n, '\0');
        Int count = 0;
        for (Int i = 0; i < n; ++i) {
            a[i] = rng() % 256;
            b[i] = rng() % 256;
            count += a[i] != b[i];
        }
        assert(hammingDistance(a, b) == count);
    }
}
int main() {
    std::mt19937_64 rng(13579);
    testAvl(rng);
    testQueuesAndDistinct(rng);
    testStrings(rng);
}
