#include <cplib/collections/merge_sort_tree.hpp>
#include <array>
#include <bit>
#include <cassert>
#include <functional>
#include <span>
#include <string>
#include <utility>
#include <vector>

using namespace cplib;

int main() {
    for (Int n = 0; n <= 65; ++n) {
        std::string expected;
        std::vector<std::string> leaves;
        for (Int i = 0; i < n; ++i) {
            expected.push_back(static_cast<char>('a' + i % 26));
            leaves.emplace_back(1, expected.back());
        }
        Int merges = 0;
        const auto tree = initMergeSortTree(leaves, [&](const auto &left, const auto &right) {
            ++merges;
            return left + right;
        });
        assert(tree.size() == n);
        assert(merges == std::max<Int>(0, n - 1));
        for (Int l = 0; l <= n; ++l) {
            for (Int r = l; r <= n; ++r) {
                Int visits = 0;
                std::string result;
                forEachNode(tree, l, r, [&](const std::string &value) {
                    result += value;
                    ++visits;
                });
                assert(result == expected.substr(l, r - l));
                assert(visits <= 2 * static_cast<Int>(std::bit_width(static_cast<UInt>(n))));
                if (l == r) {
                    assert(visits == 0);
                } else if (l == 0 && r == n) {
                    assert(visits == 1);
                }
            }
        }
    }

    std::array<std::string, 5> leaves{"a", "b", "c", "d", "e"};
    auto tree = initMergeSortTree<std::string>(std::span(leaves), std::plus<>{});
    std::vector<std::string> parts;
    tree.forEachNode(1, 5, [&](std::string &value) {
        parts.push_back(value);
        value += "!";
    });
    assert((parts == std::vector<std::string>{"b", "cde"}));
    parts.clear();
    forEachNode(tree, 1, 5, [&](const std::string &value) { parts.push_back(value); });
    assert((parts == std::vector<std::string>{"b!", "cde!"}));
    std::as_const(tree).forEachNode(0, 5, [](const auto &value) { assert(value == "abcde"); });
    assert(leaves[1] == "b");

    MergeSortTree flags(std::vector<bool>{false, true}, std::logical_or<>{});
    flags.forEachNode(0, 2, [](bool &value) {
        assert(value);
        value = false;
    });
    flags.forEachNode(0, 2, [](const bool &value) { assert(!value); });
}
