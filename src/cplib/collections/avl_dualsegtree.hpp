#pragma once
#include <cplib/collections/private/avl_sequence.hpp>

namespace cplib {
template <class S, class F> using AvlDualSegmentTree = detail::AvlSequence<S, F, false, true>;
// 区間作用付きAVL列を配列から作る。合成は新しい作用を左に置く。時間・領域 O(N)。
template <class S, class F, class Map, class Comp>
auto initAvlDualSegmentTree(const std::vector<S> &values, Map mapping, Comp composition, F id) {
    return AvlDualSegmentTree<S, F>(values, {}, S{}, mapping, composition, id);
}
// 同じ初期値をn個持つ区間作用付きAVL列を作る。時間・領域 O(n)。
template <class S, class F, class Map, class Comp>
auto initAvlDualSegmentTree(Int n, S initialValue, Map mapping, Comp composition, F id) {
    assert(n >= 0);
    return initAvlDualSegmentTree(std::vector<S>(n, initialValue), mapping, composition, id);
}
// 配列から区間作用付きAVL列を作る。時間・領域 O(N)。
template <class S, class Map, class Comp, class F>
auto newAvlDualSegWith(const std::vector<S> &values, Map mapping, Comp composition, F id) {
    return initAvlDualSegmentTree(values, mapping, composition, id);
}
// 同じ初期値をn個持つ区間作用付きAVL列を作る。時間・領域 O(n)。
template <class S, class Map, class Comp, class F>
auto newAvlDualSegWith(Int n, S initialValue, Map mapping, Comp composition, F id) {
    return initAvlDualSegmentTree(n, initialValue, mapping, composition, id);
}
}
