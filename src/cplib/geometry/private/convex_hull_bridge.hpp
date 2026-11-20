#pragma once
#include <cplib/geometry/base.hpp>

namespace cplib::detail {

// y で分割した左右の左側凸包の共通接線を探す。各節点は両子をつなぐ接線を持つ。
// 毎回どちらかの子へ降りるので、両木の高さの和に比例する。
// 同一直線上の点も残る側を選ぶ。座標幅 D に対し中間値は最大 4D^3。
template <class Node, class Leaf, class Child, class Ends>
std::pair<Node, Node> convexHullBridge(Node left, Node right, __int128 split,
                                       const Leaf &leaf, const Child &child, const Ends &ends) {
    auto turn = [](Point<__int128> a, Point<__int128> b, Point<__int128> c) {
        return cross(b - a, c - a);
    };
    while (!leaf(left) || !leaf(right)) {
        const auto [a, b] = ends(left);
        const auto [c, d] = ends(right);
        if (!leaf(left) && turn(a, b, c) > 0) {
            left = child(left, false);
        } else if (!leaf(right) && turn(b, c, d) > 0) {
            right = child(right, true);
        } else if (leaf(left)) {
            right = child(right, false);
        } else if (leaf(right)) {
            left = child(left, true);
        } else {
            const __int128 s = turn(a, b, c), t = turn(b, a, d);
            // 両子の接線の交点と分割線との位置関係を、除算せず比較する。
            if (s + t == 0 || s * (d.y - split) + t * (c.y - split) < 0)
                left = child(left, true);
            else
                right = child(right, false);
        }
    }
    return {left, right};
}

}
