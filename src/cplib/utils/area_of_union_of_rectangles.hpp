#pragma once
#include <cplib/common.hpp>
#include <tuple>

namespace cplib {
namespace detail {
template <class T> void rectangleRadixSort(std::vector<T> &a) {
    if (a.size() < 2)
        return;
    constexpr UInt signBit = UInt(1) << 63, mask = 2047;
    UInt first = UInt(a[0].coordinate), varying = 0;
    for (const auto &e : a)
        varying |= UInt(e.coordinate) ^ first;
    std::vector<T> scratch(a.size());
    for (int shift = 0; shift < 64 && (varying >> shift) != 0; shift += 11) {
        if (((varying >> shift) & mask) == 0)
            continue;
        std::array<Int, 2048> offsets{};
        for (const auto &e : a)
            ++offsets[((UInt(e.coordinate) ^ signBit) >> shift) & mask];
        Int total = 0;
        for (auto &x : offsets) {
            Int count = x;
            x = total;
            total += count;
        }
        for (const auto &e : a)
            scratch[offsets[((UInt(e.coordinate) ^ signBit) >> shift) & mask]++] = e;
        a.swap(scratch);
    }
}

template <class Length, class Index>
Int rectangleSweep(std::span<const std::tuple<Int, Int, Int, Int>> rectangles,
                   const std::vector<Int> &positions, const std::vector<Int> &seps) {
    struct Event {
        Int coordinate;
        Index left, right;
    };

    struct Node {
        Index cover{};
        Length length{}, width{};
    };

    std::vector<Event> events;
    events.reserve(2 * rectangles.size());
    for (Int i = 0; i < Int(rectangles.size()); ++i) {
        auto [l, d, r, u] = rectangles[i];
        if (l == r || d == u)
            continue;
        Index x = positions[2 * i], z = positions[2 * i + 1];
        events.push_back({d, x, z});
        events.push_back({u, Index(~x), z});
    }
    rectangleRadixSort(events);
    Int n = seps.size() - 1, size = 1;
    while (size < n)
        size *= 2;
    std::vector<Node> nodes(2 * size);
    for (Int i = 0; i < n; ++i)
        nodes[size + i].width = Length(seps[i + 1] - seps[i]);
    for (Int i = size - 1; i >= 1; --i)
        nodes[i].width = nodes[2 * i].width + nodes[2 * i + 1].width;
    auto change = [&](Int node, Int delta) {
        Length difference{};
        Index oldCover = nodes[node].cover;
        nodes[node].cover += Index(delta);
        if (nodes[node].cover == 0 || oldCover == 0) {
            Length old = nodes[node].length;
            if (nodes[node].cover > 0)
                nodes[node].length = nodes[node].width;
            else if (node < size)
                nodes[node].length = nodes[node * 2].length + nodes[node * 2 + 1].length;
            else
                nodes[node].length = 0;
            difference = nodes[node].length - old;
        }
        return difference;
    };
    Int previousY = events[0].coordinate, result = 0;
    for (auto event : events) {
        result += Int(nodes[1].length) * (event.coordinate - previousY);
        previousY = event.coordinate;
        Int delta = event.left < 0 ? -1 : 1,
            l = (event.left < 0 ? ~Int(event.left) : Int(event.left)) + size,
            r = Int(event.right) + size, leftNode = l, rightNode = r - 1;
        Length leftDiff{}, rightDiff{};
        // 子の被覆長の差分だけを伝播し、変化が消えた時点で打ち切る。
        while (leftNode > 0) {
            if (leftNode == rightNode) {
                leftDiff += rightDiff;
                rightDiff = 0;
            }
            if (leftDiff != 0) {
                if (nodes[leftNode].cover == 0)
                    nodes[leftNode].length += leftDiff;
                else
                    leftDiff = 0;
            }
            if (rightDiff != 0) {
                if (nodes[rightNode].cover == 0)
                    nodes[rightNode].length += rightDiff;
                else
                    rightDiff = 0;
            }
            if (l < r) {
                if (l & 1)
                    leftDiff += change(l++, delta);
                if (r & 1)
                    rightDiff += change(--r, delta);
            }
            if (l == r && leftDiff == 0 && rightDiff == 0)
                break;
            l >>= 1;
            r >>= 1;
            leftNode >>= 1;
            rightNode >>= 1;
        }
    }
    return result;
}
}

// (左,下,右,上)の長方形群。走査線O(N log N)、基数ソートO(N)、追加領域O(N)。
// 軸に平行な長方形の和集合の面積を返す。時間O(N log N)、追加空間O(N)。
// 各要素は(l,d,r,u)=(左端,下端,右端,上端)。l <= r、d <= uを満たすこと。
// 空入力・面積0の長方形・負の座標にも対応する。
// 全体のx座標幅、y座標の差、面積はIntに収まること。
inline Int area_of_union_of_rectangles(std::span<const std::tuple<Int, Int, Int, Int>> rectangles) {
    struct Endpoint {
        Int coordinate, index;
    };

    std::vector<Endpoint> endpoints;
    endpoints.reserve(2 * rectangles.size());
    for (Int i = 0; i < Int(rectangles.size()); ++i) {
        auto [l, d, r, u] = rectangles[i];
        if (l == r || d == u)
            continue;
        endpoints.push_back({l, 2 * i});
        endpoints.push_back({r, 2 * i + 1});
    }
    if (endpoints.empty())
        return 0;
    detail::rectangleRadixSort(endpoints);
    std::vector<Int> positions(2 * rectangles.size()), seps;
    seps.reserve(endpoints.size());
    for (auto e : endpoints) {
        if (seps.empty() || seps.back() != e.coordinate)
            seps.push_back(e.coordinate);
        positions[e.index] = seps.size() - 1;
    }
    std::vector<Endpoint>().swap(endpoints);
    // 通常は12byteノード。座標幅とノード数に応じて64bit版を選ぶ。
    if (rectangles.size() <= std::size_t(std::numeric_limits<std::int32_t>::max() / 2)) {
        if (seps.back() - seps.front() <= std::numeric_limits<std::int32_t>::max())
            return detail::rectangleSweep<std::int32_t, std::int32_t>(rectangles, positions, seps);
        return detail::rectangleSweep<Int, std::int32_t>(rectangles, positions, seps);
    }
    return detail::rectangleSweep<Int, Int>(rectangles, positions, seps);
}
}
