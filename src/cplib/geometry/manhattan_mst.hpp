#pragma once
#include <cplib/geometry/base.hpp>
#include <cplib/collections/unionfind.hpp>
#include <map>
#include <tuple>

namespace cplib {
// マンハッタン距離の最小全域木を元の頂点番号の組で返す。時間 O(N log N)、領域 O(N)。
// 座標の絶対値はIntの最大値/4以下。重複点は別頂点として扱う。
template <class T, class Id = void> requires (std::is_integral_v<T> && std::is_signed_v<T>)
// 整数幅: 内部でIntへ昇格し、B=4M（符号反転・座標和・辺の重みを含む）。
// M=10^9なら内部は33ビット必要だが、入力はPoint<int>/pair<int,int>で可。総重みは計算しない。
std::vector<std::pair<Int, Int>> manhattan_mst(const std::vector<Point<T, Id>> &points) {
    constexpr Int bound = std::numeric_limits<Int>::max() / 4;
    std::vector<Point<Int>> sites;
    for (const auto &p : points) {
        assert(-bound <= p.x && p.x <= bound && -bound <= p.y && p.y <= bound);
        sites.push_back({Int(p.x), Int(p.y)});
    }
    std::vector<Int> order(sites.size());
    std::iota(order.begin(), order.end(), 0);
    std::vector<std::tuple<Int, Int, Int>> candidates;
    for (Int direction = 0; direction < 4; ++direction) {
        std::sort(order.begin(), order.end(), [&](Int a, Int b) {
            return sites[a].x + sites[a].y < sites[b].x + sites[b].y;
        });
        std::map<Int, Int> sweep;
        for (Int i : order) {
            auto it = sweep.lower_bound(-sites[i].y);
            while (it != sweep.end()) {
                Int j = it->second;
                if (sites[i].x - sites[j].x < sites[i].y - sites[j].y)
                    break;
                Int cost = sites[i].x - sites[j].x + sites[i].y - sites[j].y;
                candidates.emplace_back(cost, i, j);
                it = sweep.erase(it);
            }
            sweep[-sites[i].y] = i;
        }
        for (auto &p : sites) {
            if (direction & 1)
                p.x = -p.x;
            else
                std::swap(p.x, p.y);
        }
    }
    std::sort(candidates.begin(), candidates.end());
    UnionFind uf(points.size());
    std::vector<std::pair<Int, Int>> result;
    for (auto [cost, u, v] : candidates) {
        if (!uf.issame(u, v)) {
            uf.unite(u, v);
            result.emplace_back(std::min(u, v), std::max(u, v));
        }
    }
    return result;
}
// 座標の組からマンハッタン最小全域木を返す。時間 O(N log N)、領域 O(N)。
template <class T> requires (std::is_integral_v<T> && std::is_signed_v<T>)
// 整数幅: 内部でIntへ昇格し、B=4M（符号反転・座標和・辺の重みを含む）。
// M=10^9なら内部は33ビット必要だが、入力はPoint<int>/pair<int,int>で可。総重みは計算しない。
auto manhattan_mst(const std::vector<std::pair<T, T>> &points) {
    std::vector<Point<T>> converted;
    for (auto [x, y] : points)
        converted.push_back({x, y});
    return manhattan_mst(converted);
}
}
