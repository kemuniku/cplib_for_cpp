#pragma once
#include <cplib/geometry/base.hpp>
#include <optional>
#include <stdexcept>
#include <tuple>

namespace cplib {
template <class T> struct ClosestPairResult {
    T distanceSquared;
    std::pair<Int, Int> indices;
};
// 最近点対の二乗距離と元の添字(i<j)を返す。2点未満はnullopt。時間 O(N log N)、領域 O(N)。
// 整数は各軸の幅2*10^9以下、浮動小数点は有限座標・有限二乗距離と非零差の二乗が非零であること。
// 整数はInt、浮動小数点はdoubleで計算し、不正な座標はinvalid_argument。
template <class T> requires (std::is_signed_v<T>)
auto closest_pair(const std::vector<Point<T>> &points)
    -> std::optional<ClosestPairResult<std::conditional_t<std::is_integral_v<T>, Int, double>>> {
    using R = std::conditional_t<std::is_integral_v<T>, Int, double>;
    struct Site {
        R x, y;
        Int id;
    };
    std::vector<Site> sites;
    for (Int i = 0; i < Int(points.size()); ++i) {
        Site p{R(points[i].x), R(points[i].y), i};
        if constexpr (!std::is_integral_v<T>) {
            if (!std::isfinite(p.x) || !std::isfinite(p.y))
                throw std::invalid_argument("座標は有限である必要があります");
        }
        sites.push_back(p);
    }
    if (sites.empty())
        return std::nullopt;
    auto xOrder = [](Site a, Site b) { return std::tie(a.x, a.y, a.id) < std::tie(b.x, b.y, b.id); };
    auto yOrder = [](Site a, Site b) { return std::tie(a.y, a.x, a.id) < std::tie(b.y, b.x, b.id); };
    std::sort(sites.begin(), sites.end(), xOrder);
    auto scratch = sites;
    std::sort(scratch.begin(), scratch.end(), yOrder);
    if constexpr (std::is_integral_v<T>) {
        if (__int128(sites.back().x) - sites.front().x > 2000000000 ||
            __int128(scratch.back().y) - scratch.front().y > 2000000000)
            throw std::invalid_argument("各軸の座標幅は2*10^9以下である必要があります");
    } else {
        R dx = sites.back().x - sites.front().x, dy = scratch.back().y - scratch.front().y;
        if (!std::isfinite(dx * dx + dy * dy))
            throw std::invalid_argument("座標範囲の二乗距離が有限である必要があります");
        for (std::size_t i = 1; i < sites.size(); ++i) {
            dx = sites[i].x - sites[i - 1].x;
            dy = scratch[i].y - scratch[i - 1].y;
            if ((dx != 0 && dx * dx == 0) || (dy != 0 && dy * dy == 0))
                throw std::invalid_argument("座標差の二乗がアンダーフローしています");
        }
    }
    if (sites.size() < 2)
        return std::nullopt;
    for (std::size_t i = 1; i < sites.size(); ++i)
        if (sites[i - 1].x == sites[i].x && sites[i - 1].y == sites[i].y)
            return ClosestPairResult<R>{R(0), std::minmax(sites[i - 1].id, sites[i].id)};
    auto squared = [](Site a, Site b) {
        R dx = a.x - b.x, dy = a.y - b.y;
        return dx * dx + dy * dy;
    };
    ClosestPairResult<R> best{squared(sites[0], sites[1]), std::minmax(sites[0].id, sites[1].id)};
    auto consider = [&](Site a, Site b) {
        R distance = squared(a, b);
        if (distance < best.distanceSquared)
            best = {distance, std::minmax(a.id, b.id)};
    };
    auto solve = [&](auto &&self, Int first, Int last) -> void {
        if (last - first <= 3) {
            for (Int i = first; i < last; ++i)
                for (Int j = i + 1; j < last; ++j)
                    consider(sites[i], sites[j]);
            std::sort(sites.begin() + first, sites.begin() + last, yOrder);
            return;
        }
        Int mid = first + (last - first) / 2;
        R splitX = sites[mid].x;
        self(self, first, mid);
        self(self, mid, last);
        std::merge(sites.begin() + first, sites.begin() + mid, sites.begin() + mid, sites.begin() + last,
                   scratch.begin() + first, yOrder);
        std::copy(scratch.begin() + first, scratch.begin() + last, sites.begin() + first);
        Int count = 0;
        for (Int i = first; i < last; ++i) {
            Site p = sites[i];
            R dx = p.x - splitX;
            if (dx * dx >= best.distanceSquared)
                continue;
            for (Int j = count; j-- > 0;) {
                R dy = p.y - scratch[first + j].y;
                if (dy * dy >= best.distanceSquared)
                    break;
                consider(p, scratch[first + j]);
            }
            scratch[first + count++] = p;
        }
    };
    solve(solve, 0, sites.size());
    return best;
}
}
