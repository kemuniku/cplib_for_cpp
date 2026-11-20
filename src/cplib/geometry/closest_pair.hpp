#pragma once
#include <cplib/geometry/base.hpp>
#include <stdexcept>
#include <tuple>

namespace cplib {
namespace detail {
template <class R> struct ClosestPairResult {
    R distanceSquared;
    std::pair<Int, Int> indices;
};

template <class T, class Id = void> requires (std::is_signed_v<T> || cplib::isFloatingPoint<T>)
auto closest_pair_result(const std::vector<Point<T, Id>> &points)
    -> ClosestPairResult<std::conditional_t<std::is_integral_v<T>, Int, std::common_type_t<T, double>>> {
    using R = std::conditional_t<std::is_integral_v<T>, Int, std::common_type_t<T, double>>;
    assert(points.size() >= 2);
    struct Site {
        R x, y;
        Int id;
    };
    std::vector<Site> sites;
    for (Int i = 0; i < Int(points.size()); ++i) {
        Site p{R(points[i].x), R(points[i].y), i};
        if constexpr (!std::is_integral_v<T>) {
            if (!isfinite(p.x) || !isfinite(p.y))
                throw std::invalid_argument("座標は有限である必要があります");
        }
        sites.push_back(p);
    }
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
        if (!isfinite(dx * dx + dy * dy))
            throw std::invalid_argument("座標範囲の二乗距離が有限である必要があります");
        for (std::size_t i = 1; i < sites.size(); ++i) {
            dx = sites[i].x - sites[i - 1].x;
            dy = scratch[i].y - scratch[i - 1].y;
            if ((dx != 0 && dx * dx == 0) || (dy != 0 && dy * dy == 0))
                throw std::invalid_argument("座標差の二乗がアンダーフローしています");
        }
    }
    for (std::size_t i = 1; i < sites.size(); ++i)
        if (sites[i - 1].x == sites[i].x && sites[i - 1].y == sites[i].y)
            return {R(0), {std::min(sites[i - 1].id, sites[i].id),
                           std::max(sites[i - 1].id, sites[i].id)}};
    auto squared = [](Site a, Site b) {
        R dx = a.x - b.x, dy = a.y - b.y;
        return dx * dx + dy * dy;
    };
    ClosestPairResult<R> best{squared(sites[0], sites[1]),
                              {std::min(sites[0].id, sites[1].id), std::max(sites[0].id, sites[1].id)}};
    auto consider = [&](Site a, Site b) {
        R distance = squared(a, b);
        if (distance < best.distanceSquared)
            best = {distance, {std::min(a.id, b.id), std::max(a.id, b.id)}};
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

// 最近点対の元の添字(i<j)を返す。2点未満は不可。時間 O(N log N)、領域 O(N)。
// 整数は各軸の幅2*10^9以下、浮動小数点は有限座標・有限二乗距離と非零差の二乗が非零であること。
// 整数はInt、浮動小数点はdouble以上で計算し、long doubleの精度を保つ。不正な座標はinvalid_argument。
template <class T, class Id = void> requires (std::is_signed_v<T> || cplib::isFloatingPoint<T>)
// 整数幅: 座標をIntへ変換後に計算。各軸幅D<=2*10^9なら二乗和<=2D^2<=8*10^18（64ビット）。
// |座標|<=10^9ならPoint<int>のまま渡せる。入力自体はIntに収まること。
std::pair<Int, Int> closest_pair(const std::vector<Point<T, Id>> &points) {
    return detail::closest_pair_result(points).indices;
}

// 最近点対の距離を返す。2点未満は不可。時間 O(N log N)、領域 O(N)。
template <class T, class Id = void> requires (std::is_signed_v<T> || cplib::isFloatingPoint<T>)
// 整数幅: 座標をIntへ変換後に計算。各軸幅D<=2*10^9なら二乗和<=2D^2<=8*10^18（64ビット）。
// |座標|<=10^9ならPoint<int>のまま渡せる。入力自体はIntに収まること。
auto closest_pair_distance(const std::vector<Point<T, Id>> &points) {
    using F = std::conditional_t<std::is_integral_v<T>, double, std::common_type_t<T, double>>;
    return sqrt(F(detail::closest_pair_result(points).distanceSquared));
}
}
