#pragma once
#include <cplib/geometry/polygon.hpp>
#include <cplib/math/bigint.hpp>
#include <stdexcept>

namespace cplib {

struct LatticePolygonPoints {
    BigInt interior;
    BigInt boundary;
};

// 単純な格子多角形の内点数と境界点数を厳密に返す。周回方向は自由。O(N) 回の gcd と算術。
// 相異なる頂点を周回順に並べる（終点に始点を再掲してもよい）。面積 0 の多角形は対象外。
template <class T, class Id = void>
    requires(std::is_integral_v<T> && sizeof(T) <= sizeof(Int))
LatticePolygonPoints lattice_polygon_points(const Polygon<T, Id> &polygon) {
    BigInt twice_area{0}, boundary{0};
    for (std::size_t i = 0; i < polygon.size(); ++i) {
        const auto &a = polygon[i], &b = polygon[(i + 1) % polygon.size()];
        BigInt ax(a.x), ay(a.y), bx(b.x), by(b.y);
        twice_area += ax * by - ay * bx;

        // 絶対値を符号なし差で求めると、座標型の全域でも差が収まる。
        auto difference = [](T x, T y) -> UInt {
            return x < y ? UInt(y) - UInt(x) : UInt(x) - UInt(y);
        };
        boundary += BigInt(std::gcd(difference(a.x, b.x), difference(a.y, b.y)));
    }
    if (twice_area == 0)
        throw std::invalid_argument("lattice_polygon_points requires nonzero area");
    if (twice_area < 0)
        twice_area = -twice_area;
    return {(twice_area - boundary + BigInt(2)) / BigInt(2), std::move(boundary)};
}

}
