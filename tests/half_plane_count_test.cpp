#include <cplib/geometry/half_plane_count.hpp>
#include <array>
#include <cassert>
#include <random>

using cplib::Int;
using cplib::Point;
using Wide = __int128;

struct HalfPlane {
    Int a, b;
    Wide c;
};

template <class T>
void check(const std::vector<Point<T>> &points, std::span<const HalfPlane> queries) {
    const cplib::PartitionTree<T> partition(points);
    for (const auto [a, b, c] : queries) {
        std::size_t expected = 0;
        for (const auto &p : points) {
            const Wide lhs = Wide(a) * p.x, rhs = Wide(b) * p.y;
            if (rhs > 0 && c < std::numeric_limits<Wide>::min() + rhs)
                continue;
            if (rhs < 0 && c > std::numeric_limits<Wide>::max() + rhs) {
                ++expected;
                continue;
            }
            expected += lhs <= c - rhs;
        }
        assert(partition.count(a, b, c) == expected);
    }
}

int main() {
    const std::array small_queries{
        HalfPlane{0, 0, -1}, HalfPlane{0, 0, 0}, HalfPlane{1, 0, 0}, HalfPlane{0, -1, 0},
        HalfPlane{1, -1, 0}, HalfPlane{-1, 1, 0}, HalfPlane{1, -1, -1}, HalfPlane{7, 2, 14},
        HalfPlane{1, 1, std::numeric_limits<Wide>::min()}, HalfPlane{1, 1, std::numeric_limits<Wide>::max()},
    };
    check<Int>({}, small_queries);
    check<Int>({{0, 0}}, small_queries);
    check<Int>(std::vector<Point<Int>>(1025, {0, 0}), small_queries);

    std::vector<Point<Int>> vertical, diagonal, convex, repeated_column;
    for (int i = -128; i <= 128; ++i) {
        vertical.push_back({0, i});
        diagonal.push_back({i, i});
        convex.push_back({i, i * i});
        repeated_column.push_back({i % 2 ? 0 : i, i});
    }
    check(vertical, small_queries);
    check(diagonal, small_queries);
    check(repeated_column, small_queries);
    std::vector<HalfPlane> tangencies;
    for (int i = -139; i <= 139; ++i) {
        tangencies.push_back({-2 * i, 1, -i * i});
        tangencies.push_back({-2 * i, 1, -i * i - 1});
        tangencies.push_back({2 * i, -1, i * i - 1});
    }
    check(convex, tangencies);

    constexpr Int lo = std::numeric_limits<Int>::min(), hi = std::numeric_limits<Int>::max();
    std::vector<Point<Int>> extreme;
    for (Int x : {lo, lo + 1, Int(-2), Int(0), Int(1), hi - 1, hi})
        for (Int y : {lo, lo + 1, Int(-2), Int(0), Int(1), hi - 1, hi})
            extreme.push_back({x, y});
    const std::array wide_queries{
        HalfPlane{lo, lo, std::numeric_limits<Wide>::max()}, HalfPlane{lo, lo, 0},
        HalfPlane{lo, hi, std::numeric_limits<Wide>::min()}, HalfPlane{hi, hi, Wide(hi) * hi},
        HalfPlane{hi, lo, -1}, HalfPlane{lo, 0, Wide(lo) * lo - 1},
        HalfPlane{0, hi, Wide(hi) * lo}, HalfPlane{-1, 1, 0}, HalfPlane{1, 1, 0},
    };
    check(extreme, wide_queries);
    // 分子が UINT64_MAX になる分割方向と、その逆向き。
    for (bool upward : {false, true}) {
        std::vector<Point<Int>> steep;
        for (int i = 0; i < 17; ++i) {
            steep.push_back({0, upward ? lo : hi});
            steep.push_back({1, upward ? hi : lo});
        }
        check(steep, wide_queries);
    }

    std::mt19937_64 rng(536423);
    std::uniform_int_distribution<Int> large(lo, hi);
    std::uniform_int_distribution<int> small(-12, 12);
    for (int trial = 0; trial < 320; ++trial) {
        std::vector<Point<Int>> points(31 + rng() % 550);
        for (auto &p : points) {
            if (trial % 4 == 0)
                p = {large(rng), large(rng)};
            else if (trial % 4 == 1)
                p = {small(rng) / 6, small(rng)};
            else if (trial % 4 == 2)
                p = {small(rng), small(rng)};
            else
                p = {large(rng), small(rng)};
        }
        std::vector<HalfPlane> queries;
        for (int i = 0; i < 80; ++i) {
            const Int a = (i % 2) ? small(rng) : large(rng), b = (i % 3) ? small(rng) : large(rng);
            auto p = points[rng() % points.size()];
            if (p.x == lo && p.y == lo)
                p.y = 0;
            const Wide c = i % 5 == 0 ? large(rng) : Wide(a) * p.x + Wide(b) * p.y;
            queries.push_back({a, b, c});
        }
        check(points, queries);
    }

    const std::vector<Point<int>> small_coordinates{{-2, 0}, {1, 1}, {1, 1}, {0, 2}};
    check(small_coordinates, small_queries);
}
