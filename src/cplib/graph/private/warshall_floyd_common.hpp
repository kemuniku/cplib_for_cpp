#pragma once
#include <cplib/graph/graph.hpp>
#include <cplib/graph/warshall_floyd_negative.hpp>
#include <cplib/utils/constants.hpp>

namespace cplib::detail {
template <class T> constexpr T distance_inf() {
    if constexpr (std::is_same_v<T, float>)
        return T(1e30);
    else if constexpr (std::is_floating_point_v<T>)
        return T(1e100);
    else if constexpr (sizeof(T) <= 4)
        return T(INF32);
    else
        return T(INF64);
}

struct ScalarWarshall {
    // 通常版の全点対最短路を求める。O(V^3)。
    template <bool Check, class T> static bool run(std::vector<std::vector<T>> &d, T zero, T inf) {
        for (std::size_t k = 0; k < d.size(); ++k) {
            for (std::size_t i = 0; i < d.size(); ++i)
                for (std::size_t j = 0; j < d.size(); ++j)
                    if (d[i][k] != inf && d[k][j] != inf)
                        d[i][j] = std::min(d[i][j], T(d[i][k] + d[k][j]));
            if constexpr (Check)
                for (std::size_t i = 0; i < d.size(); ++i)
                    if (d[i][i] < zero)
                        return true;
        }
        return false;
    }
};

template <class Kernel, bool Check, class T>
bool warshall_run(std::vector<std::vector<T>> &d, T zero, T inf) {
    for (const auto &row : d)
        assert(row.size() == d.size());
    for (std::size_t i = 0; i < d.size(); ++i) {
        if constexpr (Check)
            d[i][i] = std::min(d[i][i], zero);
        else
            d[i][i] = zero;
    }
    if constexpr (Check)
        for (std::size_t i = 0; i < d.size(); ++i)
            if (d[i][i] < zero)
                return true;
    if (d.empty())
        return false;
    return Kernel::template run<Check>(d, zero, inf);
}

template <class Kernel, bool Check, class T>
void warshall_inplace(std::vector<std::vector<T>> &d, T zero, T inf) {
    if (warshall_run<Kernel, Check>(d, zero, inf))
        warshall_floyd_negative_finish(
            d, zero, inf, [](auto &a, T z, T f) { return warshall_run<Kernel, true>(a, z, f); });
}

template <GraphTypes G, class T> auto graph_distance_matrix(const G &g, T zero, T inf) {
    std::vector<std::vector<T>> d(g.len, std::vector<T>(g.len, inf));
    for (Int i = 0; i < g.len; ++i) {
        d[i][i] = zero;
        for (auto [j, cost] : g.to_and_cost(i))
            d[i][j] = std::min(d[i][j], T(cost));
    }
    return d;
}
}
