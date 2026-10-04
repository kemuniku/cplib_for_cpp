#pragma once
#include <cplib/utils/constants.hpp>
#include <tuple>

namespace cplib {
using KnapsackItem = std::pair<Int, Int>;
using BoundedKnapsackItem = std::tuple<Int, Int, Int>;

// 0/1ナップサックを容量DPで解く。O(NW)。
// 各品物(v,w)を最大一回選び、総重量がW以下となる最大価値を返す。
inline Int solve_01knapsack_NW(std::span<const KnapsackItem> items, Int W) {
    std::vector<Int> dp(W + 1, -INF64);
    dp[0] = 0;
    for (auto [v, w] : items)
        for (Int j = W - w; j >= 0; --j)
            dp[j + w] = std::max(dp[j + w], dp[j] + v);
    return *std::max_element(dp.begin(), dp.end());
}

// 0/1ナップサックを価値DPで解く。O(NΣv)。
inline Int solve_01knapsack_NV(std::span<const KnapsackItem> items, Int W) {
    Int V = 0;
    for (auto [v, w] : items)
        V += v;
    std::vector<Int> dp(V + 1, INF64);
    dp[0] = 0;
    for (auto [v, w] : items)
        for (Int j = V - v; j >= 0; --j)
            dp[j + v] = std::min(dp[j + v], dp[j] + w);
    for (Int i = V; i >= 0; --i)
        if (dp[i] <= W)
            return i;
    return 0;
}

// 半分全列挙で0/1ナップサックを解く。O(N 2^(N/2))。
inline Int solve_01knapsack_meet_in_middle(std::span<const KnapsackItem> items, Int W) {
    auto enumerate = [](std::span<const KnapsackItem> part) {
        std::vector<KnapsackItem> r(std::size_t(1) << part.size());
        for (std::size_t bit = 1; bit < r.size(); ++bit) {
            int i = std::bit_width(bit) - 1;
            auto [v, w] = r[bit ^ (std::size_t(1) << i)];
            r[bit] = {v + part[i].first, w + part[i].second};
        }
        return r;
    };
    auto a = enumerate(items.first(items.size() / 2)),
         b = enumerate(items.subspan(items.size() / 2));
    for (auto &e : b)
        std::swap(e.first, e.second);
    std::sort(b.begin(), b.end());
    for (std::size_t i = 1; i < b.size(); ++i)
        b[i].second = std::max(b[i].second, b[i - 1].second);
    Int ans = -INF64;
    for (auto [v, w] : a)
        if (w <= W) {
            auto it = std::lower_bound(b.begin(), b.end(), KnapsackItem{W - w, INF64});
            --it;
            ans = std::max(ans, it->second + v);
        }
    return ans;
}

// 個数無制限のナップサックを容量DPで解く。O(NW)。
inline Int solve_UBknapsack_NW(std::span<const KnapsackItem> items, Int W) {
    std::vector<Int> dp(W + 1, -INF64);
    dp[0] = 0;
    for (auto [v, w] : items)
        for (Int j = 0; j <= W - w; ++j)
            dp[j + w] = std::max(dp[j + w], dp[j] + v);
    return *std::max_element(dp.begin(), dp.end());
}

// 個数制限付きナップサックをブロック内の前後走査で解く。O(NW)。
// 各品物(v,w,m)を最大m個選ぶ。個数0の品物は無視し、総重量W以下の最大価値を返す。
inline Int solve_BoundedKnapsack(std::span<const BoundedKnapsackItem> items, Int W) {
    std::vector<Int> dp(W + 1);
    for (auto [v, w, m0] : items) {
        if (m0 == 0 || w > W)
            continue;
        if (w == 0) {
            for (Int &x : dp)
                x += v * m0;
            continue;
        }
        Int m = std::min(m0, W / w);
        auto buf = dp;
        for (Int s = 0; s * w <= W; s += m) {
            Int l = s * w, r = std::min(W + 1, (s + m) * w);
            for (Int i = l; i < r - w; ++i)
                dp[i + w] = std::max(dp[i + w], dp[i] + v);
            for (Int i = r - w - 1; i >= l; --i)
                buf[i] = std::max(buf[i], buf[i + w] - v);
        }
        for (Int i = w * m; i <= W; ++i)
            dp[i] = std::max(dp[i], buf[i - w * m] + v * m);
    }
    return dp[W];
}
}
