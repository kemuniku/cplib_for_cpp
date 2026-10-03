#pragma once
#include <cplib/common.hpp>

namespace cplib {
namespace detail {
template <class A, class B> auto lcs_table(const A &a, const B &b) {
    std::vector<std::vector<Int>> dp(b.size() + 1, std::vector<Int>(a.size()));
    for (std::size_t i = 0; i < b.size(); ++i) {
        Int now = 0;
        for (std::size_t j = 0; j < a.size(); ++j) {
            if (a[j] == b[i]) {
                Int tmp = dp[i][j];
                dp[i + 1][j] = now + 1;
                if (tmp > now)
                    now = tmp;
            } else {
                dp[i + 1][j] = dp[i][j];
                if (dp[i][j] > now)
                    now = dp[i][j];
            }
        }
    }
    return dp;
}
}

// 最長共通部分列長を求める。O(|a||b|)時間・領域。
template <class A, class B> Int LCS(const A &a, const B &b) {
    if (a.empty() || b.empty())
        return 0;
    auto dp = detail::lcs_table(a, b);
    return *std::max_element(dp.back().begin(), dp.back().end());
}

// 元実装と同じ後方走査でLCSを復元する。O(|a||b|)時間・領域。
template <class A, class B> auto restoreLCS(const A &a, const B &b) {
    std::vector<typename A::value_type> ans;
    if (a.empty() || b.empty())
        return ans;
    auto dp = detail::lcs_table(a, b);
    Int now = std::max_element(dp.back().begin(), dp.back().end()) - dp.back().begin();
    for (Int i = b.size(); i >= 1; --i) {
        if (dp[i - 1][now] == dp[i][now])
            continue;
        for (Int j = now - 1; j >= 0; --j)
            if (dp[i - 1][j] == dp[i][now] - 1) {
                now = j;
                break;
            }
        ans.push_back(b[i - 1]);
    }
    std::reverse(ans.begin(), ans.end());
    return ans;
}
}
