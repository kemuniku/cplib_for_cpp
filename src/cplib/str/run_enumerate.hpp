#pragma once
#include <cplib/str/zalgorithm.hpp>
#include <tuple>

namespace cplib {
// 分割統治と Z 配列により最大周期区間を列挙。(最小周期,l,r)、半開区間。
template <class Range> std::vector<std::tuple<Int, Int, Int>> run_enumerate(const Range &a) {
    using T = typename Range::value_type;
    Int n = a.size();
    std::vector<std::tuple<Int, Int, Int>> raw;
    auto dfs = [&](auto &&self, Int l, Int r) -> void {
        if (r - l <= 1)
            return;
        Int m = (l + r) / 2;
        self(self, l, m);
        self(self, m, r);
        std::vector<T> sl, sr;
        for (Int i = m; i-- > l;)
            sl.push_back(a[i]);
        for (Int i = r; i-- > l;)
            sl.push_back(a[i]);
        for (Int i = m; i < r; ++i)
            sr.push_back(a[i]);
        for (Int i = l; i < r; ++i)
            sr.push_back(a[i]);
        auto zsl = zalgorithm(sl), zsr = zalgorithm(sr);
        auto record = [&](Int ml, Int mr, Int p) {
            if (mr - ml >= 2 * p && (ml == 0 || a[ml - 1] != a[ml + p - 1]) &&
                (mr == n || a[mr] != a[mr - p]))
                raw.emplace_back(ml, mr, p);
        };
        for (Int p = 1; p <= m - l; ++p)
            record(std::max(l, m - p - zsl[p]), std::min(r, m + zsr[r - l - p]), p);
        for (Int p = 1; p <= r - m; ++p)
            record(std::max(l, m - zsl[r - l - p]), std::min(r, m + p + zsr[p]), p);
    };
    dfs(dfs, 0, n);
    std::sort(raw.begin(), raw.end());
    Int previous_l = -1, previous_r = -1;
    std::vector<std::tuple<Int, Int, Int>> result;
    for (auto [l, r, p] : raw) {
        if (l == previous_l && r == previous_r)
            continue;
        result.emplace_back(p, l, r);
        previous_l = l;
        previous_r = r;
    }
    std::sort(result.begin(), result.end());
    return result;
}

template <class Range> auto RunEnumerate(const Range &s) {
    return run_enumerate(s);
}
}
