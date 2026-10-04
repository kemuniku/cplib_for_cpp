#pragma once
#include <cplib/collections/staticRMQ.hpp>
#include <cplib/str/suffix_array.hpp>

namespace cplib {
// 編集距離を返す。k超過は-1。N=|s|+|t|、O(N log N+k²)時間・O(N log N+k)領域。
// 挿入・削除・置換は各コスト1。k >= 0が必要で、各バイトを1文字として扱う。ハッシュは使わない。
inline Int editDistance(const std::string &s, const std::string &t, Int k) {
    assert(k >= 0);
    Int n = s.size(), m = t.size();
    if (std::abs(n - m) > k)
        return -1;
    if (!n || !m)
        return std::max(n, m);
    if (!k)
        return s == t ? 0 : -1;
    auto joined = s + t;
    auto sa = suffix_array(joined);
    std::vector<Int> rank(joined.size());
    for (Int i = 0; i < Int(sa.size()); ++i)
        rank[sa[i]] = i;
    auto rmq = initRMQ(lcp_array(joined, sa));
    auto extend = [&](Int x, Int y) {
        if (x == n || y == m)
            return Int(0);
        return std::min(
            {n - x, m - y,
             rmq.query(std::min(rank[x], rank[n + y]), std::max(rank[x], rank[n + y]))});
    };
    Int limit = std::min(k, std::max(n, m)), offset = limit + 1;
    std::vector<Int> previous(2 * limit + 3, -1), current(previous.size());
    previous[offset] = extend(0, 0);
    if (n == m && previous[offset] == n)
        return 0;
    for (Int edits = 1; edits <= limit; ++edits) {
        std::fill(current.begin(), current.end(), -1);
        for (Int d = -std::min(edits, n); d <= std::min(edits, m); ++d) {
            Int idx = offset + d, x = previous[idx];
            if (x >= 0 && x < n && x + d < m)
                ++x;
            Int insertion = previous[idx - 1];
            if (insertion >= 0 && insertion + d - 1 < m)
                x = std::max(x, insertion);
            Int deletion = previous[idx + 1];
            if (deletion >= 0 && deletion < n)
                x = std::max(x, deletion + 1);
            if (x < 0)
                continue;
            Int y = x + d;
            x += extend(x, y);
            current[idx] = x;
            if (d == m - n && x == n)
                return edits;
        }
        std::swap(previous, current);
    }
    return -1;
}
}
