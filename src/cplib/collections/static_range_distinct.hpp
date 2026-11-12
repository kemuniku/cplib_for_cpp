#pragma once
#include <cplib/collections/waveletmatrix.hpp>
#include <functional>
#include <unordered_map>

namespace cplib {
class StaticRangeDistinct {
    WaveletMatrix wm_;
    Int n_;
    template <class T, class Compare>
    static std::vector<Int> previousByComparison(const std::vector<T> &values, Compare cmp) {
        std::vector<Int> order(values.size()), previous(values.size());
        std::iota(order.begin(), order.end(), 0);
        std::sort(order.begin(), order.end(), [&](Int a, Int b) {
            auto c = cmp(values[a], values[b]);
            return c < 0 || (c == 0 && a < b);
        });
        for (Int i = 1; i < Int(order.size()); ++i)
            if (cmp(values[order[i - 1]], values[order[i]]) == 0)
                previous[order[i]] = order[i - 1] + 1;
        return previous;
    }
    template <class T> static std::vector<Int> previous(const std::vector<T> &values) {
        if constexpr (requires(T value) { std::hash<T>{}(value); value == value; }) {
            std::unordered_map<T, Int> last;
            std::vector<Int> result(values.size());
            for (Int i = 0; i < Int(values.size()); ++i) {
                auto [it, inserted] = last.try_emplace(values[i], 0);
                result[i] = it->second;
                it->second = i + 1;
            }
            return result;
        } else {
            return previousByComparison(values, [](const T &a, const T &b) {
                return a < b ? -1 : b < a ? 1 : 0;
            });
        }
    }
public:
    // 静的な区間種類数の索引を構築する。期待時間 O(N log(N+1))、領域 O(N log(N+1)) ビットとO(N)ワード。
    template <class T>
    explicit StaticRangeDistinct(const std::vector<T> &values) : wm_(previous(values)), n_(values.size()) {
    }
    // 負・0・正を返す比較関数だけで構築する。時間 O(N log(N+1))、構築用領域 O(N) ワードを追加する。
    template <class T, class Compare>
    StaticRangeDistinct(const std::vector<T> &values, Compare cmp)
        : wm_(previousByComparison(values, cmp)), n_(values.size()) {
    }
    // 元の列の長さを返す。計算量 O(1)。
    Int size() const {
        return n_;
    }
    // 半開区間[l,r)の異なる値の数を返す。時間 O(log(N+1))、空区間は O(1)。
    Int count(Int l, Int r) const {
        assert(0 <= l && l <= r && r <= n_);
        return l == r ? 0 : wm_.range_lowerbound(l, r, l + 1);
    }
};
// 区間種類数の索引を構築する。期待時間 O(N log(N+1))、領域 O(N log(N+1)) ビットとO(N)ワード。
template <class T> auto initStaticRangeDistinct(const std::vector<T> &values) {
    return StaticRangeDistinct(values);
}
// 比較関数で区間種類数の索引を構築する。時間 O(N log(N+1))、領域 O(N log(N+1)) ビットとO(N)ワード。
template <class T, class Compare> auto initStaticRangeDistinct(const std::vector<T> &values, Compare cmp) {
    return StaticRangeDistinct(values, cmp);
}
}
