#pragma once
#include <cplib/str/static_string.hpp>
#include <cplib/collections/waveletmatrix.hpp>
#include <cplib/utils/itertools.hpp>
#include <unordered_map>

namespace cplib {
namespace detail {
struct StaticSearchEntry {
    virtual ~StaticSearchEntry() = default;
};

inline std::unordered_map<const void *, std::shared_ptr<StaticSearchEntry>> staticSearchCache;
inline const void *staticSearchLastKey = nullptr;
inline std::shared_ptr<StaticSearchEntry> staticSearchLast;
}

template <class T> class StaticStringSearch {
    struct Cache {
        Int first = 0, last = 0, l = 0, r = 0, length = 0, count = 0;
    };

    struct Index {
        WaveletMatrix wm;
        std::vector<Cache> cache;
        Int shift = 0;

        explicit Index(const std::vector<Int> &positions) : wm(positions) {
            Int slots = positions.size();
            while (slots > 4096) {
                slots = (slots + 1) >> 1;
                ++shift;
            }
            cache.resize(slots);
        }
    };

    std::shared_ptr<Index> index;

    std::pair<Int, Int> suffixRange(const StaticString<T> &pattern) const {
        Int rank = base->RSA[pattern.l];
        return {base->RMQ.minLeft(rank, std::int32_t(pattern.len())),
                base->RMQ.maxRight(rank, std::int32_t(pattern.len())) + 1};
    }

public:
    std::shared_ptr<StaticStringBase<T>> base;

    explicit StaticStringSearch(std::shared_ptr<StaticStringBase<T>> b) : base(std::move(b)) {
        std::vector<Int> positions(base->SA.begin(), base->SA.end());
        index = std::make_shared<Index>(positions);
    }

    // 同じ基底のs内で重複を許したtの出現回数をO(log(N+2))で返す。空のtはs.len()+1回。
    // 最大4096件の結果を保持し、同じ対象・内容のキャッシュに当たればO(1)。Nは基底の長さ。
    Int count(const StaticString<T> &s, const StaticString<T> &t) const {
        assert(s.base == base && t.base == base);
        Int m = t.len();
        if (!m)
            return s.len() + 1;
        if (m > s.len())
            return 0;
        Int rank = base->RSA[t.l];
        auto &c = index->cache[rank >> index->shift];
        if (c.length == m && c.l == s.l && c.r == s.r && c.first <= rank && rank < c.last)
            return c.count;
        auto [first, last] = suffixRange(t);
        Int out = index->wm.range_freq(first, last, s.l, Int(s.r) - m + 1);
        c = {first, last, s.l, s.r, m, out};
        return out;
    }

    // 同じ基底のsにtが含まれるかO(log(N+2))で判定する。空のtは常に含まれる。
    bool contains(const StaticString<T> &s, const StaticString<T> &t) const {
        return count(s, t) > 0;
    }

    // Copy the small search handle into the coroutine so temporary handles are safe.
    static Generator<Int> positions(StaticStringSearch search, StaticString<T> s,
                                    StaticString<T> t) {
        assert(s.base == search.base && t.base == search.base);
        Int m = t.len();
        if (!m) {
            for (Int p = 0; p <= s.len(); ++p)
                co_yield p;
        } else if (m <= s.len()) {
            auto [first, last] = search.suffixRange(t);
            auto &wm = search.index->wm;
            Int begin = wm.range_lowerbound(first, last, s.l),
                end = wm.range_lowerbound(first, last, Int(s.r) - m + 1);
            for (Int k = begin; k < end; ++k)
                co_yield wm.kth_smallest(first, last, k) - s.l;
        }
    }

    // 同じ基底のs内で重複を許したtの出現位置を、sの先頭を0とする昇順で列挙する。
    // 準備と各要素の取得はO(log(N+2))、追加空間O(1)。空のtは0..s.len()を列挙する。
    auto findAll(StaticString<T> s, StaticString<T> t) const {
        return positions(*this, s, t);
    }

    struct View {
        StaticStringSearch search;
        StaticString<T> target;

        // 対象内で重複を許したtの出現回数をO(log(N+2))で返す。空のtは対象長+1回。
        Int count(const StaticString<T> &t) const {
            return search.count(target, t);
        }

        // 対象にtが含まれるかO(log(N+2))で判定する。空のtは常に含まれる。
        bool contains(const StaticString<T> &t) const {
            return count(t) > 0;
        }

        // 対象の先頭を0とする出現位置を昇順で列挙する。準備と各要素O(log(N+2))、追加空間O(1)。
        auto findAll(StaticString<T> t) const {
            return search.findAll(target, t);
        }
    };

    // 同じ基底のtargetと索引をO(1)で組にし、検索用ビューを作る。
    View operator[](StaticString<T> target) const {
        assert(target.base == base);
        return {*this, target};
    }
};

template <class T> struct StaticStringSearchCacheEntry : detail::StaticSearchEntry {
    StaticStringSearch<T> search;

    explicit StaticStringSearchCacheEntry(std::shared_ptr<StaticStringBase<T>> b)
        : search(std::move(b)) {
    }
};

// 全ての型・基底の索引とキャッシュの領域を手放す。保持済みの検索オブジェクトは引き続き使える。
inline void clearStaticStringSearchCache() {
    detail::staticSearchLastKey = nullptr;
    detail::staticSearchLast.reset();
    detail::staticSearchCache.clear();
}

// 指定した基底の索引をキャッシュから除く。保持済みの検索オブジェクトは引き続き使える。
template <class T>
void clearStaticStringSearchCache(const std::shared_ptr<StaticStringBase<T>> &b) {
    if (detail::staticSearchLastKey == b.get()) {
        detail::staticSearchLastKey = nullptr;
        detail::staticSearchLast.reset();
    }
    detail::staticSearchCache.erase(b.get());
}

// 基底ごとの索引を取得する。Nはb->S.size()。初回O(N log(N+2))、再取得は期待O(1)。
// 同じ基底の連続取得はテーブル検索を省きO(1)。
// clearStaticStringSearchCacheを呼ぶまでキャッシュが索引と基底を保持する。
template <class T>
StaticStringSearch<T> initStaticStringSearch(const std::shared_ptr<StaticStringBase<T>> &b) {
    auto key = b.get();
    if (detail::staticSearchLastKey == key && detail::staticSearchLast)
        return std::static_pointer_cast<StaticStringSearchCacheEntry<T>>(detail::staticSearchLast)
            ->search;
    auto it = detail::staticSearchCache.find(key);
    std::shared_ptr<detail::StaticSearchEntry> entry;
    if (it == detail::staticSearchCache.end()) {
        entry = std::make_shared<StaticStringSearchCacheEntry<T>>(b);
        detail::staticSearchCache[key] = entry;
    } else
        entry = it->second;
    detail::staticSearchLastKey = key;
    detail::staticSearchLast = entry;
    return std::static_pointer_cast<StaticStringSearchCacheEntry<T>>(entry)->search;
}

// 同じ基底のs内で重複を許したtの出現回数を返す。空のtはs.len()+1回。
// 索引が必要な初回はO(N log(N+2))、構築後は期待O(log(N+2))。Nは基底の長さ。
template <class T> Int count(const StaticString<T> &s, const StaticString<T> &t) {
    assert(s.base == t.base);
    if (!t.len())
        return s.len() + 1;
    if (t.len() > s.len())
        return 0;
    return initStaticStringSearch(s.base).count(s, t);
}

// 同じ基底のsにtが含まれるか判定する。空のtは常に含まれる。
// 索引が必要な初回はO(N log(N+2))、構築後は期待O(log(N+2))。
template <class T> bool contains(const StaticString<T> &s, const StaticString<T> &t) {
    return count(s, t) > 0;
}

// sの先頭を0とする出現位置を重複を許して昇順に列挙する。空のtは0..s.len()。
// 初回は必要なら索引を構築する。構築後の準備は期待O(log(N+2))、各要素O(log(N+2))、追加空間O(1)。
template <class T> auto findAll(StaticString<T> s, StaticString<T> t) {
    assert(s.base == t.base);
    return initStaticStringSearch(s.base).findAll(s, t);
}
}
