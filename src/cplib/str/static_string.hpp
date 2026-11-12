#pragma once
#include <cplib/utils/backwards_index.hpp>
#include <cplib/str/suffix_array.hpp>
#include <cplib/collections/staticRMQ.hpp>
#include <memory>
#include <sstream>

namespace cplib {
namespace detail {
template <class T> bool static_value_less(const T &a, const T &b) {
    if constexpr (std::is_same_v<T, char>)
        return static_cast<unsigned char>(a) < static_cast<unsigned char>(b);
    else
        return a < b;
}
}

// N要素の接尾辞配列を作る。バイト型は O(N+256)、その他は O(N log(N+1))。
template <class T> std::vector<Int> genericSuffixArray(const std::vector<T> &s) {
    return suffix_array(s);
}

template <class T> struct StaticStringBase {
    std::vector<T> S;
    StaticRMQ<std::int32_t> RMQ;
    std::vector<std::int32_t> SA, RSA, LCP;
    std::int32_t size = 0;
    bool reversible = false;

    // 長さNの基底と検索用索引を作る。比較を O(1) として計算量 O(N log(N+1)+1)。
    StaticStringBase(std::vector<T> s, bool rev = false)
        : S(std::move(s)), size(S.size()), reversible(rev) {
        if (rev) {
            auto copy = S;
            S.insert(S.end(), copy.rbegin(), copy.rend());
        }
        auto sa = genericSuffixArray(S);
        SA.assign(sa.begin(), sa.end());
        RSA.resize(S.size());
        for (Int i = 0; i < Int(S.size()); ++i)
            RSA[SA[i]] = i;
        if (!S.empty()) {
            auto lcp = cplib::lcp_array(S, sa);
            LCP.assign(lcp.begin(), lcp.end());
        }
        RMQ = initRMQ(LCP);
    }

    // 検索列に対する接尾辞配列の下限または上限を返す。基底長N、検索列長Mとして O((M+1) log(N+2))。
    // 基底は反転非対応で構築し、要素の比較を O(1) とする。
    template <class Range> Int suffix_bound(const Range &s, bool upper) const {
        assert(!reversible);
        auto compare = [&](Int x) {
            for (Int i = 0; i < Int(s.size()); ++i) {
                if (i + x >= size)
                    return -1;
                if (detail::static_value_less(S[i + x], s[i]))
                    return -1;
                if (detail::static_value_less(s[i], S[i + x]))
                    return 1;
            }
            return 0;
        };
        Int l = 0, r = SA.size();
        while (l < r) {
            Int m = (l + r) / 2, c = compare(SA[m]);
            if (c < 0 || (upper && c == 0))
                l = m + 1;
            else
                r = m;
        }
        return l;
    }

    // 検索列に対する接尾辞配列の下限を返す。基底長N、検索列長Mとして O((M+1) log(N+2))。
    // 基底は反転非対応で構築し、要素の比較を O(1) とする。
    template <class Range> Int suffix_lowerbound(const Range &s) const {
        return suffix_bound(s, false);
    }

    // 検索列に対する接尾辞配列の上限を返す。基底長N、検索列長Mとして O((M+1) log(N+2))。
    // 基底は反転非対応で構築し、要素の比較を O(1) とする。
    template <class Range> Int suffix_upperbound(const Range &s) const {
        return suffix_bound(s, true);
    }

    // 検索列を接頭辞に持つ接尾辞数を返す。基底長N、検索列長Mとして O((M+1) log(N+2))。
    // 基底は反転非対応で構築し、要素の比較を O(1) とする。
    template <class Range> Int count(const Range &s) const {
        return suffix_upperbound(s) - suffix_lowerbound(s);
    }
};

// 長さNの基底と索引を作る。比較を O(1) として計算量 O(N log(N+1)+1)。
template <class T> auto initStaticStringBase(const std::vector<T> &s, bool reversible = false) {
    return std::make_shared<StaticStringBase<T>>(s, reversible);
}

// 長さNの文字列から基底と索引を作る。計算量 O(N log(N+1)+1)。
inline auto initStaticStringBase(const std::string &s, bool reversible = false) {
    return initStaticStringBase(std::vector<char>(s.begin(), s.end()), reversible);
}

namespace detail {
template <class T>
Int static_lcp_range(const std::shared_ptr<StaticStringBase<T>> &base, Int sl, Int sr, Int tl,
                     Int tr) {
    Int result = std::min(sr - sl, tr - tl);
    if (!result)
        return 0;
    Int l = base->RSA[sl], r = base->RSA[tl];
    if (l > r)
        std::swap(l, r);
    else if (l == r)
        return result;
    return std::min(result, Int(base->RMQ.query(l, r)));
}
}

template <class T> struct StaticString {
    std::shared_ptr<StaticStringBase<T>> base;
    std::int32_t l = 0, r = 0;

    // 部分列の長さを返す。計算量 O(1)。
    Int len() const {
        return Int(r) - l;
    }

    // 部分列の長さを返す。計算量 O(1)。
    Int size() const {
        return len();
    }

    // 先頭からi番目の要素を返す。計算量 O(1)。
    T operator[](Int i) const {
        assert(0 <= i && i < len());
        return base->S[l + i];
    }

    // 末尾基準の添字で要素を返す。計算量 O(1)。
    T operator[](BackwardsIndex i) const {
        return (*this)[len() - i.value];
    }

    // 半開区間の部分列ビューを返す。計算量 O(1)。
    StaticString substr(Int a, Int b) const {
        assert(0 <= a && a <= b && b <= len());
        return {base, std::int32_t(l + a), std::int32_t(l + b)};
    }

    // 両端を含む区間の部分列ビューを返す。計算量 O(1)。
    template <class L, class R> StaticString operator[](ClosedSlice<L, R> s) const {
        return substr(resolve_index(len(), s.a), resolve_index(len(), s.b) + 1);
    }

    // 反転した部分列のビューを返す。基底は反転対応で構築すること。計算量 O(1)。
    StaticString reversed() const {
        assert(base->reversible);
        return {base, 2 * base->size - r, 2 * base->size - l};
    }

    // 部分列が回文か判定する。基底は反転対応で構築すること。計算量 O(1)。
    bool isPalindrome() const {
        auto other = reversed();
        return detail::static_lcp_range(base, l, r, other.l, other.r) == len();
    }

    // 部分列を文字列に変換する。要素の出力時間の合計と O(長さ+1) が必要。
    std::string to_string() const {
        std::ostringstream out;
        for (Int i = 0; i < len(); ++i) {
            if constexpr (!std::is_same_v<T, char>) {
                if (i)
                    out << ' ';
            }
            out << (*this)[i];
        }
        return out.str();
    }
};

// 長さNの列から静的文字列と索引を作る。比較を O(1) として O(N log(N+1)+1)。
template <class T> auto toStaticString(const std::vector<T> &s, bool reversible = false) {
    return StaticString<T>{initStaticStringBase(s, reversible), 0, std::int32_t(s.size())};
}

// 長さNの文字列から静的文字列と索引を作る。計算量 O(N log(N+1)+1)。
inline auto toStaticString(const std::string &s, bool reversible = false) {
    return StaticString<char>{initStaticStringBase(s, reversible), 0, std::int32_t(s.size())};
}

// 部分列の長さを返す。計算量 O(1)。
template <class T> Int len(const StaticString<T> &s) {
    return s.len();
}

// 同じ基底の部分列の最長共通接頭辞長を返す。計算量 O(1)。
template <class T> Int lcp(const StaticString<T> &s, const StaticString<T> &t) {
    assert(s.base == t.base);
    return detail::static_lcp_range(s.base, s.l, s.r, t.l, t.r);
}

// 反転した部分列ビューを返す。基底は反転対応で構築すること。計算量 O(1)。
template <class T> auto reversed(const StaticString<T> &s) {
    return s.reversed();
}

// 部分列が回文か判定する。基底は反転対応で構築すること。計算量 O(1)。
template <class T> bool isPalindrome(const StaticString<T> &s) {
    return s.isPalindrome();
}

// 同じ基底の部分列の最長共通接尾辞長を返す。基底は反転対応で構築すること。計算量 O(1)。
template <class T> Int lcs(const StaticString<T> &s, const StaticString<T> &t) {
    assert(s.base == t.base && s.base->reversible);
    return lcp(s.reversed(), t.reversed());
}

// 同じ基底の部分文字列を RMQ と接尾辞順位で O(1) 比較。
template <class T> int cmp(const StaticString<T> &s, const StaticString<T> &t) {
    assert(s.base == t.base);
    Int n = std::min(s.len(), t.len());
    if (!n || s.l == t.l)
        return (s.len() > t.len()) - (s.len() < t.len());
    Int a = s.base->RSA[s.l], b = s.base->RSA[t.l];
    if (s.base->RMQ.query(std::min(a, b), std::max(a, b)) >= n)
        return (s.len() > t.len()) - (s.len() < t.len());
    return a < b ? -1 : 1;
}

// 同じ基底の部分列が等しいか判定する。計算量 O(1)。
template <class T> bool operator==(const StaticString<T> &s, const StaticString<T> &t) {
    return s.len() == t.len() && lcp(s, t) == s.len();
}

// 同じ基底の部分列を辞書順に比較する。計算量 O(1)。
template <class T> auto operator<=>(const StaticString<T> &s, const StaticString<T> &t) {
    return cmp(s, t) <=> 0;
}

// 同じ基底の部分列が指定接頭辞で始まるか判定する。計算量 O(1)。
template <class T> bool startsWith(const StaticString<T> &s, const StaticString<T> &prefix) {
    return lcp(s, prefix) == prefix.len();
}

// 部分列を文字列に変換する。要素の出力時間の合計と O(長さ+1) が必要。
template <class T> std::string to_string(const StaticString<T> &s) {
    return s.to_string();
}

// 基底長 M、要素数 N に対し O(N log(M+2))。順位境界の探索後に安定基数ソート。
template <class T> void sortStaticStrings(std::span<StaticString<T>> strings) {
    if (strings.empty())
        return;
    auto base = strings[0].base;

    struct Key {
        std::int32_t rank, length;
        std::size_t index;
    };

    std::vector<Key> keys(strings.size()), buffer(keys.size());
    for (std::size_t i = 0; i < strings.size(); ++i) {
        const auto &s = strings[i];
        assert(s.base == base);
        Int left = 0;
        if (s.len()) {
            Int rank = base->RSA[s.l], right = rank;
            while (left < right) {
                Int mid = (left + right) / 2;
                if (base->RMQ.query(mid, rank) >= s.len())
                    right = mid;
                else
                    left = mid + 1;
            }
            ++left;
        }
        keys[i] = {std::int32_t(left), std::int32_t(s.len()), i};
    }
    for (Int field = 0; field < 2; ++field)
        for (Int shift = 0; shift <= 24; shift += 8) {
            std::array<std::size_t, 256> counts{};
            for (auto key : keys)
                ++counts[((field ? key.rank : key.length) >> shift) & 255];
            std::size_t total = 0;
            for (auto &c : counts) {
                auto count = c;
                c = total;
                total += count;
            }
            for (auto key : keys) {
                auto digit = ((field ? key.rank : key.length) >> shift) & 255;
                buffer[counts[digit]++] = key;
            }
            keys.swap(buffer);
        }
    std::vector<StaticString<T>> out;
    out.reserve(keys.size());
    for (auto key : keys)
        out.push_back(strings[key.index]);
    std::copy(out.begin(), out.end(), strings.begin());
}

// 同じ基底のN個の部分列を辞書順に並べる。基底長Mとして O(N log(M+2))。
template <class T> void sortStaticStrings(std::vector<StaticString<T>> &s) {
    sortStaticStrings<T>(std::span<StaticString<T>>(s));
}

// 基底の全接尾辞を辞書順に返す。長さNとして反転対応なら O(N log(N+1)+1)、それ以外は O(N+1)。
template <class T> auto initSuffixArray(const std::shared_ptr<StaticStringBase<T>> &base) {
    std::vector<Int> sa;
    if (base->reversible)
        sa = genericSuffixArray(std::vector<T>(base->S.begin(), base->S.begin() + base->size));
    else
        sa.assign(base->SA.begin(), base->SA.end());
    std::vector<StaticString<T>> out;
    for (Int x : sa)
        out.push_back({base, std::int32_t(x), base->size});
    return out;
}

// 長さNの部分列の全接尾辞を辞書順に返す。比較を O(1) として O(N log(N+1)+1)。
template <class T> auto initSuffixArray(const StaticString<T> &s) {
    auto sa = genericSuffixArray(std::vector<T>(s.base->S.begin() + s.l, s.base->S.begin() + s.r));
    std::vector<StaticString<T>> out;
    for (Int x : sa)
        out.push_back({s.base, std::int32_t(x + s.l), s.r});
    return out;
}

// 複数の文字列を共通基底にまとめる。文字列数Kと全長Lに対し O((L+K) log(L+K+1)+1)。
inline auto toStaticStrings(const std::vector<std::string> &strings, bool reversible = false) {
    std::string tmp;
    for (const auto &s : strings) {
        tmp += s;
        tmp += '$';
    }
    auto base = initStaticStringBase(tmp, reversible);
    std::vector<StaticString<char>> out;
    std::int32_t now = 0;
    for (const auto &s : strings) {
        out.push_back({base, now, std::int32_t(now + s.size())});
        now += s.size() + 1;
    }
    return out;
}

// 反転非対応の基底について、接尾辞配列の下限を求める。
// 基底長N、検索列長M、要素の比較を O(1) として O((M+1) log(N+2))。
template <class T, class Range>
Int suffix_lowerbound(const std::shared_ptr<StaticStringBase<T>> &base, const Range &s) {
    return base->suffix_lowerbound(s);
}

// 反転非対応の基底について、接尾辞配列の上限を求める。
// 基底長N、検索列長M、要素の比較を O(1) として O((M+1) log(N+2))。
template <class T, class Range>
Int suffix_upperbound(const std::shared_ptr<StaticStringBase<T>> &base, const Range &s) {
    return base->suffix_upperbound(s);
}

// 反転非対応の基底について、検索列を接頭辞に持つ接尾辞数を返す。
// 基底長N、検索列長M、要素の比較を O(1) として O((M+1) log(N+2))。
template <class T, class Range>
Int count(const std::shared_ptr<StaticStringBase<T>> &base, const Range &s) {
    return base->count(s);
}
}
