#pragma once
#include <cplib/str/fixedlength_merged_static_string.hpp>

namespace cplib {
template <class T> struct RepeatedStaticString {
    StaticString<T> period;
    Int length_ = 0;

    // 長さを返す。計算量 O(1)。
    Int size() const {
        return length_;
    }

    // i番目の要素を返す。計算量 O(1)。
    T operator[](Int i) const {
        assert(0 <= i && i < length_);
        return period[i % period.size()];
    }

    // 末尾基準の添字で要素を返す。計算量 O(1)。
    T operator[](BackwardsIndex i) const {
        return (*this)[length_ - i.value];
    }

    // 文字列に変換する。要素の出力時間の合計と O(長さ+1) が必要。
    std::string to_string() const {
        std::ostringstream out;
        for (Int i = 0; i < length_; ++i) {
            if constexpr (!std::is_same_v<T, char>) {
                if (i)
                    out << ' ';
            }
            out << (*this)[i];
        }
        return out.str();
    }
};

// sを無限に繰り返した列の先頭k要素を表す。kは反復回数ではなく要素数。
// kは非負。k > 0ならsは空でないこと。計算量 O(1)。
template <class T> auto initRepeatedStaticString(const StaticString<T> &s, Int k) {
    assert(k >= 0 && (k == 0 || s.size() > 0));
    return RepeatedStaticString<T>{s, k};
}

namespace detail {
// 無限反復列の最初の不一致は ST と TS の最初の不一致に等しい。O(1)。
template <class T> Int infinite_lcp(const StaticString<T> &s, const StaticString<T> &t) {
    assert(s.base == t.base && s.size() > 0 && t.size() > 0);
    Int result = lcp(s & t, t & s);
    return result == s.size() + t.size() ? std::numeric_limits<Int>::max() : result;
}

template <class S, class T> int repeated_cmp(const S &s, const T &t, Int common) {
    if (common == std::min(s.size(), t.size()))
        return (s.size() > t.size()) - (s.size() < t.size());
    return static_value_less(s[common], t[common]) ? -1 : 1;
}
}

// 同じ基底の反復列同士の最長共通接頭辞長を返す。計算量 O(1)。
template <class T> Int lcp(const RepeatedStaticString<T> &s, const RepeatedStaticString<T> &t) {
    assert(s.period.base == t.period.base);
    Int n = std::min(s.size(), t.size());
    return n ? std::min(n, detail::infinite_lcp(s.period, t.period)) : 0;
}

// 同じ基底の反復列と部分列の最長共通接頭辞長を返す。計算量 O(1)。
template <class T> Int lcp(const RepeatedStaticString<T> &s, const StaticString<T> &t) {
    assert(s.period.base == t.base);
    Int n = std::min(s.size(), t.size());
    return n ? std::min(n, detail::infinite_lcp(s.period, t)) : 0;
}

// 同じ基底の部分列と反復列の最長共通接頭辞長を返す。計算量 O(1)。
template <class T> Int lcp(const StaticString<T> &s, const RepeatedStaticString<T> &t) {
    return lcp(t, s);
}

// 同じ基底の反復列同士を辞書順に比較する。計算量 O(1)。
template <class T> int cmp(const RepeatedStaticString<T> &s, const RepeatedStaticString<T> &t) {
    return detail::repeated_cmp(s, t, lcp(s, t));
}

// 同じ基底の反復列と部分列を辞書順に比較する。計算量 O(1)。
template <class T> int cmp(const RepeatedStaticString<T> &s, const StaticString<T> &t) {
    return detail::repeated_cmp(s, t, lcp(s, t));
}

// 同じ基底の部分列と反復列を辞書順に比較する。計算量 O(1)。
template <class T> int cmp(const StaticString<T> &s, const RepeatedStaticString<T> &t) {
    return detail::repeated_cmp(s, t, lcp(s, t));
}

// 同じ基底の反復列同士が等しいか判定する。計算量 O(1)。
template <class T>
bool operator==(const RepeatedStaticString<T> &s, const RepeatedStaticString<T> &t) {
    return s.size() == t.size() && lcp(s, t) == s.size();
}

// 同じ基底の反復列同士を辞書順に比較する。計算量 O(1)。
template <class T>
auto operator<=>(const RepeatedStaticString<T> &s, const RepeatedStaticString<T> &t) {
    return cmp(s, t) <=> 0;
}

// 同じ基底の反復列と部分列が等しいか判定する。計算量 O(1)。
template <class T> bool operator==(const RepeatedStaticString<T> &s, const StaticString<T> &t) {
    return s.size() == t.size() && lcp(s, t) == s.size();
}

// 同じ基底の反復列と部分列を辞書順に比較する。計算量 O(1)。
template <class T> auto operator<=>(const RepeatedStaticString<T> &s, const StaticString<T> &t) {
    return cmp(s, t) <=> 0;
}

// 反復列を文字列に変換する。要素の出力時間の合計と O(長さ+1) が必要。
template <class T> std::string to_string(const RepeatedStaticString<T> &s) {
    return s.to_string();
}
}
