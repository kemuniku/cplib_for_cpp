#pragma once
#include <cplib/str/private/merged_string_base.hpp>

namespace cplib {
template <class T> struct MergedStaticString {
    using value_type = T;
    std::shared_ptr<StaticStringBase<T>> base;
    std::vector<std::int32_t> L, R;

    // 同じ基底の半開区間を末尾に追加する。償却 O(1)。
    void addRange(std::shared_ptr<StaticStringBase<T>> b, std::int32_t l, std::int32_t r) {
        assert(l <= r);
        if (L.empty())
            base = b;
        else
            assert(base == b);
        L.push_back(l);
        R.push_back(r);
    }

    // 同じ基底の部分列を末尾に追加する。償却 O(1)。
    MergedStaticString &operator&=(const StaticString<T> &value) {
        addRange(value.base, value.l, value.r);
        return *this;
    }

    // 結合した区間の全長を返す。結合数Nとして O(N+1)。
    Int size() const {
        return detail::merged_length(*this);
    }

    // i番目の要素を返す。結合数Nとして O(N+1)。
    T operator[](Int i) const {
        return detail::merged_at(*this, i);
    }

    // 末尾基準の添字で要素を返す。結合数Nとして O(N+1)。
    T operator[](BackwardsIndex i) const {
        return (*this)[size() - i.value];
    }

    // 両端を含む区間のビューを返す。結合数Nとして O(N+1)。
    template <class A, class B> MergedStaticString operator[](ClosedSlice<A, B> slice) const {
        Int a = resolve_index(size(), slice.a), b = resolve_index(size(), slice.b), tmp = 0;
        MergedStaticString out;
        for (std::size_t i = 0; i < L.size(); ++i) {
            Int next = tmp + R[i] - L[i];
            if (tmp < a) {
                if (b < next)
                    out.addRange(base, L[i] + a - tmp, L[i] + b - tmp + 1);
                else if (a < next)
                    out.addRange(base, L[i] + a - tmp, R[i]);
            } else if (next <= b)
                out.addRange(base, L[i], R[i]);
            else if (tmp <= b)
                out.addRange(base, L[i], L[i] + b - tmp + 1);
            tmp = next;
        }
        return out;
    }

    // 文字列に変換する。結合数N、全長Lとして O(N+L+1) と要素の出力費用が必要。
    std::string to_string() const {
        return detail::merged_string(*this);
    }
};

// 同じ基底のN個の部分列を結合したビューを作る。O(N+1)。
template <class T> auto initMergedStaticString(const std::vector<StaticString<T>> &s) {
    MergedStaticString<T> out;
    for (const auto &v : s)
        out &= v;
    return out;
}

// 同じ基底のN個の区間を結合したビューを作る。O(N+1)。
template <class T>
auto initMergedStaticString(const StaticString<T> &s,
                            const std::vector<std::pair<Int, Int>> &ranges) {
    MergedStaticString<T> out;
    out.base = s.base;
    for (auto [l, r] : ranges) {
        assert(0 <= l && l <= r && r <= s.size());
        out.addRange(s.base, s.l + l, s.l + r);
    }
    return out;
}

// 部分列を末尾に結合した新しいビューを返す。元の結合数Nとして O(N+1)。
template <class T> auto operator&(MergedStaticString<T> s, const StaticString<T> &t) {
    s &= t;
    return s;
}

// StaticString & StaticString は固定長版と戻り値だけが違うため、動的版は名前付き関数。
// 同じ基底の2つの部分列を結合したビューを返す。O(1)。
template <class T> auto mergeStaticStrings(const StaticString<T> &s, const StaticString<T> &t) {
    assert(s.base == t.base);
    MergedStaticString<T> out;
    out &= s;
    out &= t;
    return out;
}

// 同じ基底の結合列の最長共通接頭辞長を返す。結合数の和をNとして O(N+1)。
template <class T> Int lcp(const MergedStaticString<T> &s, const MergedStaticString<T> &t) {
    return detail::merged_lcp(s, t);
}

// 同じ基底の結合列を辞書順に比較する。結合数の和をNとして O(N+1)。
template <class T> int cmp(const MergedStaticString<T> &s, const MergedStaticString<T> &t) {
    return detail::merged_cmp(s, t);
}

// 同じ基底の結合列の等値を判定する。結合数の和をNとして O(N+1)。
template <class T> bool operator==(const MergedStaticString<T> &s, const MergedStaticString<T> &t) {
    return s.size() == t.size() && lcp(s, t) == s.size();
}

// 同じ基底の結合列を辞書順に比較する。結合数の和をNとして O(N+1)。
template <class T>
auto operator<=>(const MergedStaticString<T> &s, const MergedStaticString<T> &t) {
    return cmp(s, t) <=> 0;
}

// 文字列に変換する。結合数N、全長Lとして O(N+L+1) と要素の出力費用が必要。
template <class T> std::string to_string(const MergedStaticString<T> &s) {
    return s.to_string();
}
}
