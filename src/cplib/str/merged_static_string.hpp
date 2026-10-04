#pragma once
#include <cplib/str/private/merged_string_base.hpp>

namespace cplib {
template <class T> struct MergedStaticString {
    using value_type = T;
    std::shared_ptr<StaticStringBase<T>> base;
    std::vector<std::int32_t> L, R;

    void addRange(std::shared_ptr<StaticStringBase<T>> b, std::int32_t l, std::int32_t r) {
        assert(l <= r);
        if (L.empty())
            base = b;
        else
            assert(base == b);
        L.push_back(l);
        R.push_back(r);
    }

    MergedStaticString &operator&=(const StaticString<T> &value) {
        addRange(value.base, value.l, value.r);
        return *this;
    }

    // 計算量が O(結合数) である点に注意！
    Int len() const {
        return detail::merged_length(*this);
    }

    // 要素取得の計算量はO(結合数)である点に注意。
    T operator[](Int i) const {
        return detail::merged_at(*this, i);
    }

    T operator[](BackwardsIndex i) const {
        return (*this)[len() - i.value];
    }

    template <class A, class B> MergedStaticString operator[](ClosedSlice<A, B> slice) const {
        Int a = resolve_index(len(), slice.a), b = resolve_index(len(), slice.b), tmp = 0;
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

    std::string to_string() const {
        return detail::merged_string(*this);
    }
};

template <class T> auto initMergedStaticString(const std::vector<StaticString<T>> &s) {
    MergedStaticString<T> out;
    for (const auto &v : s)
        out &= v;
    return out;
}

template <class T>
auto initMergedStaticString(const StaticString<T> &s,
                            const std::vector<std::pair<Int, Int>> &ranges) {
    MergedStaticString<T> out;
    out.base = s.base;
    for (auto [l, r] : ranges) {
        assert(0 <= l && l <= r && r <= s.len());
        out.addRange(s.base, s.l + l, s.l + r);
    }
    return out;
}

template <class T> auto operator&(MergedStaticString<T> s, const StaticString<T> &t) {
    s &= t;
    return s;
}

// StaticString & StaticString は固定長版と戻り値だけが違うため、動的版は名前付き関数。
template <class T> auto mergeStaticStrings(const StaticString<T> &s, const StaticString<T> &t) {
    assert(s.base == t.base);
    MergedStaticString<T> out;
    out &= s;
    out &= t;
    return out;
}

template <class T> Int lcp(const MergedStaticString<T> &s, const MergedStaticString<T> &t) {
    return detail::merged_lcp(s, t);
}

template <class T> int cmp(const MergedStaticString<T> &s, const MergedStaticString<T> &t) {
    return detail::merged_cmp(s, t);
}

template <class T> bool operator==(const MergedStaticString<T> &s, const MergedStaticString<T> &t) {
    return s.len() == t.len() && lcp(s, t) == s.len();
}

template <class T>
auto operator<=>(const MergedStaticString<T> &s, const MergedStaticString<T> &t) {
    return cmp(s, t) <=> 0;
}

template <class T> Int len(const MergedStaticString<T> &s) {
    return s.len();
}

template <class T> std::string to_string(const MergedStaticString<T> &s) {
    return s.to_string();
}
}
