#pragma once
#include <cplib/str/private/merged_string_base.hpp>

namespace cplib {
template <class T, std::size_t N> struct FixedLengthMergedStaticString {
    using value_type = T;
    std::shared_ptr<StaticStringBase<T>> base;
    std::array<std::int32_t, N> L{}, R{};

    void setRange(std::size_t i, std::shared_ptr<StaticStringBase<T>> b, std::int32_t l,
                  std::int32_t r) {
        assert(i < N && l <= r);
        if (i == 0)
            base = b;
        else
            assert(base == b);
        L[i] = l;
        R[i] = r;
    }

    Int len() const {
        return detail::merged_length(*this);
    }

    T operator[](Int i) const {
        assert(0 <= i && i < len());
        return detail::merged_at(*this, i);
    }

    T operator[](BackwardsIndex i) const {
        return (*this)[len() - i.value];
    }

    template <class A, class B>
    FixedLengthMergedStaticString operator[](ClosedSlice<A, B> slice) const {
        Int a = resolve_index(len(), slice.a), b = resolve_index(len(), slice.b);
        assert(0 <= a && a <= b + 1 && b < len());
        FixedLengthMergedStaticString out;
        out.base = base;
        Int previous = 0;
        for (std::size_t i = 0; i < N; ++i) {
            Int current = previous + R[i] - L[i], left = std::max(previous, a),
                right = std::min(current, b + 1);
            if (left < right) {
                out.L[i] = L[i] + left - previous;
                out.R[i] = L[i] + right - previous;
            } else
                out.L[i] = out.R[i] = L[i];
            previous = current;
        }
        return out;
    }

    std::string to_string() const {
        return detail::merged_string(*this);
    }
};

template <class T, std::size_t N>
auto initFixedLengthMergedStaticString(const std::array<StaticString<T>, N> &s) {
    FixedLengthMergedStaticString<T, N> out;
    for (std::size_t i = 0; i < N; ++i)
        out.setRange(i, s[i].base, s[i].l, s[i].r);
    return out;
}

template <class T, std::size_t N>
auto initFixedLengthMergedStaticString(const StaticString<T> &s,
                                       const std::array<std::pair<Int, Int>, N> &ranges) {
    FixedLengthMergedStaticString<T, N> out;
    for (std::size_t i = 0; i < N; ++i) {
        auto [l, r] = ranges[i];
        assert(0 <= l && l <= r && r <= s.len());
        out.setRange(i, s.base, s.l + l, s.l + r);
    }
    return out;
}

template <class T> auto operator&(const StaticString<T> &s, const StaticString<T> &t) {
    return initFixedLengthMergedStaticString(std::array{s, t});
}

template <class T, std::size_t N>
auto operator&(const FixedLengthMergedStaticString<T, N> &s, const StaticString<T> &t) {
    FixedLengthMergedStaticString<T, N + 1> out;
    if constexpr (N > 0) {
        assert(s.base == t.base);
        out.base = s.base;
        std::copy(s.L.begin(), s.L.end(), out.L.begin());
        std::copy(s.R.begin(), s.R.end(), out.R.begin());
    } else
        out.base = t.base;
    out.L[N] = t.l;
    out.R[N] = t.r;
    return out;
}

template <class T, std::size_t N, std::size_t M>
auto operator&(const FixedLengthMergedStaticString<T, N> &s,
               const FixedLengthMergedStaticString<T, M> &t) {
    FixedLengthMergedStaticString<T, N + M> out;
    if constexpr (N > 0 && M > 0)
        assert(s.base == t.base);
    if constexpr (N > 0) {
        out.base = s.base;
        std::copy(s.L.begin(), s.L.end(), out.L.begin());
        std::copy(s.R.begin(), s.R.end(), out.R.begin());
    } else if constexpr (M > 0)
        out.base = t.base;
    std::copy(t.L.begin(), t.L.end(), out.L.begin() + N);
    std::copy(t.R.begin(), t.R.end(), out.R.begin() + N);
    return out;
}

template <class T, std::size_t N, std::size_t M>
Int lcp(const FixedLengthMergedStaticString<T, N> &s,
        const FixedLengthMergedStaticString<T, M> &t) {
    return detail::merged_lcp(s, t);
}

template <class T, std::size_t N, std::size_t M>
int cmp(const FixedLengthMergedStaticString<T, N> &s,
        const FixedLengthMergedStaticString<T, M> &t) {
    return detail::merged_cmp(s, t);
}

template <class T, std::size_t N, std::size_t M>
bool operator==(const FixedLengthMergedStaticString<T, N> &s,
                const FixedLengthMergedStaticString<T, M> &t) {
    return s.len() == t.len() && lcp(s, t) == s.len();
}

template <class T, std::size_t N, std::size_t M>
auto operator<=>(const FixedLengthMergedStaticString<T, N> &s,
                 const FixedLengthMergedStaticString<T, M> &t) {
    return cmp(s, t) <=> 0;
}

template <class T, std::size_t N> Int len(const FixedLengthMergedStaticString<T, N> &s) {
    return s.len();
}

template <class T, std::size_t N>
std::string to_string(const FixedLengthMergedStaticString<T, N> &s) {
    return s.to_string();
}
}
