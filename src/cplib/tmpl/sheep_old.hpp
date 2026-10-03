#pragma once
#include <cplib/tmpl/private/template_support.hpp>

namespace cplib {
using namespace detail::template_support;
inline constexpr Int INF = INF64;

// 旧版のscanf末尾空白消費と、空白を読み飛ばさない文字列入力を維持する。
inline Int ii() {
    long long value = 0;
    [[maybe_unused]] int count = std::scanf("%lld\n", &value);
    return Int(value);
}

inline std::vector<Int> lii(Int n) {
    assert(n >= 0);
    std::vector<Int> out(n);
    for (auto &x : out)
        x = ii();
    return out;
}

inline std::string si() {
    std::string out;
    while (true) {
        int c = getchar_unlocked();
        if (c == ' ' || c == '\n' || c == -1 || c == 255)
            break;
        out += char(c);
    }
    return out;
}

// 旧joinは全要素を個別の文字列に変換してから連結する。
template <class R> std::string join(const R &values, std::string_view separator = "") {
    std::vector<std::string> parts;
    std::size_t length = 0;
    for (const auto &value : values) {
        parts.push_back(stringValue(value));
        length += parts.back().size();
    }
    if (!parts.empty())
        length += (parts.size() - 1) * separator.size();
    std::string out;
    out.reserve(length);
    for (std::size_t i = 0; i < parts.size(); ++i) {
        if (i)
            out += separator;
        out += parts[i];
    }
    return out;
}

inline void yes(bool b = true) {
    detail::template_support::print(b ? "Yes" : "No");
}

inline void no(bool b = true) {
    yes(!b);
}

inline void takahashi(bool b = true) {
    detail::template_support::print(b ? "Takahashi" : "Aoki");
}

inline void aoki(bool b = true) {
    takahashi(!b);
}
}
