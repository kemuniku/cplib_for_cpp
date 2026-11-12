#pragma once
#include <cplib/tmpl/private/template_support.hpp>

namespace cplib {
using namespace detail::template_support;
inline constexpr Int INF = INF64;

// 整数を一つ読み取る。旧版のscanf末尾空白消費を維持する。読み取る文字数に比例。
inline Int ii() {
    long long value = 0;
    [[maybe_unused]] int count = std::scanf("%lld\n", &value);
    return Int(value);
}

// 整数をn個読み取る。計算量は入力した文字数とnに比例する。
inline std::vector<Int> lii(Int n) {
    assert(n >= 0);
    std::vector<Int> out(n);
    for (auto &x : out)
        x = ii();
    return out;
}

// 空白を読み飛ばさず次の空白まで文字列を読む。計算量は読み取った文字数に比例。
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

// 各要素を個別の文字列に変換して区切り文字を挟んで連結する。
// 要素数と出力長に比例する時間に、各要素の文字列変換の費用を加える。
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

// 真ならYes、偽ならNoを出力する。計算量 O(1) と出力費用。
inline void yes(bool b = true) {
    detail::template_support::print(b ? "Yes" : "No");
}

// 真ならNo、偽ならYesを出力する。計算量 O(1) と出力費用。
inline void no(bool b = true) {
    yes(!b);
}

// 真ならTakahashi、偽ならAokiを出力する。計算量 O(1) と出力費用。
inline void takahashi(bool b = true) {
    detail::template_support::print(b ? "Takahashi" : "Aoki");
}

// 真ならAoki、偽ならTakahashiを出力する。計算量 O(1) と出力費用。
inline void aoki(bool b = true) {
    takahashi(!b);
}
}
