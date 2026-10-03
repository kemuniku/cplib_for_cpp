#pragma once
#include <cplib/common.hpp>

namespace cplib {
// 二次元配列の90度単位の回転・転置。O(HW)。
template <class T> auto rotated(const std::vector<std::vector<T>> &a, Int turns = 1) {
    if (a.empty())
        return a;
    Int h = a.size(), w = a[0].size();
    turns = (turns % 4 + 4) % 4;
    if (!turns)
        return a;
    std::vector<std::vector<T>> out;
    if (turns == 1 || turns == 3) {
        out.assign(w, std::vector<T>(h));
        for (Int i = 0; i < w; ++i)
            for (Int j = 0; j < h; ++j)
                out[i][j] = turns == 1 ? a[h - 1 - j][i] : a[j][w - 1 - i];
    } else {
        out.assign(h, std::vector<T>(w));
        for (Int i = 0; i < h; ++i)
            for (Int j = 0; j < w; ++j)
                out[i][j] = a[h - 1 - i][w - 1 - j];
    }
    return out;
}

template <class T> auto transposed(const std::vector<std::vector<T>> &a) {
    if (a.empty())
        return a;
    Int h = a.size(), w = a[0].size();
    std::vector<std::vector<T>> out(w, std::vector<T>(h));
    for (Int i = 0; i < w; ++i)
        for (Int j = 0; j < h; ++j)
            out[i][j] = a[j][i];
    return out;
}

template <class T> auto trimmed(std::vector<std::vector<T>> a, const T &zero) {
    for (int turn = 0; turn < 4; ++turn) {
        a = rotated(a);
        while (!a.empty() && std::all_of(a.back().begin(), a.back().end(),
                                         [&](const T &value) { return value == zero; }))
            a.pop_back();
    }
    return a;
}

namespace detail {
inline auto matops_chars(const std::vector<std::string> &a) {
    std::vector<std::vector<char>> out;
    for (const auto &row : a)
        out.emplace_back(row.begin(), row.end());
    return out;
}

inline auto matops_strings(const std::vector<std::vector<char>> &a) {
    std::vector<std::string> out;
    for (const auto &row : a)
        out.emplace_back(row.begin(), row.end());
    return out;
}
}

inline auto rotated(const std::vector<std::string> &a, Int turns = 1) {
    return detail::matops_strings(rotated(detail::matops_chars(a), turns));
}

inline auto transposed(const std::vector<std::string> &a) {
    return detail::matops_strings(transposed(detail::matops_chars(a)));
}

inline auto trimmed(const std::vector<std::string> &a, char zero = '0') {
    return detail::matops_strings(trimmed(detail::matops_chars(a), zero));
}

template <class T> void rotate(std::vector<std::vector<T>> &a, Int turns = 1) {
    a = rotated(a, turns);
}

inline void rotate(std::vector<std::string> &a, Int turns = 1) {
    a = rotated(a, turns);
}

template <class T> void transpose(std::vector<std::vector<T>> &a) {
    a = transposed(a);
}

inline void transpose(std::vector<std::string> &a) {
    a = transposed(a);
}

template <class T> void trim(std::vector<std::vector<T>> &a, const T &zero) {
    a = trimmed(a, zero);
}

inline void trim(std::vector<std::string> &a, char zero = '0') {
    a = trimmed(a, zero);
}
}
