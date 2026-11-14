#pragma once
#include <cplib/str/private/hash_string_base.hpp>
#include <concepts>

namespace cplib::reversible_hash {
using HashString = BasicHashString<true>;
using RollingHash = BasicRollingHash<true>;

// 1要素のハッシュを返す。計算量 O(1)。
template <std::integral T>
    requires(!std::is_same_v<T, char>)
HashString tohash(T s) {
    return HashString::character(Int(s));
}

// 1文字のハッシュを返す。計算量 O(1)。
inline HashString tohash(char s) {
    return HashString::character(static_cast<unsigned char>(s));
}

// 長さNの列のハッシュを返す。計算量 O(N+1)。
template <class T> HashString tohash(const std::vector<T> &s) {
    return HashString::sequence(s);
}

// 長さNの文字列のハッシュを返す。計算量 O(N+1)。
inline HashString tohash(std::string_view s) {
    return HashString::sequence(s);
}

// 空文字列のハッシュを返す。計算量 O(1)。
inline HashString get_emptystring_hash() {
    return {};
}

// 長さNの文字列から累積ハッシュを作る。計算量 O(N+1)。
inline RollingHash initRollingHash(std::string_view s) {
    return RollingHash(s);
}

// N文字から累積ハッシュを作る。計算量 O(N+1)。
inline RollingHash initRollingHash(const std::vector<char> &s) {
    return RollingHash(std::string_view(s.data(), s.size()));
}

using cplib::isPalindrome;
using cplib::reversed;
using cplib::LCP;
using cplib::cmp;
using cplib::toHashString;
}
