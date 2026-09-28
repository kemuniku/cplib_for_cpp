#pragma once
#include <cplib/str/private/hash_string_base.hpp>
#include <concepts>
namespace cplib {
using HashString=BasicHashString<true>;using RollingHash=BasicRollingHash<true>;
template<std::integral T> requires (!std::is_same_v<T,char>) HashString tohash(T s){return HashString::character(Int(s));}inline HashString tohash(char s){return HashString::character(static_cast<unsigned char>(s));}template<class T> HashString tohash(const std::vector<T>& s){return HashString::sequence(s);}inline HashString tohash(std::string_view s){return HashString::sequence(s);}
inline HashString get_emptystring_hash(){return {};}
inline RollingHash initRollingHash(std::string_view s){return RollingHash(s);}inline RollingHash initRollingHash(const std::vector<char>& s){return RollingHash(std::string_view(s.data(),s.size()));}
using cplib::isPalindrome;using cplib::reversed;using cplib::LCP;using cplib::cmp;using cplib::toHashString;
}
