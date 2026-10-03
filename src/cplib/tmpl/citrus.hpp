#pragma once
#include <cplib/tmpl/private/template_support.hpp>
#include <cplib/math/isqrt.hpp>
#include <iostream>
#include <cctype>

namespace cplib {
using detail::template_support::PrintOptions;
using detail::template_support::chmin;
using detail::template_support::chmax;
using detail::template_support::min_assign;
using detail::template_support::max_assign;
using detail::template_support::floor_mod;
using detail::template_support::floor_div;
using detail::template_support::floor_mod_assign;
using detail::template_support::floor_div_assign;
using detail::template_support::bit;
using detail::template_support::set_bit;
using detail::template_support::at;
using cplib::isqrt;
inline constexpr Int MODINT998244353 = 998244353, MODINT1000000007 = 1000000007, INFL = INF64;
inline constexpr std::array<std::pair<Int, Int>, 4> DXY{{{0, -1}, {0, 1}, {-1, 0}, {1, 0}}};
inline constexpr std::array<std::pair<Int, Int>, 8> DDXY{
    {{1, -1}, {1, 0}, {1, 1}, {0, -1}, {0, 1}, {-1, -1}, {-1, 0}, {-1, 1}}};

// 行を分割する旧入力イテレータ。文字入力中は次の呼び出しも同じトークンの残りを読む。
class Input {
    std::istream *stream_;
    std::vector<std::string> tokens_;
    std::size_t token_ = 0, character_ = 0;
    bool characters_ = false;

public:
    explicit Input(std::istream &stream = std::cin) : stream_(&stream) {
    }

    std::string readNext(bool getsChar = false) {
        while (true) {
            if (token_ < tokens_.size()) {
                auto &word = tokens_[token_];
                if (character_ == 0)
                    characters_ = getsChar;
                if (characters_) {
                    std::string out(1, word[character_++]);
                    if (character_ == word.size()) {
                        character_ = 0;
                        ++token_;
                    }
                    return out;
                }
                ++token_;
                return word;
            }
            std::string line;
            if (!std::getline(*stream_, line))
                return "";
            tokens_.clear();
            token_ = character_ = 0;
            for (std::size_t i = 0; i < line.size();) {
                while (i < line.size() && std::isspace(static_cast<unsigned char>(line[i])))
                    ++i;
                std::size_t first = i;
                while (i < line.size() && !std::isspace(static_cast<unsigned char>(line[i])))
                    ++i;
                if (first < i)
                    tokens_.push_back(line.substr(first, i - first));
            }
        }
    }

    template <class T> T input() {
        if constexpr (std::is_same_v<T, std::string>)
            return readNext();
        else if constexpr (std::is_same_v<T, char>) {
            auto s = readNext(true);
            if (s.empty())
                throw std::out_of_range("入力の終端");
            return s[0];
        } else if constexpr (std::integral<T>) {
            auto s = readNext();
            std::size_t end;
            auto value = std::stoll(s, &end);
            if (end != s.size())
                throw std::invalid_argument("整数入力が不正です");
            return T(value);
        } else if constexpr (std::floating_point<T>) {
            auto s = readNext();
            std::size_t end;
            auto value = std::stod(s, &end);
            if (end != s.size())
                throw std::invalid_argument("実数入力が不正です");
            return T(value);
        } else
            static_assert(sizeof(T) == 0, "unsupported citrus input type");
    }
};

inline Input reader;

template <class... T> auto input() {
    static_assert(sizeof...(T) > 0);
    if constexpr (sizeof...(T) == 1)
        return reader.template input<std::tuple_element_t<0, std::tuple<T...>>>();
    else
        return std::tuple<T...>{reader.template input<T>()...};
}

template <class... T> auto input(Int n) {
    assert(n >= 0);
    using V = decltype(input<T...>());
    std::vector<V> out(n);
    for (auto &x : out)
        x = input<T...>();
    return out;
}

template <class T, class... D> auto input(Int n, Int m, D... dimensions) {
    assert(n >= 0);
    using V = decltype(input<T>(m, dimensions...));
    std::vector<V> out(n);
    for (auto &x : out)
        x = input<T>(m, dimensions...);
    return out;
}

template <class T> std::string fmtprint(const T &x) {
    if constexpr (std::floating_point<T>) {
        std::ostringstream os;
        os << std::fixed << std::setprecision(16) << x;
        return os.str();
    } else if constexpr (requires {
                             typename T::mapped_type;
                             x.begin();
                         }) {
        std::string out;
        for (const auto &[key, value] : x) {
            if (!out.empty())
                out += ' ';
            out += detail::template_support::stringValue(key) + ": " +
                   detail::template_support::stringValue(value);
        }
        return out;
    } else if constexpr (!std::is_convertible_v<T, std::string_view> && requires {
                             x.begin();
                             x.end();
                         }) {
        std::string out;
        bool first = true;
        for (const auto &value : x) {
            if (!first)
                out += ' ';
            first = false;
            out += detail::template_support::stringValue(value);
        }
        return out;
    } else if constexpr (!std::is_convertible_v<T, std::string_view> &&
                         requires { std::tuple_size<T>::value; }) {
        std::string out;
        std::size_t i = 0;
        std::apply(
            [&](const auto &...value) {
                ((out += (i++ ? " " : "") + detail::template_support::stringValue(value)), ...);
            },
            x);
        return out;
    } else
        return detail::template_support::stringValue(x);
}

template <class T, class C, class Compare>
std::string fmtprint(std::priority_queue<T, C, Compare> q) {
    std::string out;
    while (!q.empty()) {
        out += detail::template_support::stringValue(q.top());
        q.pop();
        if (!q.empty())
            out += ' ';
    }
    return out;
}

template <class... T> void print(const PrintOptions &options, const T &...values) {
    std::array<std::string, sizeof...(T)> formatted{fmtprint(values)...};
    detail::template_support::write(options, formatted);
}

template <class... T> void print(const T &...values) {
    cplib::print(PrintOptions{}, values...);
}

inline Int pow(Int a, Int n, Int m = INF64) {
    Int result = 1;
    while (n > 0) {
        if (n % 2)
            result = (result * a) % m;
        if (n > 1)
            a = (a * a) % m;
        n >>= 1;
    }
    return result;
}

inline void Yes(bool b = true) {
    cplib::print(b ? "Yes" : "No");
}

inline void No(bool b = true) {
    Yes(!b);
}

inline void YES_upper(bool b = true) {
    cplib::print(b ? "YES" : "NO");
}

inline void NO_upper(bool b = true) {
    YES_upper(!b);
}

template <class T> auto initHashSet() {
    return std::unordered_set<T>{};
}

template <class F> [[noreturn]] void exit(F &&statement) {
    std::forward<F>(statement)();
    std::exit(0);
}
}
#if defined(CPLIB_LOCAL_DEBUG) && CPLIB_LOCAL_DEBUG
#define CPLIB_CITRUS_DEBUG(...)                                                                    \
    ::cplib::detail::template_support::debugPrintValues(#__VA_ARGS__, __VA_ARGS__)
#else
#define CPLIB_CITRUS_DEBUG(...) ((void)0)
#endif
