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
    // 指定したストリームを入力元にする。計算量 O(1)。
    explicit Input(std::istream &stream = std::cin) : stream_(&stream) {
    }

    // 次のトークン、またはトークンの次の1文字を読む。終端では空文字列を返す。
    // 計算量は今回走査した入力文字数と返す文字数に比例する。
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

    // 指定型の値を1つ読む。計算量は readNext と文字列の型変換にかかる時間の和。
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

// 指定した型の値を順に読み、1個なら値、多数なら組で返す。
// 計算量は各値の読み込みと型変換の合計。
template <class... T> auto input() {
    static_assert(sizeof...(T) > 0);
    if constexpr (sizeof...(T) == 1)
        return reader.template input<std::tuple_element_t<0, std::tuple<T...>>>();
    else
        return std::tuple<T...>{reader.template input<T>()...};
}

// 指定型の値または組を n 個読み、配列で返す。
// 計算量は O(n) と各値の読み込み・型変換の合計。
template <class... T> auto input(Int n) {
    assert(n >= 0);
    using V = decltype(input<T...>());
    std::vector<V> out(n);
    for (auto &x : out)
        x = input<T...>();
    return out;
}

// 指定された各次元の配列に値を読み込む。
// 計算量は作成する配列・要素数と、各値の読み込み・型変換の合計。
template <class T, class... D> auto input(Int n, Int m, D... dimensions) {
    assert(n >= 0);
    using V = decltype(input<T>(m, dimensions...));
    std::vector<V> out(n);
    for (auto &x : out)
        x = input<T>(m, dimensions...);
    return out;
}

// 値やコンテナを出力用文字列に変換する。
// 計算量は要素ごとの文字列化と連結に要する時間の合計。
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

// 優先度付きキューを取り出し順に文字列化する。比較・要素のコピーを O(1) とすると、
// n要素なら O(n log(n+1)) と各要素の文字列化・連結の時間。
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

// 出力設定に従い各値を文字列化して書き出す。計算量は各値の変換と出力に要する時間の合計。
template <class... T> void print(const PrintOptions &options, const T &...values) {
    std::array<std::string, sizeof...(T)> formatted{fmtprint(values)...};
    detail::template_support::write(options, formatted);
}

// 各値を既定の区切りで出力する。計算量は設定付き print と同じ。
template <class... T> void print(const T &...values) {
    cplib::print(PrintOptions{}, values...);
}

// aの非負整数n乗を法mで計算する。計算量 O(log(n+1))。
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

// 真なら Yes、偽なら No を出力する。計算量 O(1)。
inline void Yes(bool b = true) {
    cplib::print(b ? "Yes" : "No");
}

// 真なら No、偽なら Yes を出力する。計算量 O(1)。
inline void No(bool b = true) {
    Yes(!b);
}

// 真なら YES、偽なら NO を出力する。計算量 O(1)。
inline void YES_upper(bool b = true) {
    cplib::print(b ? "YES" : "NO");
}

// 真なら NO、偽なら YES を出力する。計算量 O(1)。
inline void NO_upper(bool b = true) {
    YES_upper(!b);
}

// 空のハッシュ集合を作る。計算量 O(1)。
template <class T> auto initHashSet() {
    return std::unordered_set<T>{};
}

// 関数を実行してから正常終了する。計算量は関数と終了処理の実行時間の合計。
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
