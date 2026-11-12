#pragma once
#include <cplib/tmpl/private/template_support.hpp>
#include <cplib/common.hpp>
#include <cplib/tmpl/private/fastio_native.hpp>
#include <span>
#include <sstream>
#include <array>

namespace cplib {
struct PrintOptions {
    FILE *f = stdout;
    std::string sepc = " ", endc = "\n";
    bool flush = false;
};

// 値を文字列に変換する。浮動小数点や独自型はストリーム出力を利用する。計算量は値の変換と出力長に比例。
template <class T> std::string stringify(const T &value) {
    return detail::template_support::stringValue(value);
}

// 型と符号をコンパイル時に選ぶ。uint32単独/配列は元のSIMD専用経路を使用する。
// 指定型の値を一つ読む。計算量は読み取る文字数に比例。
template <class T = Int> __attribute__((target("avx2"))) T input() {
    if constexpr (std::is_same_v<T, std::string>) {
        std::size_t length;
        auto p = detail::fastio_native::cplib_fio_read_token(&length);
        return std::string(p, length);
    } else {
        static_assert(std::integral<T> && !std::is_same_v<T, bool> && !std::is_same_v<T, char>);
        if constexpr (std::is_signed_v<T>)
            return T(detail::fastio_native::cplib_fio_read_int());
        else if constexpr (std::is_same_v<T, std::uint32_t>)
            return detail::fastio_native::cplib_fio_read_u32();
        else
            return T(detail::fastio_native::cplib_fio_read_uint());
    }
}

// 指定型の値をn個読む。計算量はnと読み取る文字数に比例。
template <class T = Int> __attribute__((target("avx2"))) std::vector<T> input(Int n) {
    assert(n >= 0);
    std::vector<T> out(n);
    if (n) {
        if constexpr (std::is_same_v<T, std::string>) {
            for (auto &x : out)
                x = input<T>();
        } else {
            static_assert(std::integral<T> && !std::is_same_v<T, bool> && !std::is_same_v<T, char>);
            if constexpr (std::is_signed_v<T>) {
                if constexpr (sizeof(T) == 8)
                    detail::fastio_native::cplib_fio_read_i64_array(out.data(), n);
                else if constexpr (sizeof(T) == 4)
                    detail::fastio_native::cplib_fio_read_i32_array(out.data(), n);
                else if constexpr (sizeof(T) == 2)
                    detail::fastio_native::cplib_fio_read_i16_array(out.data(), n);
                else
                    detail::fastio_native::cplib_fio_read_i8_array(out.data(), n);
            } else {
                if constexpr (sizeof(T) == 8)
                    detail::fastio_native::cplib_fio_read_u64_array(out.data(), n);
                else if constexpr (sizeof(T) == 4)
                    detail::fastio_native::cplib_fio_read_u32_array(
                        reinterpret_cast<std::uint32_t *>(out.data()), n);
                else if constexpr (sizeof(T) == 2)
                    detail::fastio_native::cplib_fio_read_u16_array(out.data(), n);
                else
                    detail::fastio_native::cplib_fio_read_u8_array(out.data(), n);
            }
        }
    }
    return out;
}

// 整数を一つ読む。計算量は読み取る文字数に比例。
inline Int ii() {
    return input<Int>();
}

// 整数をn個読む。計算量はnと読み取る文字数に比例。
inline std::vector<Int> lii(Int n) {
    return input<Int>(n);
}

// 文字列を一つ読む。計算量は読み取る文字数に比例。
inline std::string si() {
    return input<std::string>();
}

// 入力から1バイトを読み取る。償却 O(1) と入力費用。
inline int fastioGetChar() {
    return detail::fastio_native::cplib_fio_get_char();
}

// 文字列群を指定の区切りで出力する。文字列数と出力文字数に比例する時間と出力費用。
inline void print_internal(const PrintOptions &prop, std::span<const std::string> args) {
    for (std::size_t i = 0; i < args.size(); ++i) {
        std::fwrite(args[i].data(), 1, args[i].size(), prop.f);
        const auto &sep = i + 1 == args.size() ? prop.endc : prop.sepc;
        std::fwrite(sep.data(), 1, sep.size(), prop.f);
    }
    if (prop.flush)
        std::fflush(prop.f);
}

// 設定に従って値を変換・出力する。各値の変換と出力文字数に比例する時間と出力費用。
template <class... T> void print(const PrintOptions &prop, const T &...args) {
    std::array<std::string, sizeof...(T)> values{stringify(args)...};
    print_internal(prop, values);
}

// 値を区切り文字を挟んで連結する。要素数と出力長に比例する時間に各値の変換費用を加える。
template <class T> std::string join(std::span<const T> values, std::string_view sep = "") {
    if (values.empty())
        return {};
    if constexpr (std::integral<T> && (sizeof(T) == 4 || sizeof(T) == 8)) {
        constexpr std::size_t digits = sizeof(T) == 8 ? 20 : std::is_signed_v<T> ? 11 : 10;
        std::string out(values.size() * digits + (values.size() - 1) * sep.size(), '\0');
        std::size_t written;
        if constexpr (std::is_signed_v<T>) {
            if constexpr (sizeof(T) == 8)
                written = detail::fastio_native::cplib_fio_join_i64(
                    values.data(), values.size(), out.data(), sep.data(), sep.size());
            else
                written = detail::fastio_native::cplib_fio_join_i32(
                    values.data(), values.size(), out.data(), sep.data(), sep.size());
        } else {
            if constexpr (sizeof(T) == 8)
                written = detail::fastio_native::cplib_fio_join_u64(
                    values.data(), values.size(), out.data(), sep.data(), sep.size());
            else
                written = detail::fastio_native::cplib_fio_join_u32(
                    values.data(), values.size(), out.data(), sep.data(), sep.size());
        }
        out.resize(written);
        return out;
    } else if constexpr (requires(T x) {
                             T::umod();
                             x.val();
                         }) {
        std::vector<std::uint32_t> canonical(values.size());
        for (std::size_t i = 0; i < values.size(); ++i)
            canonical[i] = values[i].val();
        return join<std::uint32_t>(canonical, sep);
    } else {
        std::string out;
        for (std::size_t i = 0; i < values.size(); ++i) {
            if (i)
                out += sep;
            out += stringify(values[i]);
        }
        return out;
    }
}

// vectorの値を区切り文字を挟んで連結する。計算量はspan版joinと同じ。
template <class T, class A>
std::string join(const std::vector<T, A> &values, std::string_view sep = "") {
    if constexpr (std::is_same_v<T, bool>) {
        std::string out;
        for (std::size_t i = 0; i < values.size(); ++i) {
            if (i)
                out += sep;
            out += values[i] ? "true" : "false";
        }
        return out;
    } else
        return join<T>(std::span<const T>(values.data(), values.size()), sep);
}

// arrayの値を区切り文字を挟んで連結する。計算量はspan版joinと同じ。
template <class T, std::size_t N>
std::string join(const std::array<T, N> &values, std::string_view sep = "") {
    return join<T>(std::span<const T>(values), sep);
}

// 配列出力は固定長バッファを使用する。巨大な区切り文字も分割して書き出す。
// 列の要素を指定区切りで出力する。各値の変換と出力文字数に比例する時間と出力費用。
template <class T>
void print_array(std::span<const T> values, std::string_view sep = " ", FILE *file = stdout) {
    if (values.empty()) {
        std::fwrite("\n", 1, 1, file);
        return;
    }
    if constexpr (std::integral<T> && (sizeof(T) == 4 || sizeof(T) == 8)) {
        if constexpr (std::is_signed_v<T>) {
            if constexpr (sizeof(T) == 8)
                detail::fastio_native::cplib_fio_print_array_i64(file, values.data(), values.size(),
                                                                 sep.data(), sep.size());
            else
                detail::fastio_native::cplib_fio_print_array_i32(file, values.data(), values.size(),
                                                                 sep.data(), sep.size());
        } else {
            if constexpr (sizeof(T) == 8)
                detail::fastio_native::cplib_fio_print_array_u64(file, values.data(), values.size(),
                                                                 sep.data(), sep.size());
            else
                detail::fastio_native::cplib_fio_print_array_u32(file, values.data(), values.size(),
                                                                 sep.data(), sep.size());
        }
    } else if constexpr (requires(T x) {
                             T::umod();
                             x.val();
                         }) {
        std::vector<std::uint32_t> canonical(values.size());
        for (std::size_t i = 0; i < values.size(); ++i)
            canonical[i] = values[i].val();
        print_array<std::uint32_t>(canonical, sep, file);
    } else {
        auto text = join<T>(values, sep);
        std::fwrite(text.data(), 1, text.size(), file);
        std::fwrite("\n", 1, 1, file);
    }
}

// vectorの要素を指定区切りで出力する。計算量はspan版print_arrayと同じ。
template <class T, class A>
void print_array(const std::vector<T, A> &values, std::string_view sep = " ", FILE *file = stdout) {
    if constexpr (std::is_same_v<T, bool>) {
        auto text = join(values, sep);
        std::fwrite(text.data(), 1, text.size(), file);
        std::fwrite("\n", 1, 1, file);
    } else
        print_array<T>(std::span<const T>(values.data(), values.size()), sep, file);
}

// arrayの要素を指定区切りで出力する。計算量はspan版print_arrayと同じ。
template <class T, std::size_t N>
void print_array(const std::array<T, N> &values, std::string_view sep = " ", FILE *file = stdout) {
    print_array<T>(std::span<const T>(values), sep, file);
}

// 1値の整数出力はstdio自身のバッファへ直接書き込む経路を維持する。
// 値を一つ出力する。値の変換と出力に要する時間に比例。
template <class T> void print_one(const T &value, FILE *file = stdout) {
    if constexpr (std::integral<T> && (sizeof(T) == 4 || sizeof(T) == 8)) {
        if constexpr (std::is_signed_v<T>) {
            if constexpr (sizeof(T) == 8)
                detail::fastio_native::cplib_fio_print_one_i64(file, value);
            else
                detail::fastio_native::cplib_fio_print_one_i32(file, value);
        } else {
            if constexpr (sizeof(T) == 8)
                detail::fastio_native::cplib_fio_print_one_u64(file, value);
            else
                detail::fastio_native::cplib_fio_print_one_u32(file, value);
        }
    } else if constexpr (requires {
                             T::umod();
                             value.val();
                         })
        detail::fastio_native::cplib_fio_print_one_u32(file, value.val());
    else
        print(PrintOptions{file}, value);
}

struct Separator {
    std::string value;
};

// 出力用の区切り文字列を作る。文字列の移動に加え O(1)。
inline Separator sep(std::string value) {
    return {std::move(value)};
}

// 出力用の区切り文字を作る。計算量 O(1)。
inline Separator sep(char c) {
    return {std::string(1, c)};
}

template <class R> struct Splat {
    const R &values;
};

// 列を出力するための参照を作る。計算量 O(1)。
template <class R> Splat<R> splat(const R &values) {
    return {values};
}

// 列の要素を空白区切りで出力する。計算量はprint_arrayと同じ。
template <class R> void print(Splat<R> values) {
    print_array(values.values);
}

// 列の要素を指定区切りで出力する。計算量はprint_arrayと同じ。
template <class R> void print(Separator separator, Splat<R> values) {
    print_array(values.values, separator.value);
}

// 値を一つ出力する。計算量はprint_oneと同じ。
template <class T>
    requires(!std::is_same_v<T, Separator> && !std::is_same_v<T, PrintOptions>)
void print(const T &value) {
    print_one(value);
}

// 値を指定区切りで出力する。値の変換と出力に要する時間に比例。
template <class... T> void print(Separator separator, const T &...values) {
    if constexpr (sizeof...(T) > 0 && (std::integral<T> && ...)) {
        using V = std::common_type_t<T...>;
        if constexpr ((std::is_same_v<V, T> && ...)) {
            std::array<V, sizeof...(T)> a{values...};
            print_array(a, separator.value);
        } else
            print(PrintOptions{stdout, separator.value, "\n", false}, values...);
    } else
        print(PrintOptions{stdout, separator.value, "\n", false}, values...);
}

// 引数なしでは何も出力しない。計算量 O(1)。
inline void print() {
} // 元の0引数版は改行も出力しない。

// 複数の値を空白区切りで出力する。値の変換と出力に要する時間に比例。
template <class T, class U, class... V>
    requires(!std::is_same_v<T, Separator> && !std::is_same_v<T, PrintOptions>)
void print(const T &a, const U &b, const V &...rest) {
    print(sep(" "), a, b, rest...);
}

struct InputSource {
    // 指定型の値を一つ読む。計算量は読み取る文字数に比例。
    template <class T> T input() {
        return cplib::input<T>();
    }
};
}
