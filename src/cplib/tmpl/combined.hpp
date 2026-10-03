#pragma once
#include <cplib/tmpl/private/template_support.hpp>
#include <cstdio>
#include <cstdint>
#include <cstring>
#include <sys/mman.h>
#include <sys/stat.h>
#include <cstddef>
#include <type_traits>
// 展開済みファイル固有の旧版: 64 MiB 入力と独立した4桁テーブル。
#if defined(__clang__)
#pragma clang attribute push(__attribute__((target("avx2"))), apply_to = function)
#elif defined(__GNUC__)
#pragma GCC push_options
#pragma GCC target("avx2")
#pragma GCC optimize("O3", "unroll-loops")
#endif
namespace cplib {

namespace detail::combined_input_native {
constexpr std::size_t buffer_size = 1U << 26;
inline char buffer[buffer_size];
inline std::size_t cursor = 0;
inline std::size_t length = 0;
inline const char *mapped = nullptr;
inline bool initialized = false;

inline void initialize() {
    if (initialized)
        return;
    initialized = true;

    struct stat st;
    const int fd = fileno(stdin);
    if (fstat(fd, &st) == 0 && S_ISREG(st.st_mode) && st.st_size > 0) {
        void *p =
            mmap(nullptr, static_cast<std::size_t>(st.st_size), PROT_READ, MAP_PRIVATE, fd, 0);
        if (p != MAP_FAILED) {
            mapped = static_cast<const char *>(p);
            length = static_cast<std::size_t>(st.st_size);
            madvise(const_cast<char *>(mapped), length, MADV_SEQUENTIAL);
        }
    }
}

inline int get_char() {
    initialize();
    if (mapped != nullptr) {
        if (cursor == length)
            return -1;
        return static_cast<unsigned char>(mapped[cursor++]);
    }

    if (cursor == length) {
        length = fread_unlocked(buffer, 1, buffer_size, stdin);
        cursor = 0;
        if (length == 0)
            return -1;
    }
    return static_cast<unsigned char>(buffer[cursor++]);
}

inline bool refill() {
    length = fread_unlocked(buffer, 1, buffer_size, stdin);
    cursor = 0;
    return length != 0;
}

inline bool has_eight_digits(const char *source) {
    std::uint64_t bytes;
    std::memcpy(&bytes, source, sizeof(bytes));
    constexpr std::uint64_t high_nibbles = 0xf0f0f0f0f0f0f0f0ULL;
    return (bytes & high_nibbles) == 0x3030303030303030ULL &&
           ((bytes + 0x0606060606060606ULL) & high_nibbles) == 0x3030303030303030ULL;
}

inline unsigned parse_eight_digits(const char *source) {
#if defined(__BYTE_ORDER__) && __BYTE_ORDER__ == __ORDER_LITTLE_ENDIAN__
    std::uint64_t digits;
    std::memcpy(&digits, source, sizeof(digits));
    digits -= 0x3030303030303030ULL;
    digits = (digits * 10 + (digits >> 8)) & 0x00ff00ff00ff00ffULL;
    digits = (digits * 100 + (digits >> 16)) & 0x0000ffff0000ffffULL;
    return static_cast<unsigned>((digits * 10000 + (digits >> 32)) & 0xffffffffULL);
#else
    unsigned result = 0;
    for (int i = 0; i < 8; ++i) {
        result = result * 10U + static_cast<unsigned>(source[i] - '0');
    }
    return result;
#endif
}

inline long long read_int() {
    initialize();

    if (mapped != nullptr) {
        while (cursor < length && mapped[cursor] <= ' ')
            ++cursor;
        if (cursor == length)
            return 0;

        const bool negative = mapped[cursor] == '-';
        if (negative) {
            ++cursor;
            if (cursor == length)
                return 0;
        }

        if (!negative && length - cursor >= 9 && has_eight_digits(mapped + cursor)) {
            const unsigned value = parse_eight_digits(mapped + cursor);
            const unsigned ninth = static_cast<unsigned>(mapped[cursor + 8] - '0');
            if (ninth >= 10U) {
                cursor += 8;
                return static_cast<long long>(value);
            }
            if (length - cursor >= 10 && static_cast<unsigned>(mapped[cursor + 9] - '0') >= 10U) {
                cursor += 9;
                return static_cast<long long>(value * 10U + ninth);
            }
        }

        long long value = 0;
        if (negative) {
            while (length - cursor >= 2) {
                const unsigned first = static_cast<unsigned>(mapped[cursor] - '0');
                const unsigned second = static_cast<unsigned>(mapped[cursor + 1] - '0');
                if (first >= 10U || second >= 10U)
                    break;
                value = value * 100 - static_cast<long long>(first * 10U + second);
                cursor += 2;
            }
            if (cursor < length) {
                const unsigned digit = static_cast<unsigned>(mapped[cursor] - '0');
                if (digit < 10U) {
                    value = value * 10 - static_cast<long long>(digit);
                    ++cursor;
                }
            }
        } else {
            while (length - cursor >= 2) {
                const unsigned first = static_cast<unsigned>(mapped[cursor] - '0');
                const unsigned second = static_cast<unsigned>(mapped[cursor + 1] - '0');
                if (first >= 10U || second >= 10U)
                    break;
                value = value * 100 + static_cast<long long>(first * 10U + second);
                cursor += 2;
            }
            if (cursor < length) {
                const unsigned digit = static_cast<unsigned>(mapped[cursor] - '0');
                if (digit < 10U) {
                    value = value * 10 + static_cast<long long>(digit);
                    ++cursor;
                }
            }
        }
        return value;
    }

    for (;;) {
        if (cursor == length && !refill())
            return 0;
        while (cursor < length && buffer[cursor] <= ' ')
            ++cursor;
        if (cursor < length)
            break;
    }

    const bool negative = buffer[cursor] == '-';
    if (negative)
        ++cursor;
    long long value = 0;

    for (;;) {
        if (!negative && length - cursor >= 9 && has_eight_digits(buffer + cursor)) {
            const unsigned first_eight = parse_eight_digits(buffer + cursor);
            const unsigned ninth = static_cast<unsigned>(buffer[cursor + 8] - '0');
            if (ninth >= 10U) {
                cursor += 8;
                return static_cast<long long>(first_eight);
            }
            if (length - cursor >= 10 && static_cast<unsigned>(buffer[cursor + 9] - '0') >= 10U) {
                cursor += 9;
                return static_cast<long long>(first_eight * 10U + ninth);
            }
        }

        while (length - cursor >= 2) {
            const unsigned first = static_cast<unsigned>(buffer[cursor] - '0');
            const unsigned second = static_cast<unsigned>(buffer[cursor + 1] - '0');
            if (first >= 10U)
                return value;
            if (second >= 10U) {
                value = negative ? value * 10 - static_cast<long long>(first)
                                 : value * 10 + static_cast<long long>(first);
                ++cursor;
                return value;
            }
            value = negative ? value * 100 - static_cast<long long>(first * 10U + second)
                             : value * 100 + static_cast<long long>(first * 10U + second);
            cursor += 2;
        }

        if (cursor < length) {
            const unsigned digit = static_cast<unsigned>(buffer[cursor] - '0');
            if (digit >= 10U)
                return value;
            value = negative ? value * 10 - static_cast<long long>(digit)
                             : value * 10 + static_cast<long long>(digit);
            ++cursor;
        }
        if (!refill())
            return value;
    }
}

template <class T> inline void read_int_array(T *output, std::size_t count) {
    for (std::size_t i = 0; i < count; ++i) {
        output[i] = static_cast<T>(read_int());
    }
}
} // namespace input_native

namespace detail::combined_output_native {
struct FourDigits {
    char data[10000][4];

    FourDigits() {
        for (unsigned i = 0; i < 10000; ++i) {
            data[i][0] = static_cast<char>('0' + i / 1000);
            data[i][1] = static_cast<char>('0' + i / 100 % 10);
            data[i][2] = static_cast<char>('0' + i / 10 % 10);
            data[i][3] = static_cast<char>('0' + i % 10);
        }
    }
};

inline const FourDigits &four_digits() {
    static const FourDigits table;
    return table;
}

inline char *write_small(char *output, unsigned value, const FourDigits &table) {
    if (value >= 1000) {
        std::memcpy(output, table.data[value], 4);
        return output + 4;
    }
    if (value >= 100) {
        std::memcpy(output, table.data[value] + 1, 3);
        return output + 3;
    }
    if (value >= 10) {
        std::memcpy(output, table.data[value] + 2, 2);
        return output + 2;
    }
    *output++ = static_cast<char>('0' + value);
    return output;
}

template <class Unsigned>
inline char *write_unsigned(char *output, Unsigned value, const FourDigits &table) {
    unsigned chunks[5];
    unsigned count = 0;
    while (value >= 10000) {
        const Unsigned quotient = value / 10000;
        chunks[count++] = static_cast<unsigned>(value - quotient * 10000);
        value = quotient;
    }
    output = write_small(output, static_cast<unsigned>(value), table);
    while (count != 0) {
        std::memcpy(output, table.data[chunks[--count]], 4);
        output += 4;
    }
    return output;
}

template <class Integer>
inline std::size_t join_signed(const Integer *values, std::size_t count, char *output,
                               const char *separator, std::size_t separator_length) {
    const FourDigits &table = four_digits();
    char *cursor = output;
    using Unsigned = typename std::make_unsigned<Integer>::type;
    for (std::size_t i = 0; i < count; ++i) {
        const Integer value = values[i];
        Unsigned magnitude = static_cast<Unsigned>(value);
        if (value < 0) {
            *cursor++ = '-';
            magnitude = Unsigned(0) - magnitude;
        }
        cursor = write_unsigned(cursor, magnitude, table);
        if (i + 1 != count) {
            std::memcpy(cursor, separator, separator_length);
            cursor += separator_length;
        }
    }
    return static_cast<std::size_t>(cursor - output);
}

template <class Integer>
inline std::size_t join_unsigned(const Integer *values, std::size_t count, char *output,
                                 const char *separator, std::size_t separator_length) {
    const FourDigits &table = four_digits();
    char *cursor = output;
    for (std::size_t i = 0; i < count; ++i) {
        cursor = write_unsigned(cursor, values[i], table);
        if (i + 1 != count) {
            std::memcpy(cursor, separator, separator_length);
            cursor += separator_length;
        }
    }
    return static_cast<std::size_t>(cursor - output);
}
} // namespace output_native

#if defined(__clang__)
#pragma clang attribute pop
#elif defined(__GNUC__)
#pragma GCC pop_options
#endif
using namespace detail::template_support;
using detail::template_support::print;
inline constexpr Int INF = INF64;

inline Int ii() {
    return detail::combined_input_native::read_int();
}

inline std::vector<Int> lii(Int n) {
    assert(n >= 0);
    std::vector<Int> out(n);
    detail::combined_input_native::read_int_array(out.data(), out.size());
    return out;
}

inline std::vector<std::uint32_t> lii2(Int n) {
    assert(n >= 0);
    std::vector<std::uint32_t> out(n);
    detail::combined_input_native::read_int_array(out.data(), out.size());
    return out;
}

inline std::string si() {
    std::string out;
    int c = detail::combined_input_native::get_char();
    while (c >= 0 && c <= ' ')
        c = detail::combined_input_native::get_char();
    while (c > ' ') {
        out += char(c);
        c = detail::combined_input_native::get_char();
    }
    return out;
}

template <class T> std::string join(std::span<const T> values, std::string_view separator = "") {
    if (values.empty())
        return {};
    if constexpr (std::integral<T> && (sizeof(T) == 4 || sizeof(T) == 8)) {
        std::string out(values.size() * 20 + (values.size() - 1) * separator.size(), '\0');
        std::size_t n;
        if constexpr (std::is_signed_v<T>)
            n = detail::combined_output_native::join_signed(
                values.data(), values.size(), out.data(), separator.data(), separator.size());
        else
            n = detail::combined_output_native::join_unsigned(
                values.data(), values.size(), out.data(), separator.data(), separator.size());
        out.resize(n);
        return out;
    } else if constexpr (requires(const T &x) {
                             x.umod();
                             x.val();
                         }) {
        std::vector<std::uint32_t> canonical;
        canonical.reserve(values.size());
        for (const auto &x : values)
            canonical.push_back(x.val());
        return join<std::uint32_t>(canonical, separator);
    } else {
        std::string out;
        bool first = true;
        for (const auto &x : values) {
            if (!first)
                out += separator;
            first = false;
            out += stringValue(x);
        }
        return out;
    }
}

template <class T, class A>
std::string join(const std::vector<T, A> &values, std::string_view sep = "") {
    return join<T>(std::span<const T>(values.data(), values.size()), sep);
}

template <class T, std::size_t N>
std::string join(const std::array<T, N> &values, std::string_view sep = "") {
    return join<T>(std::span<const T>(values), sep);
}

struct Separator {
    std::string value;
};

inline Separator sep(std::string value) {
    return {std::move(value)};
}

template <class R> struct Splat {
    const R &values;
};

template <class R> Splat<R> splat(const R &r) {
    return {r};
}

template <class... T> void print(Separator separator, const T &...values) {
    detail::template_support::print(PrintOptions{stdout, std::move(separator.value), "\n", false},
                                    values...);
}

template <class R> void print(Splat<R> values) {
    detail::template_support::print(cplib::join(values.values, " "));
}

template <class R> void print(Separator separator, Splat<R> values) {
    detail::template_support::print(cplib::join(values.values, separator.value));
}

inline void yes(bool b = true) {
    detail::template_support::print(b ? "Yes" : "No");
}
}
