#pragma once
#include <cplib/math/isqrt.hpp>

namespace cplib {
namespace detail {
inline constexpr std::array<Int, 8> wheel_residues{1, 7, 11, 13, 17, 19, 23, 29};
inline constexpr std::array<int, 30> wheel_index{-1, 0,  -1, -1, -1, -1, -1, 1,  -1, -1,
                                                 -1, 2,  -1, 3,  -1, -1, -1, 4,  -1, 5,
                                                 -1, -1, -1, 6,  -1, -1, -1, -1, -1, 7};
inline constexpr std::array<Int, 11> presieve_primes{7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43};

template <std::size_t N> consteval auto make_presieve(std::array<Int, 3> group) {
    std::array<std::uint8_t, N> pattern{};
    for (std::size_t q = 0; q < N; ++q) {
        std::uint8_t bits = 255;
        for (int i = 0; i < 8; ++i)
            for (Int p : group)
                if (p > 1 && (30 * q + wheel_residues[i]) % p == 0)
                    bits &= ~(1u << i);
        pattern[q] = bits;
    }
    return pattern;
}

inline constexpr auto presieve0 = make_presieve<7 * 11 * 13>({7, 11, 13});
inline constexpr auto presieve1 = make_presieve<17 * 19>({17, 19, 1});
inline constexpr auto presieve2 = make_presieve<23 * 29>({23, 29, 1});
inline constexpr auto presieve3 = make_presieve<31 * 37>({31, 37, 1});
inline constexpr auto presieve4 = make_presieve<41 * 43>({41, 43, 1});
inline constexpr std::array<std::span<const std::uint8_t>, 5> presieves{
    presieve0, presieve1, presieve2, presieve3, presieve4};

// 基底素数を奇数篩で列挙する。O(limit log log limit)。
inline std::vector<Int> base_primes(Int limit) {
    std::vector<bool> composite(limit / 2 + 1);
    for (Int p = 3; p <= limit / p; p += 2)
        if (!composite[p / 2])
            for (Int j = p * p / 2; j < static_cast<Int>(composite.size()); j += p)
                composite[j] = true;
    std::vector<Int> result;
    for (Int i = 1; i < static_cast<Int>(composite.size()); ++i)
        if (2 * i + 1 <= limit && !composite[i])
            result.push_back(2 * i + 1);
    return result;
}

struct SieveStrike {
    Int prime;
    std::array<Int, 8> next;
    std::array<std::uint8_t, 8> masks;
};

// 周期パターンを連続領域に合成する。O(size)。
inline void fill_sieve_block(std::uint8_t *dst, Int first_byte, Int size) {
    for (std::size_t g = 0; g < presieves.size(); ++g) {
        auto pattern = presieves[g];
        Int period = pattern.size(), phase = first_byte % period, pos = 0;
        while (pos < size) {
            Int n = std::min(period - phase, size - pos);
            if (g == 0)
                std::copy_n(pattern.data() + phase, n, dst + pos);
            else
                for (Int i = 0; i < n; ++i)
                    dst[pos + i] &= pattern[phase + i];
            pos += n;
            phase = 0;
        }
    }
}

// 倍数を消し、次のブロックに続く位置を保存する。
inline void strike_sieve_block(std::uint8_t *dst, Int size, SieveStrike &state) {
    Int p = state.prime;
    for (int r = 0; r < 8; ++r) {
        Int j = state.next[r];
        auto mask = state.masks[r];
        while (j + 3 * p < size) {
            dst[j] &= mask;
            dst[j + p] &= mask;
            dst[j + 2 * p] &= mask;
            dst[j + 3 * p] &= mask;
            j += 4 * p;
        }
        while (j < size) {
            dst[j] &= mask;
            j += p;
        }
        state.next[r] = j - size;
    }
}
}

class EratosthenesSieve {
    Int first = 0, last = 0, first_byte = 0;
    std::vector<std::uint8_t> bits;

public:
    // 閉区間を30輪とブロック篩で構築する。保持領域は約(high-low)/30 byte。
    EratosthenesSieve(Int low, Int high) : first(low), last(high), first_byte(low / 30) {
        assert(0 <= low && low <= high);
        bits.resize(high / 30 - first_byte + 1);
        std::vector<detail::SieveStrike> states;
        for (Int p : detail::base_primes(isqrt(high))) {
            if (p <= 43)
                continue;
            detail::SieveStrike state;
            state.prime = p;
            for (int i = 0; i < 8; ++i) {
                UInt multiplier = UInt(p) + UInt((detail::wheel_residues[i] - p % 30 + 30) % 30);
                UInt product = UInt(p) * multiplier, q = product / 30;
                if (q < UInt(first_byte)) {
                    UInt gap = UInt(first_byte) - q;
                    q += (gap + UInt(p) - 1) / UInt(p) * UInt(p);
                }
                state.next[i] = static_cast<Int>(q) - first_byte;
                state.masks[i] =
                    static_cast<std::uint8_t>(~(1u << detail::wheel_index[product % 30]));
            }
            states.push_back(state);
        }
        for (Int pos = 0; pos < static_cast<Int>(bits.size());) {
            Int size = std::min<Int>(32768, bits.size() - pos);
            detail::fill_sieve_block(bits.data() + pos, first_byte + pos, size);
            for (auto &state : states)
                detail::strike_sieve_block(bits.data() + pos, size, state);
            pos += size;
        }
        for (Int p : detail::presieve_primes)
            if (low <= p && p <= high)
                bits[p / 30 - first_byte] |= 1u << detail::wheel_index[p % 30];
        for (int i = 0; i < 8; ++i) {
            if (detail::wheel_residues[i] < low % 30)
                bits.front() &= ~(1u << i);
            if (detail::wheel_residues[i] > high % 30)
                bits.back() &= ~(1u << i);
        }
        if (low <= 1 && 1 <= high)
            bits.front() &= ~1u;
    }

    // 構築した閉区間内の素数判定を O(1) 時間で行う。範囲外は false。
    bool is_prime(Int n) const {
        if (n < first || n > last || n < 2)
            return false;
        if (n == 2 || n == 3 || n == 5)
            return true;
        int bit = detail::wheel_index[n % 30];
        return bit >= 0 && (bits[n / 30 - first_byte] & (1u << bit));
    }

    // 圧縮配列のバイト数を返す。O(1)。
    Int byte_size() const {
        return bits.size();
    }

    // 素数の数を数える。O(区間長/30)。
    Int count_primes() const {
        Int count = 0;
        for (Int p : {2, 3, 5})
            if (first <= p && p <= last)
                ++count;
        for (auto b : bits)
            count += std::popcount(b);
        return count;
    }

    class iterator {
        const EratosthenesSieve *sieve = nullptr;
        std::size_t q = 0;
        unsigned remaining = 0;
        int small = 0;
        Int current = 0;
        bool done = true;

        // 次の素数へ進める。列挙全体でO(配列長+素数個数)。
        void advance() {
            constexpr Int primes[]{2, 3, 5};
            while (small < 3) {
                Int p = primes[small++];
                if (sieve->first <= p && p <= sieve->last) {
                    current = p;
                    return;
                }
            }
            while (!remaining) {
                if (q == sieve->bits.size()) {
                    done = true;
                    return;
                }
                remaining = sieve->bits[q++];
            }
            int bit = std::countr_zero(remaining);
            remaining &= remaining - 1;
            current =
                (sieve->first_byte + static_cast<Int>(q) - 1) * 30 + detail::wheel_residues[bit];
        }

    public:
        using value_type = Int;
        using difference_type = std::ptrdiff_t;
        using iterator_category = std::input_iterator_tag;
        iterator() = default;

        explicit iterator(const EratosthenesSieve *s) : sieve(s), done(false) {
            advance();
        }

        Int operator*() const {
            return current;
        }

        iterator &operator++() {
            advance();
            return *this;
        }

        iterator operator++(int) {
            auto old = *this;
            advance();
            return old;
        }

        bool operator==(const iterator &rhs) const {
            return done || rhs.done ? done == rhs.done
                                    : sieve == rhs.sieve && current == rhs.current;
        }
    };

    iterator begin() const {
        return iterator(this);
    }

    iterator end() const {
        return iterator();
    }
};

// 閉区間[low,high]の篩を構築する。
inline EratosthenesSieve initSegmentedEratosthenes(Int low, Int high) {
    return {low, high};
}

// 閉区間[0,limit]の篩を構築する。
inline EratosthenesSieve initEratosthenes(Int limit) {
    return {0, limit};
}

inline bool is_prime(const EratosthenesSieve &s, Int n) {
    return s.is_prime(n);
}

inline Int byte_size(const EratosthenesSieve &s) {
    return s.byte_size();
}

inline Int count_primes(const EratosthenesSieve &s) {
    return s.count_primes();
}

// 区間内の素数を昇順に O(区間長 / 30 + 素数の個数) 時間で列挙する。
inline const EratosthenesSieve &items(const EratosthenesSieve &s) {
    return s;
}

// limit以下の素数を昇順のvectorで返す。limit < 2なら空列。
inline std::vector<Int> get_primes(Int limit) {
    if (limit < 2)
        return {};
    auto sieve = initEratosthenes(limit);
    std::vector<Int> result;
    result.reserve(sieve.count_primes());
    for (Int p : sieve)
        result.push_back(p);
    return result;
}

// 半開区間[l, r)の素数を区間篩で昇順のvectorにする。l >= rなら空列。
inline std::vector<Int> get_primes(Int l, Int r) {
    if (l >= r || r <= 2)
        return {};
    auto sieve = initSegmentedEratosthenes(std::max<Int>(l, 2), r - 1);
    std::vector<Int> result;
    result.reserve(sieve.count_primes());
    for (Int p : sieve)
        result.push_back(p);
    return result;
}
}
