#pragma once
#include <cplib/collections/bitvector.hpp>
#include <optional>

namespace cplib {
struct WaveletChildren {
    Int l0, r0, l1, r1;
};

struct WaveletSumCount {
    Int sum = 0, count = 0;
    bool operator==(const WaveletSumCount &) const = default;
};

class WaveletMatrix {
    struct Level {
        BitVector bits;
        Int zero_count = 0;
        std::vector<Int> zero_sum;
    };

    std::vector<Level> dat;
    Int H, N;
    bool with_sum;
    Int scan_h;
    std::vector<Int> scan_values;

    void validateSum(Int l, Int r) const {
        assert(with_sum && 0 <= l && l <= r && r <= N);
    }

    template <bool Sum, bool Inclusive> auto bound(Int l, Int r, Int x) const {
        WaveletSumCount out;
        if constexpr (Sum)
            validateSum(l, r);
        if (x < 0 || (!Inclusive && x == 0))
            return out;
        if (H < 64 && (UInt(x) >> H) != 0) {
            out.count = r - l;
            if constexpr (Sum)
                out.sum = sum_smallest(l, r, r - l);
            return out;
        }
        for (Int h = H - 1; h >= 0; --h) {
            if (l == r)
                return out;
            if (h == scan_h && r - l <= 64) {
                for (Int i = l; i < r; ++i)
                    if (Inclusive ? scan_values[i] <= x : scan_values[i] < x) {
                        ++out.count;
                        if constexpr (Sum)
                            out.sum += scan_values[i];
                    }
                return out;
            }
            auto [l0, r0, l1, r1] = get_child(h, l, r);
            if ((UInt(x) >> h) & 1) {
                out.count += r0 - l0;
                if constexpr (Sum)
                    out.sum += dat[h].zero_sum[r0] - dat[h].zero_sum[l0];
                l = l1;
                r = r1;
            } else {
                l = l0;
                r = r0;
            }
        }
        if constexpr (Inclusive) {
            out.count += r - l;
            if constexpr (Sum)
                out.sum += (r - l) * x;
        }
        return out;
    }

    bool share_range_path(Int low, Int high) const {
        return low > 0 && ((UInt(low) ^ UInt(high - 1)) >> std::max<Int>(H - 4, 0)) == 0;
    }

    template <bool Sum, bool Count, bool Inclusive>
    WaveletSumCount bound_from(Int l, Int r, Int x, Int start) const {
        WaveletSumCount out;
        for (Int h = start; h >= 0; --h) {
            if (l == r)
                return out;
            if (h == scan_h && r - l <= 64) {
                for (Int i = l; i < r; ++i)
                    if (Inclusive ? scan_values[i] <= x : scan_values[i] < x) {
                        if constexpr (Sum)
                            out.sum += scan_values[i];
                        if constexpr (Count)
                            ++out.count;
                    }
                return out;
            }
            auto [l0, r0, l1, r1] = get_child(h, l, r);
            if ((UInt(x) >> h) & 1) {
                if constexpr (Sum)
                    out.sum += dat[h].zero_sum[r0] - dat[h].zero_sum[l0];
                if constexpr (Count)
                    out.count += r0 - l0;
                l = l1;
                r = r1;
            } else {
                l = l0;
                r = r0;
            }
        }
        if constexpr (Inclusive) {
            if constexpr (Sum)
                out.sum += (r - l) * x;
            if constexpr (Count)
                out.count += r - l;
        }
        return out;
    }

    template <bool Sum, bool Count>
    WaveletSumCount range_query(Int l, Int r, Int low, Int high) const {
        WaveletSumCount out;
        if (low >= high || high <= 0 || l == r)
            return out;
        Int lower = std::max<Int>(low, 0), upper = high - 1;
        if (H < 63)
            upper = std::min(upper, (Int(1) << H) - 1);
        if (lower > upper)
            return out;
        Int split = lower == upper ? -1 : Int(std::bit_width(UInt(lower) ^ UInt(upper))) - 1;
        for (Int h = H - 1; h >= 0; --h) {
            if (l == r)
                return out;
            if (h == scan_h && r - l <= 64) {
                for (Int i = l; i < r; ++i)
                    if (lower <= scan_values[i] && scan_values[i] <= upper) {
                        if constexpr (Sum)
                            out.sum += scan_values[i];
                        if constexpr (Count)
                            ++out.count;
                    }
                return out;
            }
            auto [l0, r0, l1, r1] = get_child(h, l, r);
            if (h == split) {
                auto a = bound_from<Sum, Count, false>(l0, r0, lower, h - 1),
                     b = bound_from<Sum, Count, true>(l1, r1, upper, h - 1);
                if constexpr (Sum)
                    out.sum = dat[h].zero_sum[r0] - dat[h].zero_sum[l0] - a.sum + b.sum;
                if constexpr (Count)
                    out.count = r0 - l0 - a.count + b.count;
                return out;
            }
            if ((UInt(lower) >> h) & 1) {
                l = l1;
                r = r1;
            } else {
                l = l0;
                r = r0;
            }
        }
        if constexpr (Sum)
            out.sum = (r - l) * lower;
        if constexpr (Count)
            out.count = r - l;
        return out;
    }

public:
    // 非負整数列をO(NH)で構築。ワード単位のビット設定、分岐を省く振分け、
    // 64要素以下の下位層を直接走査する短縮経路を元実装どおり保持する。
    explicit WaveletMatrix(std::span<const Int> input = {}, Int height = -1, bool sums = false)
        : H(height), N(input.size()), with_sum(sums) {
        std::vector<Int> v(input.begin(), input.end());
        if (H == -1)
            H = N == 0 ? 0
                : *std::max_element(v.begin(), v.end()) == 0
                    ? 1
                    : std::bit_width(UInt(*std::max_element(v.begin(), v.end())));
        assert(0 <= H && H <= 64);
        dat.resize(H);
        for (auto &level : dat)
            level.bits = BitVector(N);
        Int width = N, depth = 0;
        while (width > 64) {
            width = (width >> 1) + (width & 1);
            ++depth;
        }
        scan_h = H - depth - 1;
        if (scan_h < 6)
            scan_h = -1;
        std::vector<Int> zero(N), one(N);
        for (Int h = H - 1; h >= 0; --h) {
            if (h == scan_h)
                scan_values = v;
            UInt mask = UInt(1) << h, word = 0;
            Int a = 0, b = 0;
            for (Int i = 0; i < N; ++i) {
                Int value = v[i], bit = (UInt(value) & mask) != 0;
                zero[a] = value;
                one[b] = value;
                a += 1 - bit;
                b += bit;
                word |= UInt(bit) << (i & 63);
                if ((i & 63) == 63) {
                    dat[h].bits.setWord(i >> 6, word);
                    word = 0;
                }
            }
            if (N & 63)
                dat[h].bits.setWord(N >> 6, word);
            if (with_sum) {
                dat[h].zero_sum.resize(a + 1);
                Int total = 0;
                for (Int i = 0; i < a; ++i) {
                    total += zero[i];
                    dat[h].zero_sum[i + 1] = total;
                }
            }
            dat[h].zero_count = a;
            std::copy_n(zero.begin(), a, v.begin());
            std::copy_n(one.begin(), b, v.begin() + a);
            dat[h].bits.build();
        }
    }

    Int size() const {
        return N;
    }

    WaveletChildren get_child(Int h, Int l, Int r) const {
        Int c0 = dat[h].zero_count, lr = dat[h].bits.rank(l), rr = dat[h].bits.rank(r);
        return {l - lr, r - rr, lr + c0, rr + c0};
    }

    // 照会はO(H)。H=0も元の空・全零列の意味を維持。
    // 半開区間[l,r)内で小さい順にk番目の値をO(H)で返す。0 <= k < r-lとし、kは0-indexed。Hは格納値のビット数。
    Int kth_smallest(Int l, Int r, Int k) const {
        UInt out = 0;
        for (Int h = H - 1; h >= 0; --h) {
            auto [l0, r0, l1, r1] = get_child(h, l, r);
            if (k < r0 - l0) {
                l = l0;
                r = r0;
            } else {
                l = l1;
                r = r1;
                k -= r0 - l0;
                out += UInt(1) << h;
            }
        }
        return std::bit_cast<Int>(out);
    }

    Int kth_largest(Int l, Int r, Int k) const {
        assert(0 <= k && k < r - l);
        return kth_smallest(l, r, r - l - 1 - k);
    }

    Int range_lowerbound(Int l, Int r, Int x) const {
        return bound<false, false>(l, r, x).count;
    }

    Int range_upperbound(Int l, Int r, Int x) const {
        return bound<false, true>(l, r, x).count;
    }

    // 半開区間[l,r)内のx未満の最大値をO(H)で返す。存在しなければ空のstd::optional。
    std::optional<Int> prev_value(Int l, Int r, Int x) const {
        Int c = range_lowerbound(l, r, x);
        if (c == 0)
            return {};
        return kth_smallest(l, r, c - 1);
    }

    // 半開区間[l,r)内のx以上の最小値をO(H)で返す。存在しなければ空のstd::optional。
    std::optional<Int> next_value(Int l, Int r, Int x) const {
        Int c = range_lowerbound(l, r, x);
        if (c == r - l)
            return {};
        return kth_smallest(l, r, c);
    }

    // 半開区間[l,r)内で値が[low,high)に入る要素数をO(H)で返す。low >= highなら0。
    Int range_freq(Int l, Int r, Int low, Int high) const {
        if (low < high && share_range_path(low, high))
            return range_query<false, true>(l, r, low, high).count;
        return low >= high ? 0 : range_lowerbound(l, r, high) - range_lowerbound(l, r, low);
    }

    Int count(Int l, Int r, Int x) const {
        if (x < 0 || (H < 64 && (UInt(x) >> H) != 0))
            return 0;
        for (Int h = H - 1; h >= 0; --h) {
            if (l == r)
                return 0;
            if (h == scan_h && r - l <= 64) {
                Int out = 0;
                for (Int i = l; i < r; ++i)
                    out += scan_values[i] == x;
                return out;
            }
            auto [l0, r0, l1, r1] = get_child(h, l, r);
            if ((UInt(x) >> h) & 1) {
                l = l1;
                r = r1;
            } else {
                l = l0;
                r = r0;
            }
        }
        return r - l;
    }

    // 半開区間[l,r)内の小さい方からk個の総和をO(H)で返す。0 <= k <= r-l。
    // 構築時にsums=trueを指定する。総和の中間値もIntに収まること。
    Int sum_smallest(Int l, Int r, Int k) const {
        validateSum(l, r);
        assert(0 <= k && k <= r - l);
        Int out = 0;
        UInt value = 0;
        for (Int h = H - 1; h >= 0; --h) {
            if (k == 0)
                return out;
            auto [l0, r0, l1, r1] = get_child(h, l, r);
            if (k < r0 - l0) {
                l = l0;
                r = r0;
            } else {
                out += dat[h].zero_sum[r0] - dat[h].zero_sum[l0];
                k -= r0 - l0;
                value += UInt(1) << h;
                l = l1;
                r = r1;
            }
        }
        return out + k * std::bit_cast<Int>(value);
    }

    Int sum_upperbound(Int l, Int r, Int x) const {
        return bound<true, true>(l, r, x).sum;
    }

    Int sum_lowerbound(Int l, Int r, Int x) const {
        validateSum(l, r);
        return x <= 0 ? 0 : sum_upperbound(l, r, x - 1);
    }

    // 半開区間[l,r)内で値が[low,high)に入る要素の総和をO(H)で返す。low >= highなら0。
    // 構築時にsums=trueを指定する。総和の中間値もIntに収まること。
    Int range_sum(Int l, Int r, Int low, Int high) const {
        validateSum(l, r);
        if (low < high && share_range_path(low, high))
            return range_query<true, false>(l, r, low, high).sum;
        return low >= high ? 0 : sum_lowerbound(l, r, high) - sum_lowerbound(l, r, low);
    }

    WaveletSumCount sum_smallest_with_count(Int l, Int r, Int k) const {
        return {sum_smallest(l, r, k), k};
    }

    WaveletSumCount sum_upperbound_with_count(Int l, Int r, Int x) const {
        return bound<true, true>(l, r, x);
    }

    WaveletSumCount sum_lowerbound_with_count(Int l, Int r, Int x) const {
        validateSum(l, r);
        return x <= 0 ? WaveletSumCount{} : sum_upperbound_with_count(l, r, x - 1);
    }

    WaveletSumCount range_sum_with_count(Int l, Int r, Int low, Int high) const {
        validateSum(l, r);
        if (low >= high)
            return {};
        if (share_range_path(low, high))
            return range_query<true, true>(l, r, low, high);
        auto a = sum_lowerbound_with_count(l, r, high), b = sum_lowerbound_with_count(l, r, low);
        return {a.sum - b.sum, a.count - b.count};
    }
};

inline auto initWaveletMatrix(std::span<const Int> v, Int H = -1, bool with_sum = false) {
    return WaveletMatrix(v, H, with_sum);
}

inline auto get_child(const WaveletMatrix &w, Int h, Int l, Int r) {
    return w.get_child(h, l, r);
}

inline auto kth_smallest(const WaveletMatrix &w, Int l, Int r, Int k) {
    return w.kth_smallest(l, r, k);
}

inline auto kth_largest(const WaveletMatrix &w, Int l, Int r, Int k) {
    return w.kth_largest(l, r, k);
}

inline auto range_lowerbound(const WaveletMatrix &w, Int l, Int r, Int x) {
    return w.range_lowerbound(l, r, x);
}

inline auto range_upperbound(const WaveletMatrix &w, Int l, Int r, Int x) {
    return w.range_upperbound(l, r, x);
}

inline auto prev_value(const WaveletMatrix &w, Int l, Int r, Int x) {
    return w.prev_value(l, r, x);
}

inline auto next_value(const WaveletMatrix &w, Int l, Int r, Int x) {
    return w.next_value(l, r, x);
}

inline auto count(const WaveletMatrix &w, Int l, Int r, Int x) {
    return w.count(l, r, x);
}

inline auto sum_smallest(const WaveletMatrix &w, Int l, Int r, Int k) {
    return w.sum_smallest(l, r, k);
}

inline auto sum_upperbound(const WaveletMatrix &w, Int l, Int r, Int x) {
    return w.sum_upperbound(l, r, x);
}

inline auto sum_lowerbound(const WaveletMatrix &w, Int l, Int r, Int x) {
    return w.sum_lowerbound(l, r, x);
}

inline auto sum_smallest_with_count(const WaveletMatrix &w, Int l, Int r, Int k) {
    return w.sum_smallest_with_count(l, r, k);
}

inline auto sum_upperbound_with_count(const WaveletMatrix &w, Int l, Int r, Int x) {
    return w.sum_upperbound_with_count(l, r, x);
}

inline auto sum_lowerbound_with_count(const WaveletMatrix &w, Int l, Int r, Int x) {
    return w.sum_lowerbound_with_count(l, r, x);
}

inline auto range_freq(const WaveletMatrix &w, Int l, Int r, Int low, Int high) {
    return w.range_freq(l, r, low, high);
}

inline auto range_sum(const WaveletMatrix &w, Int l, Int r, Int low, Int high) {
    return w.range_sum(l, r, low, high);
}

inline auto range_sum_with_count(const WaveletMatrix &w, Int l, Int r, Int low, Int high) {
    return w.range_sum_with_count(l, r, low, high);
}
}
