#pragma once
#include <cplib/modint/modint.hpp>

namespace cplib::detail {
template<unsigned M>
inline void subsetBarrettKernel(unsigned n, uint32_t * __restrict left,
                                uint32_t * __restrict right, const uint8_t *ranks,
                                const uint32_t *offsets) {
    const unsigned size = 1u << n;
    // 各積は(M-1)^2以下なので、この項数ごとに剰余を取ればuint64_tに収まる。
    constexpr uint64_t capacity = std::numeric_limits<uint64_t>::max() / (uint64_t(M - 1) * (M - 1));
    constexpr unsigned terms = capacity < 64 ? unsigned(capacity) : 64;
    const unsigned blockBits = std::min(n, 12u), blockSize = 1u << blockBits;
    auto zetaStage = [&](unsigned bit, unsigned section, unsigned limit) {
        for (unsigned base = section; base < limit; base += bit * 2) {
            const unsigned first = ranks[base];
            for (unsigned j = 0; j < bit; ++j) {
                uint32_t* llo = left + offsets[base + j];
                uint32_t* lhi = left + offsets[base + j + bit];
                uint32_t* rlo = right + offsets[base + j];
                uint32_t* rhi = right + offsets[base + j + bit];
                for (unsigned d = first; d <= unsigned(ranks[base + j]); ++d) {
                    const uint32_t lv = lhi[d] + llo[d], rv = rhi[d] + rlo[d];
                    lhi[d] = lv >= M ? lv - M : lv;
                    rhi[d] = rv >= M ? rv - M : rv;
                }
            }
        }
    };
    for (unsigned section = 0; section < size; section += blockSize)
        for (unsigned h = 0; h < blockBits; ++h) zetaStage(1u << h, section, section + blockSize);
    for (unsigned h = blockBits; h < n; ++h) zetaStage(1u << h, 0, size);
    for (unsigned mask = 0; mask < size; ++mask) {
        uint32_t* lrow = left + offsets[mask];
        const uint32_t* rrow = right + offsets[mask];
        const unsigned rank = ranks[mask];
        // 入力係数を壊さず積を求めてから、必要な次数だけを同じ行に詰める。
        uint32_t output[21];
        const unsigned lastDegree = std::min(n, 2 * rank);
        for (unsigned d = rank; d <= lastDegree; ++d) {
            const unsigned first = d - rank;
            if (rank - first + 1 <= terms) {
                uint64_t sum = 0;
                for (unsigned k = first; k <= rank; ++k) sum += uint64_t(lrow[k]) * rrow[d - k];
                output[d - rank] = uint32_t(sum % M);
            } else {
                uint64_t value = 0;
                for (unsigned k = first; k <= rank;) {
                    const unsigned last = std::min(rank, k + terms - 1);
                    uint64_t sum = 0;
                    for (unsigned j = k; j <= last; ++j) sum += uint64_t(lrow[j]) * rrow[d - j];
                    value += sum % M;
                    k = last + 1;
                }
                output[d - rank] = uint32_t(value % M);
            }
        }
        for (unsigned d = rank; d <= lastDegree; ++d) lrow[d - rank] = output[d - rank];
    }
    auto mobiusStage = [&](unsigned bit, unsigned section, unsigned limit) {
        for (unsigned base = section; base < limit; base += 2 * bit) {
            for (unsigned j = 0; j < bit; ++j) {
                const uint32_t* lo = left + offsets[base + j];
                uint32_t* hi = left + offsets[base + j + bit];
                const unsigned rank = ranks[base + j + bit];
                if (rank < 2) continue;
                const unsigned last = std::min(n - rank, rank - 2);
                for (unsigned d = 0; d <= last; ++d) {
                    // 両係数はM未満、M<=2^30なので最上位bitで負の差を判別できる。
                    const uint32_t v = hi[d] - lo[d + 1];
                    hi[d] = (v & 0x80000000u) ? v + M : v;
                }
            }
        }
    };
    for (unsigned section = 0; section < size; section += blockSize)
        for (unsigned h = 0; h < blockBits; ++h) mobiusStage(1u << h, section, section + blockSize);
    for (unsigned h = blockBits; h < n; ++h) mobiusStage(1u << h, 0, size);
}

    }
namespace cplib {
// 静的Barrett型の部分集合畳み込み。長さ2^n、n<=20、1<法<=2^30が必要。
// 時間 O(n²2^n)、領域 O(n2^n)。rankごとの係数を詰めて保持する。
template <std::uint32_t M>
std::vector<StaticBarrettModint<M>> subsetConvolutionBarrett(
    const std::vector<StaticBarrettModint<M>> &a, const std::vector<StaticBarrettModint<M>> &b) {
    static_assert(M > 1 && M <= (1u << 30));
    assert(a.size() == b.size());
    if (a.empty())
        return {};
    assert(std::has_single_bit(a.size()));
    unsigned n = std::countr_zero(a.size());
    assert(n <= 20);
    std::vector<std::uint8_t> ranks(a.size());
    std::vector<std::uint32_t> offsets(a.size() + 1);
    for (std::size_t mask = 0; mask < a.size(); ++mask) {
        ranks[mask] = std::popcount(mask);
        offsets[mask + 1] = offsets[mask] + ranks[mask] + 1;
    }
    std::vector<std::uint32_t> left(offsets.back()), right(offsets.back());
    for (std::size_t mask = 0; mask < a.size(); ++mask) {
        auto index = offsets[mask] + ranks[mask];
        left[index] = a[mask].val();
        right[index] = b[mask].val();
    }
    detail::subsetBarrettKernel<M>(n, left.data(), right.data(), ranks.data(), offsets.data());
    std::vector<StaticBarrettModint<M>> result(a.size());
    for (std::size_t mask = 0; mask < a.size(); ++mask)
        result[mask] = left[offsets[mask]];
    return result;
}
}
