#pragma once
#include <cplib/common.hpp>
#include <iterator>

namespace cplib {
namespace detail {
// すべての列挙器は一要素O(1)、補助領域O(1)。
struct BitRange {
    enum Mode {
        Combination,
        SubsetEq,
        Subset,
        SubsetEqDesc,
        SubsetDesc,
        SupersetEq,
        Superset,
        Singleton,
        Standing
    };

    Mode mode;
    UInt bits, limit, first;
    bool empty = false;

    struct Iterator {
        const BitRange *range;
        UInt value;
        bool done;

        Int operator*() const {
            if (range->mode == Standing)
                return std::countr_zero(value);
            if (range->mode == Singleton)
                return std::bit_cast<Int>(value & (0 - value));
            return std::bit_cast<Int>(value);
        }

        Iterator &operator++() {
            auto &r = *range;
            switch (r.mode) {
            case Combination: {
                if (value == 0) {
                    done = true;
                    break;
                }
                UInt t = value | (value - 1), low = (~t) & (0 - (~t));
                value = (t + 1) | ((low - 1) >> (std::countr_zero(value) + 1));
                done = value >= r.limit;
                break;
            }
            case SubsetEq:
                if (value == r.bits)
                    done = true;
                else
                    value = (value - r.bits) & r.bits;
                break;
            case Subset:
                value = (value - r.bits) & r.bits;
                done = value == r.bits;
                break;
            case SubsetEqDesc:
            case SubsetDesc:
                if (value == 0)
                    done = true;
                else
                    value = (value - 1) & r.bits;
                break;
            case SupersetEq:
            case Superset:
                value = (value + 1) | r.bits;
                done = value >= r.limit;
                break;
            case Singleton:
                value &= value - 1;
                done = value == 0;
                break;
            case Standing:
                value = r.bits & (~r.bits + (value << 1));
                done = value == 0;
                break;
            }
            return *this;
        }

        bool operator!=(std::default_sentinel_t) const {
            return !done;
        }
    };

    Iterator begin() const {
        return {this, first, empty};
    }

    std::default_sentinel_t end() const {
        return {};
    }
};
}

// n bit中 r bitが1であるようなbit列を列挙します。
inline auto bitcomb(Int n, Int r) {
    assert(0 <= r && r <= n && n <= 63);
    return detail::BitRange{detail::BitRange::Combination, 0, UInt(1) << n, (UInt(1) << r) - 1};
}

// 与えられた集合の部分集合を昇順で列挙します。与えられた集合も含みます。
inline auto bitsubseteq(Int bits) {
    return detail::BitRange{detail::BitRange::SubsetEq, UInt(bits), 0, 0};
}

// 非空のbitsの真部分集合を昇順で列挙する。bits==0の場合は0を一件返す。
inline auto bitsubset(Int bits) {
    return detail::BitRange{detail::BitRange::Subset, UInt(bits), 0, 0};
}

// 与えられた集合の部分集合を降順で列挙します。与えられた集合も含みます。
inline auto bitsubseteq_descending(Int bits) {
    return detail::BitRange{detail::BitRange::SubsetEqDesc, UInt(bits), 0, UInt(bits)};
}

// 非空のbitsの真部分集合を降順で列挙する。bits==0の場合は0を一件返す。
inline auto bitsubset_descending(Int bits) {
    return detail::BitRange{detail::BitRange::SubsetDesc, UInt(bits), 0,
                            (UInt(bits) - 1) & UInt(bits)};
}

// 与えられた集合を包含する集合(上位集合)を列挙します。与えられた集合も含みます。bit数上限をnとします。
inline auto bitsuperseteq(Int bits, Int n) {
    assert(0 <= n && n <= 63);
    return detail::BitRange{detail::BitRange::SupersetEq, UInt(bits), UInt(1) << n, UInt(bits)};
}

// 与えられた集合を包含する集合(上位集合)を列挙します。与えられた集合は含みません。bit数上限をnとします。
inline auto bitsuperset(Int bits, Int n) {
    assert(0 <= n && n <= 63);
    UInt first = (UInt(bits) + 1) | UInt(bits);
    return detail::BitRange{detail::BitRange::Superset, UInt(bits), UInt(1) << n, first,
                            first >= (UInt(1) << n)};
}

// 立っているbitを一つずつ取り出します。
inline auto bitsingleton(Int bits) {
    return detail::BitRange{detail::BitRange::Singleton, UInt(bits), 0, UInt(bits), bits == 0};
}

inline auto standingbits(Int bits) {
    UInt first = UInt(bits) & (0 - UInt(bits));
    return detail::BitRange{detail::BitRange::Standing, UInt(bits), 0, first, first == 0};
}
}
