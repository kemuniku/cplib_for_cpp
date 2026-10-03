#pragma once
#include <cplib/common.hpp>

namespace cplib::detail {
template <std::integral K> UInt compressedCoordinateKey(K x) {
    if constexpr (std::is_signed_v<K>)
        return UInt(Int(x)) ^ (UInt(1) << 63);
    else
        return UInt(x);
}

// 大きな整数配列は11bitごとの基数ソート。異なるビットを持つ桁だけ処理する。
template <class K> void sortCompressedCoordinates(std::vector<K> &xs) {
    if constexpr (std::integral<K>) {
        if (xs.size() >= 2048) {
            UInt first = compressedCoordinateKey(xs[0]), differing = 0;
            bool ordered = true;
            for (std::size_t i = 1; i < xs.size(); ++i) {
                differing |= first ^ compressedCoordinateKey(xs[i]);
                if (xs[i] < xs[i - 1])
                    ordered = false;
            }
            if (ordered)
                return;
            std::vector<K> tmp(xs.size());
            std::array<Int, 2048> counts;
            for (int shift = 0; shift <= 55; shift += 11) {
                if (((differing >> shift) & 2047) == 0)
                    continue;
                counts.fill(0);
                for (K x : xs)
                    ++counts[(compressedCoordinateKey(x) >> shift) & 2047];
                Int total = 0;
                for (auto &count : counts) {
                    Int size = count;
                    count = total;
                    total += size;
                }
                for (K x : xs) {
                    auto bucket = (compressedCoordinateKey(x) >> shift) & 2047;
                    tmp[counts[bucket]++] = x;
                }
                xs.swap(tmp);
            }
            return;
        }
    }
    std::sort(xs.begin(), xs.end());
}

template <std::integral K> UInt compressedCoordinateHash(K x) {
    UInt h = compressedCoordinateKey(x);
    h = (h ^ (h >> 30)) * 0xbf58476d1ce4e5b9ULL;
    h = (h ^ (h >> 27)) * 0x94d049bb133111ebULL;
    return h ^ (h >> 31);
}

// 8回までの線形探索索引。入り切らないキーは二分探索へフォールバック。
template <class K> std::vector<Int> initCompressedCoordinateIndex(std::span<const K> coords) {
    std::vector<Int> result;
    if constexpr (std::integral<K>) {
        if (coords.size() >= 64) {
            Int cap = 1;
            while (cap < Int(coords.size()) * 2)
                cap *= 2;
            result.resize(cap);
            for (Int i = 0; i < Int(coords.size()); ++i) {
                Int slot = compressedCoordinateHash(coords[i]) & (cap - 1);
                for (int probe = 0; probe < 8; ++probe) {
                    if (result[slot] == 0) {
                        result[slot] = i + 1;
                        break;
                    }
                    slot = (slot + 1) & (cap - 1);
                }
            }
        }
    }
    return result;
}

template <class K>
Int findCompressedCoordinate(std::span<const K> coords, std::span<const Int> slots, const K &x) {
    if constexpr (std::integral<K>) {
        if (!slots.empty()) {
            Int mask = slots.size() - 1, slot = compressedCoordinateHash(x) & mask;
            for (int probe = 0; probe < 8; ++probe) {
                Int entry = slots[slot];
                if (entry == 0)
                    return -1;
                if (coords[entry - 1] == x)
                    return entry - 1;
                slot = (slot + 1) & mask;
            }
        }
    }
    Int i = std::lower_bound(coords.begin(), coords.end(), x) - coords.begin();
    return i < Int(coords.size()) && coords[i] == x ? i : -1;
}

template <class K> auto uniqueCompressedCoordinates(std::span<const K> coords) {
    std::vector<K> xs(coords.begin(), coords.end());
    sortCompressedCoordinates(xs);
    xs.erase(std::unique(xs.begin(), xs.end()), xs.end());
    return xs;
}
}
