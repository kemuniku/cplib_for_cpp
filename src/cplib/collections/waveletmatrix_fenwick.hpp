#pragma once
#include <cplib/collections/waveletmatrix.hpp>
#include <cplib/collections/fenwick_avx2.hpp>

namespace cplib {
class WaveletMatrixFenwick {
    WaveletMatrix matrix;
    std::vector<Int> keys, weights;
    std::vector<FenwickTreeAvx2> bits;

    Int sum_less_rank(Int l, Int r, Int k) const {
        assert(0 <= l && l <= r && r <= len());
        if (k == 0 || l == r)
            return 0;
        if (k == Int(keys.size()))
            return range_sum(l, r);
        Int out = 0;
        for (Int h = Int(bits.size()) - 1; h >= 0; --h) {
            auto [l0, r0, l1, r1] = matrix.get_child(h, l, r);
            if ((UInt(k) >> h) & 1) {
                out += bits[h].get(l0, r0);
                l = l1;
                r = r1;
            } else {
                l = l0;
                r = r0;
            }
        }
        return out;
    }

public:
    // 座標圧縮した値と、各段のAVX2 Fenwick木。O(N log N+NH)で構築。
    explicit WaveletMatrixFenwick(std::span<const std::pair<Int, Int>> values) {
        weights.resize(values.size());
        for (Int i = 0; i < Int(values.size()); ++i) {
            keys.push_back(values[i].first);
            weights[i] = values[i].second;
        }
        std::sort(keys.begin(), keys.end());
        keys.erase(std::unique(keys.begin(), keys.end()), keys.end());
        Int height = keys.size() <= 1 ? 1 : std::bit_width(keys.size() - 1);
        std::vector<Int> codes(values.size());
        for (Int i = 0; i < Int(values.size()); ++i)
            codes[i] = std::lower_bound(keys.begin(), keys.end(), values[i].first) - keys.begin();
        matrix = initWaveletMatrix(codes, height);
        bits.resize(height);
        auto w = weights;
        std::vector<Int> next(values.size());
        for (Int h = height - 1; h >= 0; --h) {
            for (Int i = 0; i < Int(values.size()); ++i) {
                auto [l0, r0, l1, r1] = matrix.get_child(h, i, i + 1);
                Int p = l0 < r0 ? l0 : l1;
                next[p] = w[i];
            }
            bits[h] = initFenwickTreeAvx2(next);
            w.swap(next);
        }
    }

    Int len() const {
        return weights.size();
    }

    Int operator[](Int i) const {
        assert(0 <= i && i < len());
        return weights[i];
    }

    // 更新と値域付き区間和はO(H log N)。
    void add(Int i, Int delta) {
        assert(0 <= i && i < len());
        weights[i] += delta;
        Int p = i;
        for (Int h = Int(bits.size()) - 1; h >= 0; --h) {
            auto [l0, r0, l1, r1] = matrix.get_child(h, p, p + 1);
            p = l0 < r0 ? l0 : l1;
            bits[h].add(p, delta);
        }
    }

    void set(Int i, Int value) {
        add(i, value - std::as_const(*this)[i]);
    }

    struct Reference {
        WaveletMatrixFenwick *owner;
        Int index;

        operator Int() const {
            return std::as_const(*owner)[index];
        }

        Reference &operator=(Int value) {
            owner->set(index, value);
            return *this;
        }

        Reference &operator=(const Reference &r) {
            return *this = Int(r);
        }
    };

    Reference operator[](Int i) {
        return {this, i};
    }

    CPLIB_BACKWARDS_INDEX_OVERLOADS

    // l <= i < rを満たすb_iの総和をO(log N)で返します。
    Int range_sum(Int l, Int r) const {
        assert(0 <= l && l <= r && r <= len());
        Int h = bits.size() - 1;
        auto [l0, r0, l1, r1] = matrix.get_child(h, l, r);
        return bits[h].get(l0, r0) + bits[h].get(l1, r1);
    }

    // l <= i < rかつa_i <= xを満たすb_iの総和をO(H log N)で返します。
    Int range_sum(Int l, Int r, Int x) const {
        return sum_less_rank(l, r, std::upper_bound(keys.begin(), keys.end(), x) - keys.begin());
    }

    // l <= i < rかつlower <= a_i < upperの重み和をO(H log N)で返します。
    Int range_sum(Int l, Int r, Int lower, Int upper) const {
        assert(lower <= upper);
        return sum_less_rank(l, r,
                             std::lower_bound(keys.begin(), keys.end(), upper) - keys.begin()) -
               sum_less_rank(l, r,
                             std::lower_bound(keys.begin(), keys.end(), lower) - keys.begin());
    }
};

inline auto initWaveletMatrixFenwick(std::span<const std::pair<Int, Int>> values) {
    return WaveletMatrixFenwick(values);
}

inline Int len(const WaveletMatrixFenwick &w) {
    return w.len();
}

inline void add(WaveletMatrixFenwick &w, Int i, Int delta) {
    w.add(i, delta);
}

inline Int range_sum(const WaveletMatrixFenwick &w, Int l, Int r) {
    return w.range_sum(l, r);
}

inline Int range_sum(const WaveletMatrixFenwick &w, Int l, Int r, Int x) {
    return w.range_sum(l, r, x);
}

inline Int range_sum(const WaveletMatrixFenwick &w, Int l, Int r, Int low, Int high) {
    return w.range_sum(l, r, low, high);
}
}
