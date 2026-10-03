#pragma once
#include <cplib/utils/backwards_index.hpp>
#include <memory>
#include <random>

namespace cplib {
inline constexpr UInt RH_MOD = (UInt(1) << 61) - 1;

inline UInt calc_mod(UInt x) {
    x = (x >> 61) + (x & RH_MOD);
    if (x >= RH_MOD)
        x -= RH_MOD;
    return x;
}

inline UInt mul(UInt a, UInt b) {
    UInt au = a >> 31, al = a & ((UInt(1) << 31) - 1), bu = b >> 31, bl = b & ((UInt(1) << 31) - 1),
         mid = al * bu + au * bl;
    return au * bu * 2 + (mid >> 30) + ((mid & ((UInt(1) << 30) - 1)) << 31) + al * bl;
}

inline UInt inner_pow(UInt a, Int n) {
    UInt result = 1;
    for (; n > 0; n >>= 1) {
        if (n & 1)
            result = calc_mod(mul(result, a));
        a = calc_mod(mul(a, a));
    }
    return result;
}

struct Config {
    UInt basei, basej, inversei, inversej;
    std::vector<UInt> powsi, invpowsi, powsj, invpowsj;

    Config() : powsi(500001), invpowsi(500001), powsj(500001), invpowsj(500001) {
        std::mt19937_64 rng(std::random_device{}());
        std::uniform_int_distribution<UInt> distribution(129, UInt(1) << 30);
        basei = distribution(rng);
        inversei = inner_pow(basei, RH_MOD - 2);
        basej = distribution(rng);
        inversej = inner_pow(basej, RH_MOD - 2);
        powsi[0] = invpowsi[0] = powsj[0] = invpowsj[0] = 1;
        for (Int i = 1; i <= 500000; ++i) {
            powsi[i] = calc_mod(mul(powsi[i - 1], basei));
            invpowsi[i] = calc_mod(mul(invpowsi[i - 1], inversei));
            powsj[i] = calc_mod(mul(powsj[i - 1], basej));
            invpowsj[i] = calc_mod(mul(invpowsj[i - 1], inversej));
        }
    }
};

inline Config config;

inline UInt base_powi(Int n) {
    assert(n >= 0);
    return n <= 500000 ? config.powsi[n] : inner_pow(config.basei, n);
}

inline UInt base_powj(Int n) {
    assert(n >= 0);
    return n <= 500000 ? config.powsj[n] : inner_pow(config.basej, n);
}

inline UInt inv_base_powi(Int n) {
    assert(n >= 0);
    return n <= 500000 ? config.invpowsi[n] : inner_pow(config.inversei, n);
}

inline UInt inv_base_powj(Int n) {
    assert(n >= 0);
    return n <= 500000 ? config.invpowsj[n] : inner_pow(config.inversej, n);
}
}

namespace cplib {
template <class T> class HashMatrix;

template <class T> struct HashMatrixRow {
    HashMatrix<T> matrix;
    Int row;

    T operator[](Int col) const {
        return matrix(row, col);
    }
};

// 2方向の累乗で重みを付けた二次元累積和。区間ビューの作成・通常のハッシュ取得はO(1)。
template <class T> class HashMatrix {
    struct Base {
        std::vector<std::vector<T>> matrix;
        std::vector<std::vector<UInt>> hash;
    };

    std::shared_ptr<const Base> base_;
    Int i_ = 0, j_ = 0, height_ = 0, width_ = 0;

public:
    HashMatrix() = default;

    explicit HashMatrix(const std::vector<std::vector<T>> &values)
        : height_(values.size()), width_(values.empty() ? 0 : values[0].size()) {
        auto base = std::make_shared<Base>();
        base->matrix = values;
        base->hash.assign(height_ + 1, std::vector<UInt>(width_ + 1));
        auto &hash = base->hash;
        for (Int i = 0; i < height_; ++i)
            for (Int j = 0; j < width_; ++j) {
                UInt value;
                if constexpr (std::is_same_v<T, char>)
                    value = static_cast<unsigned char>(values[i][j]);
                else
                    value = UInt(values[i][j]);
                hash[i + 1][j + 1] =
                    calc_mod(mul(calc_mod(mul(value, base_powi(i))), base_powj(j)));
            }
        for (Int i = 0; i < height_; ++i)
            for (Int j = 0; j < width_; ++j)
                hash[i + 1][j + 1] = calc_mod(
                    hash[i + 1][j + 1] +
                    calc_mod(calc_mod(hash[i][j + 1] + hash[i + 1][j]) + RH_MOD - hash[i][j]));
        base_ = std::move(base);
    }

    Int h() const {
        return height_;
    }

    Int w() const {
        return width_;
    }

    HashMatrix submatrix(Int top, Int bottom, Int left, Int right) const {
        assert(0 <= top && top <= bottom && bottom <= h() && 0 <= left && left <= right &&
               right <= w());
        auto out = *this;
        out.i_ += top;
        out.j_ += left;
        out.height_ = bottom - top;
        out.width_ = right - left;
        return out;
    }

    template <class I, class J, class K, class L>
    HashMatrix get(ClosedSlice<I, J> rows, ClosedSlice<K, L> cols) const {
        return submatrix(resolve_index(h(), rows.a), resolve_index(h(), rows.b) + 1,
                         resolve_index(w(), cols.a), resolve_index(w(), cols.b) + 1);
    }

    UInt gethash() const {
        const auto &hash = base_->hash;
        UInt value = calc_mod(
            calc_mod(calc_mod(hash[i_ + height_][j_ + width_] + RH_MOD - hash[i_][j_ + width_]) +
                     RH_MOD - hash[i_ + height_][j_]) +
            hash[i_][j_]);
        return calc_mod(mul(value, calc_mod(mul(inv_base_powi(i_), inv_base_powj(j_)))));
    }

    // 元の版と同じく形状は比較せず、正規化したハッシュだけを比較する。
    bool operator==(const HashMatrix &other) const {
        return gethash() == other.gethash();
    }

    T operator()(Int i, Int j) const {
        return base_->matrix[i_ + i][j_ + j];
    }

    HashMatrixRow<T> operator[](Int i) const {
        return {*this, i};
    }
};

template <class T> auto initHashMatrix(const std::vector<std::vector<T>> &values) {
    return HashMatrix<T>(values);
}

template <class T> UInt gethash(const HashMatrix<T> &matrix) {
    return matrix.gethash();
}

template <class T, class I, class J, class K, class L>
auto get(const HashMatrix<T> &matrix, ClosedSlice<I, J> rows, ClosedSlice<K, L> cols) {
    return matrix.get(rows, cols);
}
}
