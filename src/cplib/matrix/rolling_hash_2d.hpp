#pragma once
#include <cplib/utils/backwards_index.hpp>
#include <memory>
#include <random>

namespace cplib {
inline constexpr UInt RH_MOD = (UInt(1) << 61) - 1;

// 61ビットの法で剰余を正規化する。計算量 O(1)。
inline UInt calc_mod(UInt x) {
    x = (x >> 61) + (x & RH_MOD);
    if (x >= RH_MOD)
        x -= RH_MOD;
    return x;
}

// ローリングハッシュ用に2数を乗算する。計算量 O(1)。
inline UInt mul(UInt a, UInt b) {
    UInt au = a >> 31, al = a & ((UInt(1) << 31) - 1), bu = b >> 31, bl = b & ((UInt(1) << 31) - 1),
         mid = al * bu + au * bl;
    return au * bu * 2 + (mid >> 30) + ((mid & ((UInt(1) << 30) - 1)) << 31) + al * bl;
}

// 法RH_MODでaの非負整数n乗を求める。計算量 O(log(n+1))。
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

    // 乱数の基数と累乗・逆累乗の表を初期化する。表長をLとして O(L)。
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

// 行方向の基数のn乗を返す。表の範囲内は O(1)、範囲外は O(log(n+1))。
inline UInt base_powi(Int n) {
    assert(n >= 0);
    return n <= 500000 ? config.powsi[n] : inner_pow(config.basei, n);
}

// 列方向の基数のn乗を返す。表の範囲内は O(1)、範囲外は O(log(n+1))。
inline UInt base_powj(Int n) {
    assert(n >= 0);
    return n <= 500000 ? config.powsj[n] : inner_pow(config.basej, n);
}

// 行方向の基数の逆数のn乗を返す。表の範囲内は O(1)、範囲外は O(log(n+1))。
inline UInt inv_base_powi(Int n) {
    assert(n >= 0);
    return n <= 500000 ? config.invpowsi[n] : inner_pow(config.inversei, n);
}

// 列方向の基数の逆数のn乗を返す。表の範囲内は O(1)、範囲外は O(log(n+1))。
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

    // 指定列の値を返す。計算量 O(1)。
    T operator[](Int col) const {
        return matrix(row, col);
    }
};

// 2方向の累乗で重みを付けた二次元累積和。区間ビューの作成はO(1)、ハッシュ取得は累乗表の範囲内でO(1)。
template <class T> class HashMatrix {
    struct Base {
        std::vector<std::vector<T>> matrix;
        std::vector<std::vector<UInt>> hash;
    };

    std::shared_ptr<const Base> base_;
    Int i_ = 0, j_ = 0, height_ = 0, width_ = 0;

public:
    // 空のビューを作る。計算量 O(1)。
    HashMatrix() = default;

    // H行W列をコピーしてハッシュを構築する。座標が累乗表の範囲内なら O(HW+H)、範囲外では各累乗に O(log(座標+1))。
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

    // ビューの行数を返す。計算量 O(1)。
    Int h() const {
        return height_;
    }

    // ビューの列数を返す。計算量 O(1)。
    Int w() const {
        return width_;
    }

    // 半開矩形を表すビューを返す。計算量 O(1)。
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

    // 両端を含む範囲の部分行列ビューを返す。計算量 O(1)。
    template <class I, class J, class K, class L>
    HashMatrix get(ClosedSlice<I, J> rows, ClosedSlice<K, L> cols) const {
        return submatrix(resolve_index(h(), rows.a), resolve_index(h(), rows.b) + 1,
                         resolve_index(w(), cols.a), resolve_index(w(), cols.b) + 1);
    }

    // ビューのハッシュを原点に正規化して返す。開始座標が累乗表の範囲内なら O(1)、範囲外は O(log(座標+1))。
    UInt gethash() const {
        const auto &hash = base_->hash;
        UInt value = calc_mod(
            calc_mod(calc_mod(hash[i_ + height_][j_ + width_] + RH_MOD - hash[i_][j_ + width_]) +
                     RH_MOD - hash[i_ + height_][j_]) +
            hash[i_][j_]);
        return calc_mod(mul(value, calc_mod(mul(inv_base_powi(i_), inv_base_powj(j_)))));
    }

    // 元の版と同じく形状は比較せず、正規化したハッシュだけを比較する。計算量はgethashの2回分。
    bool operator==(const HashMatrix &other) const {
        return gethash() == other.gethash();
    }

    // 指定位置の値を返す。計算量 O(1)。
    T operator()(Int i, Int j) const {
        return base_->matrix[i_ + i][j_ + j];
    }

    // 指定行のビューを返す。計算量 O(1)。
    HashMatrixRow<T> operator[](Int i) const {
        return {*this, i};
    }
};

// H行W列のハッシュ行列を作る。計算量はHashMatrixの構築と同じ。
template <class T> auto initHashMatrix(const std::vector<std::vector<T>> &values) {
    return HashMatrix<T>(values);
}

// ビューのハッシュを返す。計算量はHashMatrix::gethashと同じ。
template <class T> UInt gethash(const HashMatrix<T> &matrix) {
    return matrix.gethash();
}

// 両端を含む範囲の部分行列ビューを返す。計算量 O(1)。
template <class T, class I, class J, class K, class L>
auto get(const HashMatrix<T> &matrix, ClosedSlice<I, J> rows, ClosedSlice<K, L> cols) {
    return matrix.get(rows, cols);
}
}
