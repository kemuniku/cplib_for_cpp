#pragma once
#include <cplib/matrix/matrix_avx2_common.hpp>
#include <array>

namespace cplib::avx2_static_matrix {
namespace impl = ::cplib::detail::matrix_avx2;

template <Int H, Int W, impl::Element T> class StaticMatrix {
    static_assert(H >= 0 && W >= 0 && H <= INT32_MAX && W <= INT32_MAX);
    static_assert(H == 0 || W <= (std::numeric_limits<Int>::max() / 4) / H);
    std::array<T, H * W> values_{};
    std::uint32_t modulus_ = 0;
    template <Int, Int, impl::Element> friend class StaticMatrix;

    void check() const {
        [[maybe_unused]] auto p = impl::modulus<T>();
        assert(!modulus_ || modulus_ == p);
    }

    template <Int E> struct Reduction {
        std::array<std::uint32_t, H *(W + E)> values{};
        std::array<int, std::min(H, W)> pivots{};
        Int width = 0, rank = 0;
        std::uint32_t determinant = 0;
    };

    template <Int E>
    auto reduce(Int height, Int width, bool reduced, const T *rhs = nullptr,
                bool identity = false) const {
        static_assert(E >= 0 && W <= INT32_MAX - E);
        check();
        assert(height >= 0 && height <= H && width >= 0 && width <= W);
        auto r = std::make_unique<Reduction<E>>();
        Int extra = identity ? height : E;
        r->width = width + extra;
        fieldPrepareKernel(impl::raw(values_.data()), impl::raw(rhs), r->values.data(), height,
                           width, extra, impl::modulus<T>(), MontgomeryModint<T>, identity, W);
        r->rank = fieldEliminateKernel(r->values.data(), height, r->width, width, r->pivots.data(),
                                       r->determinant, impl::modulus<T>(), reduced);
        return r;
    }

public:
    using value_type = T;
    // 零行列を作る。計算量 O(HW)。
    StaticMatrix() = default;

    // 全要素を指定値で埋める。計算量 O(HW)。
    explicit StaticMatrix(T value) : modulus_(impl::modulus<T>()) {
        values_.fill(value);
    }

    // 二次元配列から行列を作る。計算量 O(HW)。
    explicit StaticMatrix(const std::array<std::array<T, W>, H> &rows)
        : modulus_(impl::modulus<T>()) {
        for (Int i = 0; i < H; ++i)
            std::copy(rows[i].begin(), rows[i].end(), values_.begin() + i * W);
    }

    // 行数を返す。計算量 O(1)。
    static constexpr Int h() {
        return H;
    }

    // 列数を返す。計算量 O(1)。
    static constexpr Int w() {
        return W;
    }

    // 指定要素への参照を返す。計算量 O(1)。
    const T &operator()(Int i, Int j) const {
        assert(i >= 0 && i < H && j >= 0 && j < W);
        return values_[i * W + j];
    }

    // 指定要素への書き込み可能な参照を返す。計算量 O(1)。
    T &operator()(Int i, Int j) {
        assert(i >= 0 && i < H && j >= 0 && j < W);
        return values_[i * W + j];
    }

    // 指定行をコピーする。O(W)。
    std::array<T, W> operator[](Int i) const {
        assert(i >= 0 && i < H);
        std::array<T, W> output_value{};
        std::copy_n(values_.begin() + i * W, W, output_value.begin());
        return output_value;
    }

    struct Row {
        StaticMatrix *a;
        Int i;

        // 行を配列にコピーする。計算量 O(W)。
        operator std::array<T, W>() const {
            return std::as_const(*a)[i];
        }

        // 行の指定要素への参照を返す。計算量 O(1)。
        T &operator[](Int j) const {
            return (*a)(i, j);
        }

        // 配列から行を代入する。計算量 O(W)。
        Row &operator=(const std::array<T, W> &row) {
            a->setRow(i, row);
            return *this;
        }

        // 別の行を代入する。計算量 O(W)。
        Row &operator=(const Row &b) {
            return *this = std::array<T, W>(b);
        }
    };

    // 書き込み可能な行ビューを返す。計算量 O(1)。
    Row operator[](Int i) {
        assert(i >= 0 && i < H);
        return {this, i};
    }

    // 配列から行を設定する。計算量 O(W)。
    void setRow(Int i, const std::array<T, W> &row) {
        assert(i >= 0 && i < H);
        std::copy(row.begin(), row.end(), values_.begin() + i * W);
    }

    // 行列が等しいか判定する。計算量 O(HW)。
    bool operator==(const StaticMatrix &b) const {
        check();
        b.check();
        for (Int i = 0; i < H * W; ++i)
            if (values_[i].val() != b.values_[i].val())
                return false;
        return true;
    }

    // 行列を文字列で返す。計算量は行数と出力文字数に比例。
    std::string str() const {
        check();
        std::string output_value;
        for (Int i = 0; i < H; ++i) {
            if (i)
                output_value += '\n';
            if constexpr (W > 0)
                output_value += matrixJoinValues(impl::raw(values_.data() + i * W), W, impl::modulus<T>(),
                                        MontgomeryModint<T>, " ");
        }
        return output_value;
    }

    // 行列積を返す。保守的な計算量 O((H+1)(W+1)(K+1))。
    template <Int K> StaticMatrix<H, K, T> operator*(const StaticMatrix<W, K, T> &b) const {
        check();
        b.check();
        StaticMatrix<H, K, T> output_value(T(0));
        if constexpr (H <= 4 && W <= 4 && K <= 4) {
            for (Int i = 0; i < H; ++i)
                for (Int j = 0; j < K; ++j)
                    for (Int k = 0; k < W; ++k)
                        output_value(i, j) += (*this)(i, k) * b(k, j);
        } else if constexpr (H > 0 && W > 0 && K > 0) {
            if constexpr (MontgomeryModint<T>)
                matrixProductMontgomeryKernel(
                    impl::raw(values_.data()), impl::raw(b.values_.data()),
                    impl::raw(output_value.values_.data()), H, W, K, impl::modulus<T>());
            else
                matrixProductKernel(impl::raw(values_.data()), impl::raw(b.values_.data()),
                                    impl::raw(output_value.values_.data()), H, W, K, impl::modulus<T>());
        }
        return output_value;
    }

    // 右から行列を掛ける。計算量は行列積と同じ。
    StaticMatrix &operator*=(const StaticMatrix<W, W, T> &b) {
        return *this = *this * b;
    }

    // 要素ごとに加算する。計算量 O(HW)。
    StaticMatrix &operator+=(const StaticMatrix &b) {
        check();
        b.check();
        for (Int i = 0; i < H * W; ++i)
            values_[i] += b.values_[i];
        return *this;
    }

    // 行列の和を返す。計算量 O(HW)。
    friend StaticMatrix operator+(StaticMatrix a, const StaticMatrix &b) {
        return a += b;
    }

    // 要素ごとに加算する。計算量 O(HW)。
    StaticMatrix &operator+=(T value) {
        check();
        for (auto &x : values_)
            x += value;
        return *this;
    }

    // スカラーと行列の和を返す。計算量 O(HW)。
    friend StaticMatrix operator+(StaticMatrix a, T value) {
        return a += value;
    }

    // 要素ごとに減算する。計算量 O(HW)。
    StaticMatrix &operator-=(const StaticMatrix &b) {
        check();
        b.check();
        for (Int i = 0; i < H * W; ++i)
            values_[i] -= b.values_[i];
        return *this;
    }

    // 行列の差を返す。計算量 O(HW)。
    friend StaticMatrix operator-(StaticMatrix a, const StaticMatrix &b) {
        return a -= b;
    }

    // 要素ごとに減算する。計算量 O(HW)。
    StaticMatrix &operator-=(T value) {
        check();
        for (auto &x : values_)
            x -= value;
        return *this;
    }

    // 行列からスカラーを引いた値を返す。計算量 O(HW)。
    friend StaticMatrix operator-(StaticMatrix a, T value) {
        return a -= value;
    }

    // 要素ごとにスカラー倍する。計算量 O(HW)。
    StaticMatrix &operator*=(T value) {
        check();
        for (auto &x : values_)
            x *= value;
        return *this;
    }

    // スカラー倍を返す。計算量 O(HW)。
    friend StaticMatrix operator*(StaticMatrix a, T value) {
        return a *= value;
    }

    // スカラー倍を返す。計算量 O(HW)。
    friend StaticMatrix operator*(T value, StaticMatrix a) {
        return a *= value;
    }

    // 各要素の符号を反転する。計算量 O(HW)。
    StaticMatrix operator-() const {
        return *this * T(-1);
    }

    // 単位行列を返す。計算量 O(HW)。
    static StaticMatrix identity(Int n = H) {
        static_assert(H == W);
        assert(n == H);
        StaticMatrix output_value(T(0));
        for (Int i = 0; i < H; ++i)
            output_value(i, i) = 1;
        return output_value;
    }

    // 正方行列の非負整数乗を返す。行列積を O(log(exponent+1)) 回行う。
    StaticMatrix pow(Int exponent) const {
        static_assert(H == W);
        assert(exponent >= 0);
        check();
        StaticMatrix output_value = identity(), base = *this;
        while (exponent) {
            if (exponent & 1)
                output_value *= base;
            exponent >>= 1;
            if (exponent)
                base *= base;
        }
        return output_value;
    }

    // 全要素の和を返す。計算量 O(HW)。
    T sum() const {
        check();
        T output_value = 0;
        for (auto x : values_)
            output_value += x;
        return output_value;
    }

    // 左上height×widthの階数をAVX2で求める。O(H*W+h*w*min(h,w))。
    Int rank(Int height = H, Int width = W) const {
        return reduce<0>(height, width, false)->rank;
    }

    // 左上n×nの行列式をAVX2で求める。空行列は1。O(H*W+n^3)。
    T determinant(Int n = H) const {
        auto r = reduce<0>(n, n, false);
        return r->rank != n ? T(0) : T(fieldCanonicalKernel(r->determinant, impl::modulus<T>()));
    }

    // 対称な左上n×nのhafnianをAVX2で求める。O(n^2*2^(n/2))。
    T hafnian(Int n = H) const {
        check();
        assert(n >= 0 && n <= std::min(H, W) && n % 2 == 0);
        for (Int i = 0; i < n; ++i)
            for (Int j = 0; j < i; ++j)
                assert((*this)(i, j).val() == (*this)(j, i).val());
        return T(fieldHafnianKernel(impl::raw(values_.data()), n, impl::modulus<T>(),
                                    MontgomeryModint<T>, W));
    }

    // 左上height×widthでAx=rhsの特殊解と核の基底を求める。解がなければstd::nullopt。
    // h=height、w=widthとして O(HW+hw min(h,w)+w²)。
    std::optional<LinearSystemSolution<T>> solveLinearSystem(std::span<const T> rhs, Int height = H,
                                                             Int width = W) const {
        assert(Int(rhs.size()) == height);
        auto r = reduce<1>(height, width, true, rhs.data());
        return impl::solution<T>(*r, height, width);
    }

    // 左上n×nの逆行列をAVX2で求める。範囲外は零、特異行列はstd::nullopt。O(H*W+n^3)。
    std::optional<StaticMatrix> inverse(Int n = H) const {
        auto r = reduce<std::min(H, W)>(n, n, true, nullptr, true);
        if (r->rank != n)
            return std::nullopt;
        StaticMatrix output_value(T(0));
        fieldInverseAdjugateKernel(r->values.data(), impl::raw(output_value.values_.data()), n, r->rank,
                                   r->pivots.data(), r->determinant, impl::modulus<T>(),
                                   MontgomeryModint<T>, false, W);
        return output_value;
    }

    // 左上n×nの余因子行列をAVX2で求める。範囲外は零。O(H*W+n^3)。
    StaticMatrix adjugate(Int n = H) const {
        auto r = reduce<std::min(H, W)>(n, n, true, nullptr, true);
        StaticMatrix output_value(T(0));
        fieldInverseAdjugateKernel(r->values.data(), impl::raw(output_value.values_.data()), n, r->rank,
                                   r->pivots.data(), r->determinant, impl::modulus<T>(),
                                   MontgomeryModint<T>, true, W);
        return output_value;
    }
};

// 指定値で固定長行列を作る。計算量 O(HW)。
template <Int H, Int W, class T> auto initMatrix(T value = T(0)) {
    return StaticMatrix<H, W, T>(value);
}

// 二次元配列から固定長行列を作る。計算量 O(HW)。
template <std::size_t H, std::size_t W, class T>
auto toMatrix(const std::array<std::array<T, W>, H> &rows) {
    return StaticMatrix<H, W, T>(rows);
}

// 単位行列を返す。計算量 O(HW)。
template <Int H, Int W, class T> auto identity_matrix(Int n = H) {
    return StaticMatrix<H, W, T>::identity(n);
}

// 行列積を返す。計算量は a * b と同じ。
template <Int H, Int W, Int K, class T>
auto matrixProduct(const StaticMatrix<H, W, T> &a, const StaticMatrix<W, K, T> &b) {
    return a * b;
}

// 行列の行数を返す。計算量は a.h と同じ。
template <Int H, Int W, class T, class... A> auto h(const StaticMatrix<H, W, T> &a, A &&...args) {
    return a.h(std::forward<A>(args)...);
}

// 行列の列数を返す。計算量は a.w と同じ。
template <Int H, Int W, class T, class... A> auto w(const StaticMatrix<H, W, T> &a, A &&...args) {
    return a.w(std::forward<A>(args)...);
}

// 行列の階数を返す。計算量は a.rank と同じ。
template <Int H, Int W, class T, class... A>
auto rank(const StaticMatrix<H, W, T> &a, A &&...args) {
    return a.rank(std::forward<A>(args)...);
}

// 行列の行列式を返す。計算量は a.determinant と同じ。
template <Int H, Int W, class T, class... A>
auto determinant(const StaticMatrix<H, W, T> &a, A &&...args) {
    return a.determinant(std::forward<A>(args)...);
}

// 行列のhafnianを返す。計算量は a.hafnian と同じ。
template <Int H, Int W, class T, class... A>
auto hafnian(const StaticMatrix<H, W, T> &a, A &&...args) {
    return a.hafnian(std::forward<A>(args)...);
}

// 行列の逆行列を返す。計算量は a.inverse と同じ。
template <Int H, Int W, class T, class... A>
auto inverse(const StaticMatrix<H, W, T> &a, A &&...args) {
    return a.inverse(std::forward<A>(args)...);
}

// 行列の余因子行列を返す。計算量は a.adjugate と同じ。
template <Int H, Int W, class T, class... A>
auto adjugate(const StaticMatrix<H, W, T> &a, A &&...args) {
    return a.adjugate(std::forward<A>(args)...);
}

// 行列の全要素の和を返す。計算量は a.sum と同じ。
template <Int H, Int W, class T, class... A> auto sum(const StaticMatrix<H, W, T> &a, A &&...args) {
    return a.sum(std::forward<A>(args)...);
}

// 行列の非負整数乗を返す。計算量は a.pow と同じ。
template <Int H, Int W, class T, class... A> auto pow(const StaticMatrix<H, W, T> &a, A &&...args) {
    return a.pow(std::forward<A>(args)...);
}

// 行列の連立一次方程式の解を返す。計算量は a.solveLinearSystem と同じ。
template <Int H, Int W, class T, class... A>
auto solveLinearSystem(const StaticMatrix<H, W, T> &a, A &&...args) {
    return a.solveLinearSystem(std::forward<A>(args)...);
}
}
