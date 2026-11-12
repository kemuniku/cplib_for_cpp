#pragma once
#include <cplib/matrix/matrix_avx2_common.hpp>
#include <cplib/modint/modint.hpp>
#include <span>

namespace cplib {
// 行優先の n×m 行列と m×k 行列を法 p で掛ける。AVX2が必要。
// 保守的な時間計算量 O((n+1)(m+1)(k+1))。
inline std::vector<std::uint32_t> matrixProduct(std::span<const std::uint32_t> a,
                                                std::span<const std::uint32_t> b, Int n, Int m,
                                                Int k, std::uint32_t p = 998244353) {
    assert(p > 0 && p < (1u << 30) && (p & 1));
    assert(n >= 0 && m >= 0 && k >= 0 && n <= INT32_MAX && m <= INT32_MAX && k <= INT32_MAX);
    assert(!n || m <= std::numeric_limits<Int>::max() / n);
    assert(!m || k <= std::numeric_limits<Int>::max() / m);
    assert(!n || k <= std::numeric_limits<Int>::max() / n);
    assert(a.size() == UInt(n * m) && b.size() == UInt(m * k));
    for ([[maybe_unused]] auto x : a)
        assert(x < p);
    for ([[maybe_unused]] auto x : b)
        assert(x < p);
    std::vector<std::uint32_t> out(n * k);
    if (n && m && k && p != 1)
        matrixProductKernel(a.data(), b.data(), out.data(), n, m, k, p);
    return out;
}

// 二次元配列の行列積を返す。計算量は行優先配列を受け取る matrixProduct と同じ。
inline std::vector<std::vector<std::uint32_t>>
matrixProduct(const std::vector<std::vector<std::uint32_t>> &a,
              const std::vector<std::vector<std::uint32_t>> &b, std::uint32_t p = 998244353) {
    Int n = a.size(), m = n ? a[0].size() : 0, k = b.empty() ? 0 : b[0].size();
    assert(m == Int(b.size()));
    std::vector<std::uint32_t> aa, bb;
    for (const auto &row : a) {
        assert(Int(row.size()) == m);
        aa.insert(aa.end(), row.begin(), row.end());
    }
    for (const auto &row : b) {
        assert(Int(row.size()) == k);
        bb.insert(bb.end(), row.begin(), row.end());
    }
    auto cc = matrixProduct(aa, bb, n, m, k, p);
    std::vector<std::vector<std::uint32_t>> out(n, std::vector<std::uint32_t>(k));
    for (Int i = 0; i < n; ++i)
        std::copy_n(cc.begin() + i * k, k, out[i].begin());
    return out;
}

}

namespace cplib {
namespace impl = ::cplib::detail::matrix_avx2;

template <impl::Element T, bool Mutable> class BasicMatrixRow {
    std::shared_ptr<impl::Storage<T>> storage_;
    Int offset_ = 0, length_ = 0;

public:
    using reference = std::conditional_t<Mutable, T &, const T &>;
    using pointer = std::conditional_t<Mutable, T *, const T *>;
    // 空の行ビューを作る。計算量 O(1)。
    BasicMatrixRow() = default;

    // 共有領域から行ビューを作る。計算量 O(1)。
    BasicMatrixRow(std::shared_ptr<impl::Storage<T>> s, Int offset, Int length)
        : storage_(std::move(s)), offset_(offset), length_(length) {
    }

    // 行の要素数を返す。計算量 O(1)。
    Int len() const {
        return length_;
    }

    // 行の要素数を返す。計算量 O(1)。
    Int size() const {
        return len();
    }

    // 指定要素への参照を返す。計算量 O(1)。
    reference operator[](Int j) const {
        assert(j >= 0 && j < len());
        return storage_->pointer[offset_ + j];
    }

    // 末尾から数えた要素への参照を返す。計算量 O(1)。
    reference operator[](BackwardsIndex j) const {
        return (*this)[resolve_index(len(), j)];
    }

    // 行の先頭を返す。計算量 O(1)。
    pointer begin() const {
        return storage_ && storage_->pointer ? storage_->pointer + offset_ : nullptr;
    }

    // 行の終端を返す。計算量 O(1)。
    pointer end() const {
        return len() ? begin() + len() : begin();
    }

    // 行の反復可能なビューを返す。計算量 O(1)。
    auto items() const {
        return *this;
    }

    // 書き込み可能な行ビューを返す。計算量 O(1)。
    auto mitems() const
        requires Mutable
    {
        return *this;
    }

    struct PairIterator {
        BasicMatrixRow row;
        Int i;

        // 現在の列番号と値を返す。計算量 O(1)。
        auto operator*() const {
            return std::pair<Int, T>(i, row[i]);
        }

        // 次の要素に進む。計算量 O(1)。
        PairIterator &operator++() {
            ++i;
            return *this;
        }

        // 位置が異なるか判定する。計算量 O(1)。
        bool operator!=(const PairIterator &b) const {
            return i != b.i;
        }
    };

    struct Pairs {
        BasicMatrixRow row;

        // 先頭位置を返す。計算量 O(1)。
        auto begin() const {
            return PairIterator{row, 0};
        }

        // 終端位置を返す。計算量 O(1)。
        auto end() const {
            return PairIterator{row, row.len()};
        }
    };

    // 列番号と値の組を列挙するビューを返す。計算量 O(1)。
    auto pairs() const {
        return Pairs{*this};
    }

    // 行ビューを独立した配列へコピーする。計算量 O(len())。
    std::vector<T> toSeq() const {
        std::vector<T> out(len());
        std::copy(begin(), end(), out.begin());
        return out;
    }

    // 行の要素を区切り文字で連結する。計算量は出力文字数に比例。
    std::string join(std::string_view sep = "") const {
        return len() ? matrixJoinValues(impl::raw(begin()), len(), impl::modulus<T>(),
                                        MontgomeryModint<T>, sep)
                     : "";
    }

    // 行の要素を出力する。計算量は出力文字数に比例。
    void writeRow(FILE *output = stdout) const {
        matrixWriteRow(impl::raw(begin()), len(), impl::modulus<T>(), MontgomeryModint<T>, output);
    }

    // 行ビューをコピーする。計算量 O(1)。
    BasicMatrixRow(const BasicMatrixRow &) = default;

    // 行の要素をコピーする。計算量 O(len())。
    BasicMatrixRow &operator=(const BasicMatrixRow &b)
        requires Mutable
    {
        assert(len() == b.len());
        std::copy(b.begin(), b.end(), begin());
        return *this;
    }

    // 読み取り専用ビューをコピーする。計算量 O(1)。
    BasicMatrixRow &operator=(const BasicMatrixRow &)
        requires(!Mutable)
    = default;

    // 行の要素をコピーする。計算量 O(len())。
    template <bool M>
    BasicMatrixRow &operator=(const BasicMatrixRow<T, M> &b)
        requires Mutable
    {
        assert(len() == b.len());
        std::copy(b.begin(), b.end(), begin());
        return *this;
    }

    // 配列から行の要素をコピーする。計算量 O(len())。
    BasicMatrixRow &operator=(std::span<const T> b)
        requires Mutable
    {
        assert(len() == Int(b.size()));
        std::copy(b.begin(), b.end(), begin());
        return *this;
    }
};
template <class T> using MatrixRow = BasicMatrixRow<T, false>;
template <class T> using MutableMatrixRow = BasicMatrixRow<T, true>;

template <impl::Element T> class Matrix {
    std::shared_ptr<impl::Storage<T>> storage_;
    Int height_ = 0, width_ = 0;
    std::uint32_t modulus_ = 0;

    void check() const {
        [[maybe_unused]] auto p = impl::modulus<T>();
        assert(!modulus_ || modulus_ == p);
    }

    T *data() {
        return storage_ ? storage_->pointer : nullptr;
    }

    const T *data() const {
        return storage_ ? storage_->pointer : nullptr;
    }

    struct Reduction {
        std::vector<std::uint32_t> values;
        std::vector<int> pivots;
        Int width, rank;
        std::uint32_t determinant;
    };

    Reduction reduce(Int extra, bool reduced, const T *rhs = nullptr, bool identity = false) const {
        check();
        assert(extra >= 0 && extra <= INT32_MAX - w());
        Reduction r;
        r.width = w() + extra;
        r.values.resize(impl::size(h(), r.width));
        r.pivots.resize(std::min(h(), w()));
        fieldPrepareKernel(impl::raw(data()), impl::raw(rhs), r.values.data(), h(), w(), extra,
                           impl::modulus<T>(), MontgomeryModint<T>, identity);
        r.rank = fieldEliminateKernel(r.values.data(), h(), r.width, w(), r.pivots.data(),
                                      r.determinant, impl::modulus<T>(), reduced);
        return r;
    }

public:
    using value_type = T;
    // 空行列を作る。計算量 O(1)。
    Matrix() = default;

    // h×w 行列を指定値で初期化する。計算量 O(hw)。
    Matrix(Int h, Int w, T value = T(0)) : height_(h), width_(w), modulus_(impl::modulus<T>()) {
        storage_ = std::make_shared<impl::Storage<T>>(std::vector<T>(impl::size(h, w), value));
    }

    // 行優先の値を移動して行列を作る。計算量 O(1)。
    Matrix(Int h, Int w, std::vector<T> values)
        : height_(h), width_(w), modulus_(impl::modulus<T>()) {
        assert(Int(values.size()) == impl::size(h, w));
        storage_ = std::make_shared<impl::Storage<T>>(std::move(values));
    }

    // 行優先の剰余列から行列を作る。計算量 O(hw)。
    Matrix(Int h, Int w, std::span<const std::uint32_t> values) : Matrix(h, w) {
        assert(Int(values.size()) == h * w);
        for ([[maybe_unused]] auto v : values)
            assert(v < modulus_);
        matrixConvertValues(values.data(), impl::raw(data()), values.size(), modulus_,
                            MontgomeryModint<T>);
    }

    // 行優先の剰余列から行列を作る。計算量 O(hw)。
    Matrix(Int h, Int w, std::vector<std::uint32_t> values)
        : height_(h), width_(w), modulus_(impl::modulus<T>()) {
        assert(Int(values.size()) == impl::size(h, w));
        for ([[maybe_unused]] auto v : values)
            assert(v < modulus_);
        storage_ = std::make_shared<impl::Storage<T>>(std::move(values));
    }

    // 二次元配列から行列を作る。計算量 O(hw)。
    explicit Matrix(const std::vector<std::vector<T>> &rows)
        : Matrix(rows.size(), rows.empty() ? 0 : rows[0].size()) {
        for (Int i = 0; i < h(); ++i) {
            assert(Int(rows[i].size()) == w());
            std::copy(rows[i].begin(), rows[i].end(), data() + i * w());
        }
    }

    // 値を行または列ベクトルとして移動する。計算量 O(1)。
    explicit Matrix(std::vector<T> values, bool vertical = false)
        : height_(vertical ? Int(values.size()) : 1), width_(vertical ? 1 : Int(values.size())),
          modulus_(impl::modulus<T>()) {
        storage_ = std::make_shared<impl::Storage<T>>(std::move(values));
    }

    // 行列を複製する。計算量 O(hw)。
    Matrix(const Matrix &b) : height_(b.height_), width_(b.width_), modulus_(b.modulus_) {
        if (b.storage_) {
            std::vector<T> copy(h() * w());
            std::copy_n(b.data(), copy.size(), copy.data());
            storage_ = std::make_shared<impl::Storage<T>>(std::move(copy));
        }
    }

    // 行列の領域を移動する。計算量 O(1)。
    Matrix(Matrix &&) = default;

    // 行列を代入する。コピー時 O(hw)、移動時 O(1)。
    Matrix &operator=(Matrix b) {
        swap(b);
        return *this;
    }

    // 行列の領域と寸法を交換する。計算量 O(1)。
    void swap(Matrix &b) noexcept {
        storage_.swap(b.storage_);
        std::swap(height_, b.height_);
        std::swap(width_, b.width_);
        std::swap(modulus_, b.modulus_);
    }

    // 行ビューとも領域を共有しない独立した行列を作る。O(h()*w())。
    Matrix clone() const {
        check();
        return *this;
    }

    // 行数を返す。計算量 O(1)。
    Int h() const {
        return height_;
    }

    // 列数を返す。行数0の場合も列数を保持する。計算量 O(1)。
    Int w() const {
        return width_;
    }

    // 指定要素への参照を返す。計算量 O(1)。
    T &operator()(Int i, Int j) {
        assert(i >= 0 && i < h() && j >= 0 && j < w());
        return data()[i * w() + j];
    }

    // 指定要素の参照を返す。計算量 O(1)。
    const T &operator()(Int i, Int j) const {
        assert(i >= 0 && i < h() && j >= 0 && j < w());
        return data()[i * w() + j];
    }

    // 元の行列領域を書き換えられる行ビューを返す。計算量 O(1)。
    MutableMatrixRow<T> operator[](Int i) {
        assert(i >= 0 && i < h());
        return {storage_, i * w(), w()};
    }

    // 行をコピーせず読み取り専用ビューとして返す。計算量 O(1)。
    MatrixRow<T> operator[](Int i) const {
        assert(i >= 0 && i < h());
        return {storage_, i * w(), w()};
    }

    // 行列積を返す。a が n×m、b が m×k のとき保守的に O((n+1)(m+1)(k+1))。
    friend Matrix operator*(const Matrix &a, const Matrix &b) {
        a.check();
        b.check();
        assert(a.w() == b.h());
        Matrix out(a.h(), b.w());
        if (a.h() && a.w() && b.w() && out.modulus_ != 1) {
            if constexpr (MontgomeryModint<T>)
                matrixProductMontgomeryKernel(impl::raw(a.data()), impl::raw(b.data()),
                                              impl::raw(out.data()), a.h(), a.w(), b.w(),
                                              out.modulus_);
            else
                matrixProductKernel(impl::raw(a.data()), impl::raw(b.data()), impl::raw(out.data()),
                                    a.h(), a.w(), b.w(), out.modulus_);
        }
        return out;
    }

    // 右から行列を掛ける。計算量は行列積と同じ。
    Matrix &operator*=(const Matrix &b) {
        auto product = *this * b;
        swap(product);
        return *this;
    }

    // 要素ごとに加算する。計算量 O(hw)。
    Matrix &operator+=(const Matrix &b) {
        check();
        b.check();
        assert(h() == b.h() && w() == b.w());
        for (Int i = 0; i < h() * w(); ++i)
            data()[i] += b.data()[i];
        return *this;
    }

    // 行列の和を返す。計算量 O(hw)。
    friend Matrix operator+(Matrix a, const Matrix &b) {
        return a += b;
    }

    // 要素ごとに加算する。計算量 O(hw)。
    Matrix &operator+=(T value) {
        check();
        for (Int i = 0; i < h() * w(); ++i)
            data()[i] += value;
        return *this;
    }

    // スカラーと行列の和を返す。計算量 O(hw)。
    friend Matrix operator+(Matrix a, T value) {
        return a += value;
    }

    // 要素ごとに減算する。計算量 O(hw)。
    Matrix &operator-=(const Matrix &b) {
        check();
        b.check();
        assert(h() == b.h() && w() == b.w());
        for (Int i = 0; i < h() * w(); ++i)
            data()[i] -= b.data()[i];
        return *this;
    }

    // 行列の差を返す。計算量 O(hw)。
    friend Matrix operator-(Matrix a, const Matrix &b) {
        return a -= b;
    }

    // 要素ごとに減算する。計算量 O(hw)。
    Matrix &operator-=(T value) {
        check();
        for (Int i = 0; i < h() * w(); ++i)
            data()[i] -= value;
        return *this;
    }

    // 行列からスカラーを引いた値を返す。計算量 O(hw)。
    friend Matrix operator-(Matrix a, T value) {
        return a -= value;
    }

    // 要素ごとにスカラー倍する。計算量 O(hw)。
    Matrix &operator*=(T value) {
        check();
        for (Int i = 0; i < h() * w(); ++i)
            data()[i] *= value;
        return *this;
    }

    // スカラー倍を返す。計算量 O(hw)。
    friend Matrix operator*(Matrix a, T value) {
        return a *= value;
    }

    // スカラー倍を返す。計算量 O(hw)。
    friend Matrix operator*(T value, Matrix a) {
        return a *= value;
    }

    // スカラーと行列の和を返す。計算量 O(hw)。
    friend Matrix operator+(T value, Matrix a) {
        return a += value;
    }

    // スカラーから各要素を引いた行列を返す。計算量 O(hw)。
    friend Matrix operator-(T value, const Matrix &a) {
        Matrix out(a.h(), a.w(), value);
        return out -= a;
    }

    // 各要素の符号を反転する。計算量 O(hw)。
    Matrix operator-() const {
        Matrix out(h(), w());
        return out -= *this;
    }

    // n次単位行列を返す。計算量 O(n²)。
    static Matrix identity(Int n, T one = T(1), T zero = T(0)) {
        Matrix out(n, n, zero);
        for (Int i = 0; i < n; ++i)
            out(i, i) = one;
        return out;
    }

    // 正方行列の非負整数乗を返す。行列積を O(log(exponent+1)) 回行う。
    Matrix pow(Int exponent) const {
        check();
        assert(h() == w() && exponent >= 0);
        if (!exponent)
            return identity(h());
        Matrix base = *this, out;
        bool initialized = false;
        while (exponent) {
            if (exponent & 1) {
                if (initialized)
                    out *= base;
                else {
                    out = base;
                    initialized = true;
                }
            }
            exponent >>= 1;
            if (exponent)
                base *= base;
        }
        return out;
    }

    // 全要素の和を返す。計算量 O(hw)。
    T sum() const {
        check();
        T out = 0;
        for (Int i = 0; i < h() * w(); ++i)
            out += data()[i];
        return out;
    }

    // 行列が等しいか判定する。計算量 O(hw)。
    bool operator==(const Matrix &b) const {
        check();
        b.check();
        if (h() != b.h() || w() != b.w())
            return false;
        for (Int i = 0; i < h() * w(); ++i)
            if (data()[i].val() % T::umod() != b.data()[i].val() % T::umod())
                return false;
        return true;
    }

    // 行列のハッシュ値を返す。計算量 O(hw)。
    std::size_t hash() const {
        check();
        std::size_t out = 0;
        auto add = [&](UInt v) {
            out ^= std::hash<UInt>{}(v) + 0x9e3779b97f4a7c15ULL + (out << 6) + (out >> 2);
        };
        add(h());
        add(w());
        add(T::umod());
        for (Int i = 0; i < h() * w(); ++i)
            add(data()[i].val() % T::umod());
        return out;
    }

    // 行列を文字列で返す。計算量は行数と出力文字数に比例。
    std::string str() const {
        check();
        std::string out;
        for (Int i = 0; i < h(); ++i) {
            if (i)
                out += '\n';
            out += (*this)[i].join(" ");
        }
        return out;
    }

    // AVX2の前進消去で階数を求める。O(h*w*min(h,w))。
    Int rank() const {
        return reduce(0, false).rank;
    }

    // AVX2の前進消去で行列式を求める。空行列は1。O(n^3)。
    T determinant() const {
        assert(h() == w());
        auto r = reduce(0, false);
        return r.rank != h() ? T(0) : T(fieldCanonicalKernel(r.determinant, impl::modulus<T>()));
    }

    // 対称な偶数次行列のhafnianをAVX2の多項式積和で求める。O(n^2*2^(n/2))。
    T hafnian() const {
        check();
        assert(h() == w() && h() % 2 == 0);
        for (Int i = 0; i < h(); ++i)
            for (Int j = 0; j < i; ++j)
                assert((*this)(i, j).val() == (*this)(j, i).val());
        return T(
            fieldHafnianKernel(impl::raw(data()), h(), impl::modulus<T>(), MontgomeryModint<T>));
    }

    // AVX2でAx=rhsを掃き出し、特殊解と核の基底を返す。解がなければstd::nullopt。
    // h行w列に対して O(hw min(h,w)+w²)。
    std::optional<LinearSystemSolution<T>> solveLinearSystem(std::span<const T> rhs) const {
        assert(Int(rhs.size()) == h());
        auto r = reduce(1, true, rhs.data());
        return impl::solution<T>(r, h(), w());
    }

    // AVX2の掃き出しで逆行列を返す。特異行列はstd::nullopt。O(n^3)。
    std::optional<Matrix> inverse() const {
        assert(h() == w());
        auto r = reduce(h(), true, nullptr, true);
        if (r.rank != h())
            return std::nullopt;
        Matrix out(h(), h());
        fieldInverseAdjugateKernel(r.values.data(), impl::raw(out.data()), h(), r.rank,
                                   r.pivots.data(), r.determinant, impl::modulus<T>(),
                                   MontgomeryModint<T>, false);
        return out;
    }

    // AVX2で特異行列を含む余因子行列を返す。O(n^3)。
    Matrix adjugate() const {
        assert(h() == w());
        auto r = reduce(h(), true, nullptr, true);
        Matrix out(h(), h());
        fieldInverseAdjugateKernel(r.values.data(), impl::raw(out.data()), h(), r.rank,
                                   r.pivots.data(), r.determinant, impl::modulus<T>(),
                                   MontgomeryModint<T>, true);
        return out;
    }
};

// 指定した引数で行列を構築する。計算量は対応する Matrix コンストラクタと同じ。
template <class T, class... Args> auto initMatrix(Args &&...args) {
    return Matrix<T>(std::forward<Args>(args)...);
}

// 二次元配列から行列を作る。計算量 O(hw)。
template <class T> auto toMatrix(const std::vector<std::vector<T>> &rows) {
    return Matrix<T>(rows);
}

// 指定した引数で行列を構築する。計算量は対応する Matrix コンストラクタと同じ。
template <class T> auto initMatrix(const std::vector<std::vector<T>> &rows) {
    return Matrix<T>(rows);
}

// 行列積を返す。計算量は a * b と同じ。
template <class T> auto matrixProduct(const Matrix<T> &a, const Matrix<T> &b) {
    return a * b;
}

// 単位行列を返す。計算量 O(n²)。
template <class T> auto identity_matrix(Int n, T one = T(1), T zero = T(0)) {
    return Matrix<T>::identity(n, one, zero);
}

// 行列の複製を返す。計算量は a.clone と同じ。
template <class T> auto clone(const Matrix<T> &a) {
    return a.clone();
}

// 行列の行数を返す。計算量は a.h と同じ。
template <class T> auto h(const Matrix<T> &a) {
    return a.h();
}

// 行列の列数を返す。計算量は a.w と同じ。
template <class T> auto w(const Matrix<T> &a) {
    return a.w();
}

// 行列の階数を返す。計算量は a.rank と同じ。
template <class T> auto rank(const Matrix<T> &a) {
    return a.rank();
}

// 行列の行列式を返す。計算量は a.determinant と同じ。
template <class T> auto determinant(const Matrix<T> &a) {
    return a.determinant();
}

// 行列のhafnianを返す。計算量は a.hafnian と同じ。
template <class T> auto hafnian(const Matrix<T> &a) {
    return a.hafnian();
}

// 行列の逆行列を返す。計算量は a.inverse と同じ。
template <class T> auto inverse(const Matrix<T> &a) {
    return a.inverse();
}

// 行列の余因子行列を返す。計算量は a.adjugate と同じ。
template <class T> auto adjugate(const Matrix<T> &a) {
    return a.adjugate();
}

// 行列の全要素の和を返す。計算量は a.sum と同じ。
template <class T> auto sum(const Matrix<T> &a) {
    return a.sum();
}

// 行列のハッシュ値を返す。計算量は a.hash と同じ。
template <class T> auto hash(const Matrix<T> &a) {
    return a.hash();
}

// 行列の非負整数乗を返す。計算量は a.pow と同じ。
template <class T> auto pow(const Matrix<T> &a, Int n) {
    return a.pow(n);
}

// 行列の連立一次方程式の解を返す。計算量は a.solveLinearSystem と同じ。
template <class T> auto solveLinearSystem(const Matrix<T> &a, const std::vector<T> &b) {
    return a.solveLinearSystem(b);
}
}

namespace std {
template <class T> struct hash<cplib::Matrix<T>> {
    // 行列のハッシュ値を返す。計算量は a.hash() と同じ。
    size_t operator()(const cplib::Matrix<T> &a) const {
        return a.hash();
    }
};
}
