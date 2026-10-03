#pragma once
#include <cplib/matrix/matrix_avx2_common.hpp>
#include <cplib/modint/modint.hpp>
#include <span>

namespace cplib {
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
    BasicMatrixRow() = default;

    BasicMatrixRow(std::shared_ptr<impl::Storage<T>> s, Int offset, Int length)
        : storage_(std::move(s)), offset_(offset), length_(length) {
    }

    Int len() const {
        return length_;
    }

    Int size() const {
        return len();
    }

    reference operator[](Int j) const {
        assert(j >= 0 && j < len());
        return storage_->pointer[offset_ + j];
    }

    reference operator[](BackwardsIndex j) const {
        return (*this)[resolve_index(len(), j)];
    }

    pointer begin() const {
        return storage_ && storage_->pointer ? storage_->pointer + offset_ : nullptr;
    }

    pointer end() const {
        return len() ? begin() + len() : begin();
    }

    auto items() const {
        return *this;
    }

    auto mitems() const
        requires Mutable
    {
        return *this;
    }

    struct PairIterator {
        BasicMatrixRow row;
        Int i;

        auto operator*() const {
            return std::pair<Int, T>(i, row[i]);
        }

        PairIterator &operator++() {
            ++i;
            return *this;
        }

        bool operator!=(const PairIterator &b) const {
            return i != b.i;
        }
    };

    struct Pairs {
        BasicMatrixRow row;

        auto begin() const {
            return PairIterator{row, 0};
        }

        auto end() const {
            return PairIterator{row, row.len()};
        }
    };

    auto pairs() const {
        return Pairs{*this};
    }

    // 行ビューを独立した配列へコピーする。
    std::vector<T> toSeq() const {
        std::vector<T> out(len());
        std::copy(begin(), end(), out.begin());
        return out;
    }

    std::string join(std::string_view sep = "") const {
        return len() ? matrixJoinValues(impl::raw(begin()), len(), impl::modulus<T>(),
                                        MontgomeryModint<T>, sep)
                     : "";
    }

    void writeRow(FILE *output = stdout) const {
        matrixWriteRow(impl::raw(begin()), len(), impl::modulus<T>(), MontgomeryModint<T>, output);
    }

    BasicMatrixRow(const BasicMatrixRow &) = default;

    BasicMatrixRow &operator=(const BasicMatrixRow &b)
        requires Mutable
    {
        assert(len() == b.len());
        std::copy(b.begin(), b.end(), begin());
        return *this;
    }

    BasicMatrixRow &operator=(const BasicMatrixRow &)
        requires(!Mutable)
    = default;

    template <bool M>
    BasicMatrixRow &operator=(const BasicMatrixRow<T, M> &b)
        requires Mutable
    {
        assert(len() == b.len());
        std::copy(b.begin(), b.end(), begin());
        return *this;
    }

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
    Matrix() = default;

    Matrix(Int h, Int w, T value = T(0)) : height_(h), width_(w), modulus_(impl::modulus<T>()) {
        storage_ = std::make_shared<impl::Storage<T>>(std::vector<T>(impl::size(h, w), value));
    }

    Matrix(Int h, Int w, std::vector<T> values)
        : height_(h), width_(w), modulus_(impl::modulus<T>()) {
        assert(Int(values.size()) == impl::size(h, w));
        storage_ = std::make_shared<impl::Storage<T>>(std::move(values));
    }

    Matrix(Int h, Int w, std::span<const std::uint32_t> values) : Matrix(h, w) {
        assert(Int(values.size()) == h * w);
        for ([[maybe_unused]] auto v : values)
            assert(v < modulus_);
        matrixConvertValues(values.data(), impl::raw(data()), values.size(), modulus_,
                            MontgomeryModint<T>);
    }

    Matrix(Int h, Int w, std::vector<std::uint32_t> values)
        : height_(h), width_(w), modulus_(impl::modulus<T>()) {
        assert(Int(values.size()) == impl::size(h, w));
        for ([[maybe_unused]] auto v : values)
            assert(v < modulus_);
        storage_ = std::make_shared<impl::Storage<T>>(std::move(values));
    }

    explicit Matrix(const std::vector<std::vector<T>> &rows)
        : Matrix(rows.size(), rows.empty() ? 0 : rows[0].size()) {
        for (Int i = 0; i < h(); ++i) {
            assert(Int(rows[i].size()) == w());
            std::copy(rows[i].begin(), rows[i].end(), data() + i * w());
        }
    }

    explicit Matrix(std::vector<T> values, bool vertical = false)
        : height_(vertical ? Int(values.size()) : 1), width_(vertical ? 1 : Int(values.size())),
          modulus_(impl::modulus<T>()) {
        storage_ = std::make_shared<impl::Storage<T>>(std::move(values));
    }

    Matrix(const Matrix &b) : height_(b.height_), width_(b.width_), modulus_(b.modulus_) {
        if (b.storage_) {
            std::vector<T> copy(h() * w());
            std::copy_n(b.data(), copy.size(), copy.data());
            storage_ = std::make_shared<impl::Storage<T>>(std::move(copy));
        }
    }

    Matrix(Matrix &&) = default;

    Matrix &operator=(Matrix b) {
        swap(b);
        return *this;
    }

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

    Int h() const {
        return height_;
    }

    // 列数を返す。行数0の場合も列数を保持する。
    Int w() const {
        return width_;
    }

    T &operator()(Int i, Int j) {
        assert(i >= 0 && i < h() && j >= 0 && j < w());
        return data()[i * w() + j];
    }

    const T &operator()(Int i, Int j) const {
        assert(i >= 0 && i < h() && j >= 0 && j < w());
        return data()[i * w() + j];
    }

    // 元の行列領域を書き換えられる行ビューを返す。
    MutableMatrixRow<T> operator[](Int i) {
        assert(i >= 0 && i < h());
        return {storage_, i * w(), w()};
    }

    // 行をコピーせず読み取り専用ビューとして返す。
    MatrixRow<T> operator[](Int i) const {
        assert(i >= 0 && i < h());
        return {storage_, i * w(), w()};
    }

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

    Matrix &operator*=(const Matrix &b) {
        auto product = *this * b;
        swap(product);
        return *this;
    }

    Matrix &operator+=(const Matrix &b) {
        check();
        b.check();
        assert(h() == b.h() && w() == b.w());
        for (Int i = 0; i < h() * w(); ++i)
            data()[i] += b.data()[i];
        return *this;
    }

    friend Matrix operator+(Matrix a, const Matrix &b) {
        return a += b;
    }

    Matrix &operator+=(T value) {
        check();
        for (Int i = 0; i < h() * w(); ++i)
            data()[i] += value;
        return *this;
    }

    friend Matrix operator+(Matrix a, T value) {
        return a += value;
    }

    Matrix &operator-=(const Matrix &b) {
        check();
        b.check();
        assert(h() == b.h() && w() == b.w());
        for (Int i = 0; i < h() * w(); ++i)
            data()[i] -= b.data()[i];
        return *this;
    }

    friend Matrix operator-(Matrix a, const Matrix &b) {
        return a -= b;
    }

    Matrix &operator-=(T value) {
        check();
        for (Int i = 0; i < h() * w(); ++i)
            data()[i] -= value;
        return *this;
    }

    friend Matrix operator-(Matrix a, T value) {
        return a -= value;
    }

    Matrix &operator*=(T value) {
        check();
        for (Int i = 0; i < h() * w(); ++i)
            data()[i] *= value;
        return *this;
    }

    friend Matrix operator*(Matrix a, T value) {
        return a *= value;
    }

    friend Matrix operator*(T value, Matrix a) {
        return a *= value;
    }

    friend Matrix operator+(T value, Matrix a) {
        return a += value;
    }

    friend Matrix operator-(T value, const Matrix &a) {
        Matrix out(a.h(), a.w(), value);
        return out -= a;
    }

    Matrix operator-() const {
        Matrix out(h(), w());
        return out -= *this;
    }

    static Matrix identity(Int n, T one = T(1), T zero = T(0)) {
        Matrix out(n, n, zero);
        for (Int i = 0; i < n; ++i)
            out(i, i) = one;
        return out;
    }

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

    T sum() const {
        check();
        T out = 0;
        for (Int i = 0; i < h() * w(); ++i)
            out += data()[i];
        return out;
    }

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

template <class T, class... Args> auto initMatrix(Args &&...args) {
    return Matrix<T>(std::forward<Args>(args)...);
}

template <class T> auto toMatrix(const std::vector<std::vector<T>> &rows) {
    return Matrix<T>(rows);
}

template <class T> auto initMatrix(const std::vector<std::vector<T>> &rows) {
    return Matrix<T>(rows);
}

template <class T> auto matrixProduct(const Matrix<T> &a, const Matrix<T> &b) {
    return a * b;
}

template <class T> auto identity_matrix(Int n, T one = T(1), T zero = T(0)) {
    return Matrix<T>::identity(n, one, zero);
}

template <class T> auto clone(const Matrix<T> &a) {
    return a.clone();
}

template <class T> auto h(const Matrix<T> &a) {
    return a.h();
}

template <class T> auto w(const Matrix<T> &a) {
    return a.w();
}

template <class T> auto rank(const Matrix<T> &a) {
    return a.rank();
}

template <class T> auto determinant(const Matrix<T> &a) {
    return a.determinant();
}

template <class T> auto hafnian(const Matrix<T> &a) {
    return a.hafnian();
}

template <class T> auto inverse(const Matrix<T> &a) {
    return a.inverse();
}

template <class T> auto adjugate(const Matrix<T> &a) {
    return a.adjugate();
}

template <class T> auto sum(const Matrix<T> &a) {
    return a.sum();
}

template <class T> auto hash(const Matrix<T> &a) {
    return a.hash();
}

template <class T> auto pow(const Matrix<T> &a, Int n) {
    return a.pow(n);
}

template <class T> auto solveLinearSystem(const Matrix<T> &a, const std::vector<T> &b) {
    return a.solveLinearSystem(b);
}
}

namespace std {
template <class T> struct hash<cplib::Matrix<T>> {
    size_t operator()(const cplib::Matrix<T> &a) const {
        return a.hash();
    }
};
}
