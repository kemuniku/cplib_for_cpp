#pragma once
#include <cplib/math/bigint.hpp>

namespace cplib {
class Factoradic;
using FactoradicDivision = QuotientRemainder<Factoradic>;
using FactoradicSmallDivision = QuotientRemainder<Factoradic, Int>;

class Factoradic {
    bool negative_ = false;
    std::vector<Int> data_;
    static constexpr Int BlockSize = 32;

    void normalize() {
        while (!data_.empty() && !data_.back())
            data_.pop_back();
        if (data_.empty())
            negative_ = false;
    }

    static UInt magnitude(Int x) {
        return x < 0 ? UInt(-(x + 1)) + 1 : UInt(x);
    }

    struct ProductNode {
        Int first, last, left = -1, right = -1;
        BigInt product;
    };

    static Int buildProducts(std::vector<ProductNode> &nodes, Int first, Int last) {
        ProductNode node{first, last, -1, -1, BigInt(1)};
        if (last - first <= BlockSize) {
            for (Int i = first; i < last; ++i)
                node.product *= BigInt(i + 1);
        } else {
            Int mid = first + (last - first) / 2;
            node.left = buildProducts(nodes, first, mid);
            node.right = buildProducts(nodes, mid, last);
            node.product = nodes[node.left].product * nodes[node.right].product;
        }
        Int output_value = nodes.size();
        nodes.push_back(std::move(node));
        return output_value;
    }

    static void restoreBlock(BigInt value, const std::vector<ProductNode> &nodes, Int index,
                             std::vector<Int> &digits) {
        if (value.isZero())
            return;
        const auto &node = nodes[index];
        if (node.left < 0) {
            for (Int i = node.first; i < node.last && !value.isZero(); ++i) {
                auto division = cplib::divmod(value, BigInt(i + 1));
                digits[i] = division.remainder.toInt();
                value = std::move(division.quotient);
            }
        } else {
            auto division = cplib::divmod(value, nodes[node.left].product);
            restoreBlock(std::move(division.remainder), nodes, node.left, digits);
            restoreBlock(std::move(division.quotient), nodes, node.right, digits);
        }
    }

    static std::pair<BigInt, BigInt> blockValue(const std::vector<Int> &digits, Int first,
                                                Int last) {
        if (last - first <= BlockSize) {
            BigInt value = 0, product = 1;
            for (Int i = last - 1; i >= first; --i) {
                value = value * BigInt(i + 1) + BigInt(digits[i]);
                product *= BigInt(i + 1);
            }
            return {std::move(value), std::move(product)};
        }
        Int mid = first + (last - first) / 2;
        auto left = blockValue(digits, first, mid), right = blockValue(digits, mid, last);
        return {left.first + left.second * right.first, left.second * right.second};
    }

    static int cmpAbs(const Factoradic &a, const Factoradic &b) {
        if (a.data_.size() != b.data_.size())
            return a.data_.size() < b.data_.size() ? -1 : 1;
        for (Int i = Int(a.data_.size()) - 1; i >= 0; --i)
            if (a.data_[i] != b.data_[i])
                return a.data_[i] < b.data_[i] ? -1 : 1;
        return 0;
    }

    static Factoradic addAbs(const Factoradic &a, const Factoradic &b) {
        Factoradic output_value;
        Int n = std::max(a.data_.size(), b.data_.size());
        output_value.data_.resize(n);
        Int carry = 0;
        for (Int i = 0; i < n; ++i) {
            Int digit = carry;
            if (i < Int(a.data_.size()))
                digit += a.data_[i];
            if (i < Int(b.data_.size()))
                digit += b.data_[i];
            output_value.data_[i] = digit % (i + 1);
            carry = digit / (i + 1);
        }
        if (carry)
            output_value.data_.push_back(carry);
        return output_value;
    }

    static Factoradic subAbs(const Factoradic &a, const Factoradic &b) {
        Factoradic output_value;
        Int n = std::max(a.data_.size(), b.data_.size());
        output_value.data_.resize(n);
        Int borrow = 0;
        for (Int i = 0; i < n; ++i) {
            Int digit = -borrow;
            if (i < Int(a.data_.size()))
                digit += a.data_[i];
            if (i < Int(b.data_.size()))
                digit -= b.data_[i];
            borrow = digit < 0;
            if (borrow)
                digit += i + 1;
            output_value.data_[i] = digit;
        }
        assert(!borrow);
        output_value.normalize();
        return output_value;
    }

    bool fitsInt() const {
        return *this >= Factoradic(std::numeric_limits<Int>::min()) &&
               *this <= Factoradic(std::numeric_limits<Int>::max());
    }

    static std::vector<Int> counts(Int n) {
        std::vector<Int> output_value(n + 1);
        for (Int i = 1; i <= n; ++i)
            output_value[i] = i & -i;
        return output_value;
    }

    static void removeValue(std::vector<Int> &counts, Int value) {
        for (Int i = value + 1; i < Int(counts.size()); i += i & -i)
            --counts[i];
    }

public:
    // 以下 N は階乗進数の桁数、C(N) = (N+1)² log⁴(N+2) は多倍長整数との変換の保守的上界。
    // 零を構築する。計算量 O(1)。
    Factoradic() = default;

    // 絶対値の下位桁から並べたi!の係数digitsと符号で構築する。
    // 0 <= digits[i] <= iとする。O(桁数)。空配列は0。
    explicit Factoradic(std::span<const Int> digits, bool negative = false)
        : negative_(negative), data_(digits.begin(), digits.end()) {
        for (Int i = 0; i < Int(data_.size()); ++i)
            assert(0 <= data_[i] && data_[i] <= i);
        normalize();
    }

    // 固定幅整数を階乗進数に変換する。計算量 O(1)。
    template <std::integral T> Factoradic(T value) {
        UInt m;
        if constexpr (std::is_signed_v<T>) {
            if (value < 0) {
                negative_ = true;
                m = UInt(-(value + 1)) + 1;
            } else
                m = UInt(value);
        } else
            m = UInt(value);
        for (UInt radix = 1; m; ++radix) {
            data_.push_back(m % radix);
            m /= radix;
        }
    }

    // 多倍長整数を積の木と分割統治で階乗進数に変換する。
    // Nを結果の階乗進数の桁数、Bを多倍長整数の桁数、M(B)をその桁数の乗除算の計算量として O(M(B) log(N+1))。
    explicit Factoradic(const BigInt &input) {
        if (BigInt(std::numeric_limits<Int>::min()) <= input &&
            input <= BigInt(std::numeric_limits<Int>::max())) {
            *this = Factoradic(input.toInt());
            return;
        }
        negative_ = input.sgn() < 0;
        auto value = input.abs();
        std::vector<ProductNode> nodes;
        Int last = BlockSize + 1, root = buildProducts(nodes, 1, last);
        while (nodes[root].product <= value) {
            Int next = 1 + 2 * (last - 1), right = buildProducts(nodes, last, next);
            nodes.push_back({1, next, root, right, nodes[root].product * nodes[right].product});
            root = nodes.size() - 1;
            last = next;
        }
        data_.resize(last);
        restoreBlock(std::move(value), nodes, root, data_);
        normalize();
    }

    // 絶対値を返す。計算量 O(N+1)。
    Factoradic abs() const {
        auto output_value = *this;
        output_value.negative_ = false;
        return output_value;
    }

    // 符号反転した値を返す。計算量 O(N+1)。
    Factoradic operator-() const {
        auto output_value = *this;
        if (!output_value.data_.empty())
            output_value.negative_ = !output_value.negative_;
        return output_value;
    }

    // 値を複製する。計算量 O(N+1)。
    Factoradic operator+() const {
        return *this;
    }

    // 符号を -1, 0, 1 で返す。計算量 O(1)。
    int sgn() const {
        return data_.empty() ? 0 : negative_ ? -1 : 1;
    }

    // n!を階乗進数で返す。O(n+1)時間・領域。0!=1、nは非負とする。
    static Factoradic factorial(Int n) {
        assert(n >= 0);
        Factoradic output_value;
        Int index = std::max(n, Int(1));
        output_value.data_.resize(index + 1);
        output_value.data_[index] = 1;
        return output_value;
    }

    // k!による非負の余り。非負値O(min(N,k)+1)、負値O(k+1)。
    Factoradic modFactorial(Int k) const {
        assert(k >= 0);
        Factoradic output_value;
        Int size = std::min(Int(data_.size()), k);
        output_value.data_.assign(data_.begin(), data_.begin() + size);
        output_value.normalize();
        if (negative_ && !output_value.data_.empty()) {
            output_value.data_.resize(k);
            Int borrow = 0;
            for (Int i = 0; i < k; ++i) {
                Int digit = -output_value.data_[i] - borrow;
                borrow = digit < 0;
                if (borrow)
                    digit += i + 1;
                output_value.data_[i] = digit;
            }
            output_value.normalize();
        }
        return output_value;
    }

    // 絶対値の下位桁から並べたi!の係数を返す。0は空配列。O(桁数)。
    std::vector<Int> digits() const {
        return data_;
    }

    // 整数型TにO(桁数)で変換する。Tの範囲外や負数の符号なし型への変換は不可。
    template <std::integral T = Int> T toInt() const {
        UInt value = 0, limit = std::numeric_limits<T>::max();
        if constexpr (std::is_signed_v<T>) {
            if (negative_)
                ++limit;
        } else
            assert(!negative_);
        for (Int i = Int(data_.size()) - 1; i >= 1; --i) {
            UInt digit = data_[i];
            assert(digit <= limit && value <= (limit - digit) / UInt(i + 1));
            value = value * UInt(i + 1) + digit;
        }
        if constexpr (std::is_signed_v<T>) {
            if (negative_) {
                if (value == limit)
                    return std::numeric_limits<T>::min();
                return -T(value);
            }
        }
        return T(value);
    }

    // 組み込み整数での余りは128bit積和による線形計算。追加領域O(1)。
    // 0方向に丸めた除算の、被除数と同符号の余りをIntで返す。
    // O(桁数+1)時間・O(1)領域。0除算はstd::domain_error。
    Int mod(Int divisor) const {
        if (!divisor)
            throw std::domain_error("階乗進数の 0 除算");
        UInt m = magnitude(divisor), rem = 0;
        for (Int i = Int(data_.size()) - 1; i >= 1; --i)
            rem = (static_cast<unsigned __int128>(rem) * UInt(i + 1) + UInt(data_[i])) % m;
        return negative_ ? -Int(rem) : Int(rem);
    }

    // 整数による床除算の余りを返す。計算量 O(N+1)。
    Int operator%(Int divisor) const {
        Int rem = mod(divisor);
        return rem && negative_ != (divisor < 0) ? rem + divisor : rem;
    }

    // 多倍長整数に変換する。計算量 O(C(N))。
    BigInt toBigInt() const {
        BigInt output_value;
        if (Int(data_.size()) <= BlockSize + 1) {
            for (Int i = Int(data_.size()) - 1; i >= 1; --i)
                output_value = output_value * BigInt(i + 1) + BigInt(data_[i]);
        } else
            output_value = blockValue(data_, 1, data_.size()).first;
        return negative_ ? -output_value : output_value;
    }

    // 十進文字列を返す。計算量 O(C(N))。
    std::string str() const {
        return toBigInt().str();
    }

    // 大小を比較して -1, 0, 1 を返す。最大桁数 N に対し O(N+1)。
    friend int cmp(const Factoradic &a, const Factoradic &b) {
        if (a.negative_ != b.negative_)
            return a.negative_ ? -1 : 1;
        Int output_value = cmpAbs(a, b);
        return a.negative_ ? -output_value : output_value;
    }

    // 値が等しいか判定する。最大桁数 N に対し O(N+1)。
    friend bool operator==(const Factoradic &a, const Factoradic &b) {
        return a.negative_ == b.negative_ && a.data_ == b.data_;
    }

    // 大小を比較する。最大桁数 N に対し O(N+1)。
    friend std::strong_ordering operator<=>(const Factoradic &a, const Factoradic &b) {
        return cmp(a, b) <=> 0;
    }

    // ハッシュ値を返す。計算量 O(N+1)。
    std::size_t hash() const {
        std::size_t output_value = negative_;
        for (auto d : data_)
            output_value ^= std::hash<Int>{}(d) + 0x9e3779b97f4a7c15ULL + (output_value << 6) + (output_value >> 2);
        return output_value;
    }

    // 和を返す。最大桁数 N に対し O(N+1)。
    friend Factoradic operator+(const Factoradic &a, const Factoradic &b) {
        Factoradic output_value;
        if (a.negative_ == b.negative_) {
            output_value = addAbs(a, b);
            output_value.negative_ = a.negative_ && !output_value.data_.empty();
        } else if (cmpAbs(a, b) >= 0) {
            output_value = subAbs(a, b);
            output_value.negative_ = a.negative_ && !output_value.data_.empty();
        } else {
            output_value = subAbs(b, a);
            output_value.negative_ = b.negative_ && !output_value.data_.empty();
        }
        return output_value;
    }

    // 差を返す。最大桁数 N に対し O(N+1)。
    friend Factoradic operator-(const Factoradic &a, const Factoradic &b) {
        return a + (-b);
    }

    // 固定幅整数との積を返す。計算量 O(N+1)。
    friend Factoradic operator*(const Factoradic &a, Int b) {
        if (!b || a.data_.empty())
            return {};
        Factoradic output_value;
        output_value.negative_ = a.negative_ != (b < 0);
        UInt m = magnitude(b), carry = 0;
        output_value.data_.resize(a.data_.size());
        for (std::size_t i = 0; i < a.data_.size(); ++i) {
            auto value = static_cast<unsigned __int128>(a.data_[i]) * m + carry;
            output_value.data_[i] = Int(value % (i + 1));
            carry = value / (i + 1);
        }
        while (carry) {
            UInt radix = output_value.data_.size() + 1;
            output_value.data_.push_back(carry % radix);
            carry /= radix;
        }
        return output_value;
    }

    // 固定幅整数との積を返す。計算量 O(N+1)。
    friend Factoradic operator*(Int a, const Factoradic &b) {
        return b * a;
    }

    // 二つの階乗進数の積を返す。最大桁数 N に対し O(C(N))。
    friend Factoradic operator*(const Factoradic &a, const Factoradic &b) {
        if (a.data_.empty() || b.data_.empty())
            return {};
        if (a.fitsInt())
            return b * a.toInt();
        if (b.fitsInt())
            return a * b.toInt();
        return Factoradic(a.toBigInt() * b.toBigInt());
    }

    // 整数による0方向の商と余りを返す。計算量 O(N+1)。
    // 整数による床除算の商と余りを返す。計算量 O(N+1)。
    static FactoradicSmallDivision divmodTrunc(const Factoradic &, Int);
    static FactoradicSmallDivision divmod(const Factoradic &, Int);
    // 階乗進数による床除算の商と余りを返す。最大桁数 N に対し O(C(N))。
    static FactoradicDivision divmod(const Factoradic &, const Factoradic &);
    // 階乗進数による0方向の商を返す。最大桁数 N に対し O(C(N))。
    friend Factoradic operator/(const Factoradic &, const Factoradic &);
    // 整数による0方向の商を返す。計算量 O(N+1)。
    friend Factoradic operator/(const Factoradic &, Int);

    // 整数を階乗進数で割り0方向の商を返す。計算量 O(C(N))。
    friend Factoradic operator/(Int a, const Factoradic &b) {
        return Factoradic(a) / b;
    }

    // 階乗進数による床除算の余りを返す。最大桁数 N に対し O(C(N))。
    friend Factoradic operator%(const Factoradic &, const Factoradic &);

    // 整数を階乗進数で割った床除算の余りを返す。計算量 O(C(N))。
    friend Factoradic operator%(Int a, const Factoradic &b) {
        return Factoradic(a) % b;
    }

    // 階乗進数による0方向の余りを返す。最大桁数 N に対し O(C(N))。
    static Factoradic mod(const Factoradic &, const Factoradic &);
    // 階乗進数による床除算の商を返す。最大桁数 N に対し O(C(N))。
    static Factoradic floorDiv(const Factoradic &, const Factoradic &);

    // 加算の結果を代入する。計算量は対応する二項演算と同じ。
    Factoradic &operator+=(const Factoradic &b) {
        return *this = *this + b;
    }

    // 減算の結果を代入する。計算量は対応する二項演算と同じ。
    Factoradic &operator-=(const Factoradic &b) {
        return *this = *this - b;
    }

    // 乗算の結果を代入する。計算量は対応する二項演算と同じ。
    Factoradic &operator*=(const Factoradic &b) {
        return *this = *this * b;
    }

    // 乗算の結果を代入する。計算量は対応する二項演算と同じ。
    Factoradic &operator*=(Int b) {
        return *this = *this * b;
    }

    // 除算の結果を代入する。計算量は対応する二項演算と同じ。
    Factoradic &operator/=(const Factoradic &b) {
        return *this = *this / b;
    }

    // 除算の結果を代入する。計算量は対応する二項演算と同じ。
    Factoradic &operator/=(Int b) {
        return *this = *this / b;
    }

    // 余り計算の結果を代入する。計算量は対応する二項演算と同じ。
    Factoradic &operator%=(const Factoradic &b) {
        return *this = *this % b;
    }

    // 余り計算の結果を代入する。計算量は対応する二項演算と同じ。
    Factoradic &operator%=(Int b) {
        return *this = Factoradic(*this % b);
    }

    // Fenwick木の個数でLehmer符号を求める。O(N log N)。
    // 0..<Nの順列の0始まりの辞書順順位を階乗進数で返す。O(N log N)時間・O(N)領域。
    static Factoradic permutationRank(std::span<const Int> permutation) {
        Int n = permutation.size();
        auto bit = counts(n);
        std::vector<bool> used(n);
        Factoradic output_value;
        output_value.data_.resize(n);
        for (Int i = 0; i < n; ++i) {
            Int value = permutation[i];
            assert(value >= 0 && value < n && !used[value]);
            used[value] = true;
            for (Int j = value; j > 0; j -= j & -j)
                output_value.data_[n - 1 - i] += bit[j];
            removeValue(bit, value);
        }
        output_value.normalize();
        return output_value;
    }

    // 辞書順でこの値を0始まりの順位とする0..<nの順列を復元する。O(n log n)時間・O(n)領域。
    // n >= 0かつ0 <= この値 < n!を必要とする。n=0の順位0は空順列。
    std::vector<Int> toPermutation(Int n) const {
        assert(n >= 0 && !negative_ && Int(data_.size()) <= n);
        auto bit = counts(n);
        std::vector<Int> output_value(n);
        Int step = 1;
        while (step <= n / 2)
            step *= 2;
        for (Int i = 0; i < n; ++i) {
            Int index = n - 1 - i, rank = index < Int(data_.size()) ? data_[index] : 0, value = 0;
            for (Int width = step; width > 0; width >>= 1) {
                Int next = value + width;
                if (next <= n && bit[next] <= rank) {
                    rank -= bit[next];
                    value = next;
                }
            }
            output_value[i] = value;
            removeValue(bit, value);
        }
        return output_value;
    }
};

// 0方向に丸めた商とIntの余りをO(入力の桁数+1)時間・領域で返す。
// 余りは被除数と同符号。0除算はstd::domain_error。
inline FactoradicSmallDivision Factoradic::divmodTrunc(const Factoradic &a, Int b) {
    if (!b)
        throw std::domain_error("階乗進数の 0 除算");
    FactoradicSmallDivision output_value;
    output_value.quotient.data_.resize(a.data_.size());
    UInt m = magnitude(b), rem = 0;
    for (Int i = Int(a.data_.size()) - 1; i >= 1; --i) {
        auto value = static_cast<unsigned __int128>(rem) * UInt(i + 1) + UInt(a.data_[i]);
        output_value.quotient.data_[i] = value / m;
        rem = value % m;
    }
    output_value.quotient.normalize();
    if (a.negative_ != (b < 0))
        output_value.quotient = -output_value.quotient;
    output_value.remainder = a.negative_ ? -Int(rem) : Int(rem);
    return output_value;
}

// Pythonと同じ床除算の商とIntの余りをO(入力の桁数+1)時間・領域で返す。
// 余りは除数と同符号。0除算はstd::domain_error。
inline FactoradicSmallDivision Factoradic::divmod(const Factoradic &a, Int b) {
    auto output_value = divmodTrunc(a, b);
    if (output_value.remainder && a.negative_ != (b < 0)) {
        output_value.quotient -= 1;
        output_value.remainder += b;
    }
    return output_value;
}

// Pythonと同じ床除算の商と余りを返す。0除算はstd::domain_error。
inline FactoradicDivision Factoradic::divmod(const Factoradic &a, const Factoradic &b) {
    if (b.data_.empty())
        throw std::domain_error("階乗進数の 0 除算");
    if (b.fitsInt()) {
        auto output_value = divmod(a, b.toInt());
        return {std::move(output_value.quotient), Factoradic(output_value.remainder)};
    }
    Int order = cmpAbs(a, b);
    if (order < 0) {
        if (!a.data_.empty() && a.negative_ != b.negative_)
            return {Factoradic(-1), a + b};
        return {Factoradic{}, a};
    }
    if (order == 0)
        return {Factoradic(a.negative_ != b.negative_ ? -1 : 1), Factoradic{}};
    auto output_value = cplib::divmod(a.toBigInt(), b.toBigInt());
    return {Factoradic(output_value.quotient), Factoradic(output_value.remainder)};
}

inline Factoradic operator/(const Factoradic &a, Int b) {
    return Factoradic::divmodTrunc(a, b).quotient;
}

inline Factoradic operator/(const Factoradic &a, const Factoradic &b) {
    if (b.data_.empty())
        throw std::domain_error("階乗進数の 0 除算");
    if (b.fitsInt())
        return Factoradic::divmodTrunc(a, b.toInt()).quotient;
    Int order = Factoradic::cmpAbs(a, b);
    if (order < 0)
        return {};
    if (order == 0)
        return Factoradic(a.negative_ != b.negative_ ? -1 : 1);
    return Factoradic(a.toBigInt() / b.toBigInt());
}

inline Factoradic Factoradic::mod(const Factoradic &a, const Factoradic &b) {
    if (b.data_.empty())
        throw std::domain_error("階乗進数の 0 除算");
    if (b.fitsInt())
        return Factoradic(a.mod(b.toInt()));
    Int order = cmpAbs(a, b);
    if (order < 0)
        return a;
    if (order == 0)
        return {};
    return Factoradic(cplib::mod(a.toBigInt(), b.toBigInt()));
}

inline Factoradic Factoradic::floorDiv(const Factoradic &a, const Factoradic &b) {
    if (b.data_.empty())
        throw std::domain_error("階乗進数の 0 除算");
    if (b.fitsInt())
        return divmod(a, b.toInt()).quotient;
    Int order = cmpAbs(a, b);
    if (order < 0)
        return Factoradic(!a.data_.empty() && a.negative_ != b.negative_ ? -1 : 0);
    if (order == 0)
        return Factoradic(a.negative_ != b.negative_ ? -1 : 1);
    return Factoradic(cplib::floor_div(a.toBigInt(), b.toBigInt()));
}

inline Factoradic operator%(const Factoradic &a, const Factoradic &b) {
    if (b.data_.empty())
        throw std::domain_error("階乗進数の 0 除算");
    if (b.fitsInt())
        return Factoradic(a % b.toInt());
    Int order = Factoradic::cmpAbs(a, b);
    if (order < 0) {
        if (!a.data_.empty() && a.negative_ != b.negative_)
            return a + b;
        return a;
    }
    if (order == 0)
        return {};
    return Factoradic(a.toBigInt() % b.toBigInt());
}

// 床除算の商と余りを返す。計算量は Factoradic::divmod(a, b) と同じ。
inline FactoradicDivision divmod(const Factoradic &a, const Factoradic &b) {
    return Factoradic::divmod(a, b);
}

// 整数による床除算の商と余りを返す。計算量 O(N+1)。
inline FactoradicSmallDivision divmod(const Factoradic &a, Int b) {
    return Factoradic::divmod(a, b);
}

// 整数を階乗進数で割った床除算の商と余りを返す。計算量は Factoradic::divmod と同じ。
inline FactoradicDivision divmod(Int a, const Factoradic &b) {
    return Factoradic::divmod(Factoradic(a), b);
}

// 整数による0方向の商と余りを返す。計算量 O(N+1)。
inline FactoradicSmallDivision divmodTrunc(const Factoradic &a, Int b) {
    return Factoradic::divmodTrunc(a, b);
}

// 0方向に丸めた商を返す。計算量は a / b と同じ。
inline Factoradic div(const Factoradic &a, const Factoradic &b) {
    return a / b;
}

// 整数による0方向に丸めた商を返す。計算量 O(N+1)。
inline Factoradic div(const Factoradic &a, Int b) {
    return a / b;
}

// 整数を階乗進数で割った0方向の商を返す。計算量は Factoradic(a) / b と同じ。
inline Factoradic div(Int a, const Factoradic &b) {
    return Factoradic(a) / b;
}

// 0方向に丸めた余りを返す。計算量は Factoradic::mod(a, b) と同じ。
inline Factoradic mod(const Factoradic &a, const Factoradic &b) {
    return Factoradic::mod(a, b);
}

// 整数による0方向に丸めた余りを返す。計算量 O(N+1)。
inline Int mod(const Factoradic &a, Int b) {
    return a.mod(b);
}

// 整数を階乗進数で割った0方向の余りを返す。計算量は Factoradic::mod と同じ。
inline Factoradic mod(Int a, const Factoradic &b) {
    return Factoradic::mod(Factoradic(a), b);
}

// 床除算の商を返す。計算量は Factoradic::floorDiv(a, b) と同じ。
inline Factoradic floor_div(const Factoradic &a, const Factoradic &b) {
    return Factoradic::floorDiv(a, b);
}

// 整数による床除算の商を返す。計算量 O(N+1)。
inline Factoradic floor_div(const Factoradic &a, Int b) {
    return Factoradic::divmod(a, b).quotient;
}

// 整数を階乗進数で割った床除算の商を返す。計算量は Factoradic::floorDiv と同じ。
inline Factoradic floor_div(Int a, const Factoradic &b) {
    return Factoradic::floorDiv(Factoradic(a), b);
}

// 0方向に丸めた商を a に代入する。計算量は a / b と同じ。
inline Factoradic &div_assign(Factoradic &a, const Factoradic &b) {
    return a = a / b;
}

// 整数による0方向に丸めた商を a に代入する。計算量 O(N+1)。
inline Factoradic &div_assign(Factoradic &a, Int b) {
    return a = a / b;
}

// 0方向に丸めた余りを a に代入する。計算量は Factoradic::mod(a, b) と同じ。
inline Factoradic &mod_assign(Factoradic &a, const Factoradic &b) {
    return a = mod(a, b);
}

// 整数による0方向に丸めた余りを a に代入する。計算量 O(N+1)。
inline Factoradic &mod_assign(Factoradic &a, Int b) {
    return a = Factoradic(mod(a, b));
}

// 値または階乗進数の桁列から構築する。計算量は対応する Factoradic のコンストラクタと同じ。
template <class T> Factoradic initFactoradic(const T &x) {
    return Factoradic(x);
}

// 値または階乗進数の桁列から構築する。計算量は対応する Factoradic のコンストラクタと同じ。
inline Factoradic initFactoradic(std::span<const Int> digits, bool negative) {
    return Factoradic(digits, negative);
}

// n! を階乗進数で返す。時間・領域 O(n+1)。
inline Factoradic factorialFactoradic(Int n) {
    return Factoradic::factorial(n);
}

// k! による非負の余りを返す。計算量は a.modFactorial(k) と同じ。
inline Factoradic modFactorial(const Factoradic &a, Int k) {
    return a.modFactorial(k);
}

// 絶対値を返す。計算量は a.abs() と同じ。
inline Factoradic abs(const Factoradic &a) {
    return a.abs();
}

// 符号を返す。計算量 O(1)。
inline int sgn(const Factoradic &a) {
    return a.sgn();
}

// 下位から階乗進数の桁列を返す。計算量は a.digits() と同じ。
inline auto digits(const Factoradic &a) {
    return a.digits();
}

// 固定幅整数に変換する。計算量は a.toInt() と同じ。
template <std::integral T = Int> T toInt(const Factoradic &a) {
    return a.toInt<T>();
}

// 多倍長整数に変換する。計算量は a.toBigInt() と同じ。
inline BigInt toBigInt(const Factoradic &a) {
    return a.toBigInt();
}

// ハッシュ値を返す。計算量は a.hash() と同じ。
inline std::size_t hashValue(const Factoradic &a) {
    return a.hash();
}

// 順列の辞書順順位を返す。順列長 n に対し O(n log(n+1))。
inline Factoradic permutationRank(std::span<const Int> p) {
    return Factoradic::permutationRank(p);
}

// 順位に対応する長さ n の順列を返す。計算量 O(n log(n+1))。
inline auto toPermutation(const Factoradic &a, Int n) {
    return a.toPermutation(n);
}

// 十進文字列として出力する。計算量は a.str() と同じ。
inline std::ostream &operator<<(std::ostream &os, const Factoradic &a) {
    return os << a.str();
}
}

namespace std {
template <> struct hash<cplib::Factoradic> {
    // ハッシュ値を返す。計算量は x.hash() と同じ。
    size_t operator()(const cplib::Factoradic &x) const {
        return x.hash();
    }
};
}
