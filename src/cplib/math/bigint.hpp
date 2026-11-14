#pragma once
#include <cplib/convolution/convolution.hpp>
#include <compare>
#include <memory>
#include <stdexcept>
#include <ostream>

namespace cplib {
class BigInt;

template <class T, class R = T> struct QuotientRemainder {
    T quotient;
    R remainder;
};

using BigIntDivision = QuotientRemainder<BigInt>;

class BigInt {
    int sign_ = 0;
    std::vector<std::uint32_t> digits_;
    static constexpr UInt Base = 1000000000, Mod1 = 754974721, Mod2 = 469762049,
                          InvMod1 = 221064492;
    static constexpr Int NttThreshold = 64, NttMaxLength = Int(1) << 24,
                         LargeBaseMaxDigits = (Mod1 * Mod2 - 1) / (UInt(999999) * 999999),
                         NewtonThreshold = 256;

    void normalize() {
        while (!digits_.empty() && !digits_.back())
            digits_.pop_back();
        if (digits_.empty())
            sign_ = 0;
    }

    static int cmpAbs(const BigInt &a, const BigInt &b) {
        if (a.digits_.size() != b.digits_.size())
            return a.digits_.size() < b.digits_.size() ? -1 : 1;
        for (Int i = Int(a.digits_.size()) - 1; i >= 0; --i)
            if (a.digits_[i] != b.digits_[i])
                return a.digits_[i] < b.digits_[i] ? -1 : 1;
        return 0;
    }

    static BigInt addAbs(const BigInt &a, const BigInt &b) {
        BigInt output_value;
        output_value.sign_ = 1;
        Int n = std::max(a.digits_.size(), b.digits_.size());
        output_value.digits_.resize(n + 1);
        UInt carry = 0;
        for (Int i = 0; i < n; ++i) {
            if (i < Int(a.digits_.size()))
                carry += a.digits_[i];
            if (i < Int(b.digits_.size()))
                carry += b.digits_[i];
            output_value.digits_[i] = carry % Base;
            carry /= Base;
        }
        output_value.digits_[n] = carry;
        output_value.normalize();
        return output_value;
    }

    static BigInt subAbs(const BigInt &a, const BigInt &b) {
        BigInt output_value;
        output_value.sign_ = 1;
        output_value.digits_.resize(a.digits_.size());
        Int borrow = 0;
        for (Int i = 0; i < Int(a.digits_.size()); ++i) {
            Int digit = Int(a.digits_[i]) - borrow;
            if (i < Int(b.digits_.size()))
                digit -= b.digits_[i];
            borrow = 0;
            if (digit < 0) {
                digit += Base;
                borrow = 1;
            }
            output_value.digits_[i] = digit;
        }
        output_value.normalize();
        return output_value;
    }

    std::vector<std::uint32_t> binaryDigits() const {
        std::vector<std::uint32_t> output_value;
        for (Int i = Int(digits_.size()) - 1; i >= 0; --i) {
            UInt carry = digits_[i];
            for (auto &d : output_value) {
                UInt value = UInt(d) * Base + carry;
                d = std::uint32_t(value);
                carry = value >> 32;
            }
            if (carry)
                output_value.push_back(carry);
        }
        return output_value;
    }

    static BigInt fromBinary(const std::vector<std::uint32_t> &digits) {
        BigInt output_value;
        output_value.sign_ = 1;
        for (Int i = Int(digits.size()) - 1; i >= 0; --i) {
            UInt carry = digits[i];
            for (auto &d : output_value.digits_) {
                UInt value = (UInt(d) << 32) + carry;
                d = value % Base;
                carry = value / Base;
            }
            while (carry) {
                output_value.digits_.push_back(carry % Base);
                carry /= Base;
            }
        }
        output_value.normalize();
        return output_value;
    }

    static void complement(std::vector<std::uint32_t> &digits) {
        UInt carry = 1;
        for (auto &d : digits) {
            UInt value = UInt(~d) + carry;
            d = std::uint32_t(value);
            carry = value >> 32;
        }
    }

    template <int Op> static BigInt bitwise(const BigInt &a, const BigInt &b) {
        auto left = a.binaryDigits(), right = b.binaryDigits();
        auto size = std::max(left.size(), right.size()) + 1;
        left.resize(size);
        right.resize(size);
        if (a.sign_ < 0)
            complement(left);
        if (b.sign_ < 0)
            complement(right);
        for (std::size_t i = 0; i < size; ++i) {
            if constexpr (Op == 0)
                left[i] &= right[i];
            else if constexpr (Op == 1)
                left[i] |= right[i];
            else
                left[i] ^= right[i];
        }
        bool negative = left.back() >> 31;
        if (negative)
            complement(left);
        auto output_value = fromBinary(left);
        if (negative)
            output_value.sign_ = -output_value.sign_;
        return output_value;
    }

    // 10^9進の筆算乗算。O(nm)。
    static BigInt mulSchoolbook(const BigInt &a, const BigInt &b) {
        BigInt output_value;
        output_value.sign_ = a.sign_ * b.sign_;
        output_value.digits_.resize(a.digits_.size() + b.digits_.size());
        for (std::size_t i = 0; i < a.digits_.size(); ++i) {
            UInt carry = 0;
            for (std::size_t j = 0; j < b.digits_.size(); ++j) {
                UInt digit = UInt(output_value.digits_[i + j]) + UInt(a.digits_[i]) * b.digits_[j] + carry;
                output_value.digits_[i + j] = digit % Base;
                carry = digit / Base;
            }
            output_value.digits_[i + b.digits_.size()] = carry;
        }
        output_value.normalize();
        return output_value;
    }

    static std::vector<std::uint32_t> convolutionDigits(const std::vector<std::uint32_t> &a,
                                                        const std::vector<std::uint32_t> &b,
                                                        std::uint32_t p, std::uint32_t root) {
        auto length = a.size() + b.size() - 1, size = std::bit_ceil(length);
        std::vector<std::uint32_t> output_value(size);
        detail::avx2_ntt::convolution_ntt_friendly(output_value.data(), a.data(), a.size(), b.data(),
                                                   b.size(), size, p, root, false);
        output_value.resize(length);
        return output_value;
    }

    template <int D> static std::vector<std::uint32_t> splitNttDigits(const BigInt &x) {
        constexpr UInt base = D == 6 ? 1000000 : 100000;
        std::vector<std::uint32_t> output_value((x.digits_.size() * 9 + D - 1) / D);
        UInt value = 0, scale = 1;
        std::size_t position = 0;
        for (auto digit : x.digits_) {
            value += UInt(digit) * scale;
            scale *= Base;
            while (scale >= base) {
                output_value[position++] = value % base;
                value /= base;
                scale /= base;
            }
        }
        if (scale > 1)
            output_value[position] = value;
        return output_value;
    }

    template <int D>
    static BigInt restoreNttProduct(const std::vector<std::uint32_t> &c1,
                                    const std::vector<std::uint32_t> &c2, Int length, int sign) {
        constexpr UInt base = D == 6 ? 1000000 : 100000;
        BigInt output_value;
        output_value.sign_ = sign;
        output_value.digits_.resize(length);
        UInt carry = 0, value = 0, scale = 1;
        Int position = 0;
        auto append = [&]() {
            value += (carry % base) * scale;
            carry /= base;
            scale *= base;
            if (scale >= Base) {
                output_value.digits_[position++] = value % Base;
                value /= Base;
                scale /= Base;
            }
        };
        for (std::size_t i = 0; i < c1.size(); ++i) {
            UInt r1 = c1[i], t2 = ((UInt(c2[i]) + 2 * Mod2 - r1) * InvMod1) % Mod2;
            carry += r1 + Mod1 * t2;
            append();
        }
        while (carry)
            append();
        if (value)
            output_value.digits_[position] = value;
        output_value.normalize();
        return output_value;
    }

    template <int D> static BigInt mulNtt(const BigInt &a, const BigInt &b) {
        auto left = splitNttDigits<D>(a), right = splitNttDigits<D>(b);
        auto c1 = convolutionDigits(left, right, Mod1, 11),
             c2 = convolutionDigits(left, right, Mod2, 3);
        return restoreNttProduct<D>(c1, c2, a.digits_.size() + b.digits_.size(), a.sign_ * b.sign_);
    }

    static BigInt mulAbsSmall(const BigInt &x, std::uint32_t y) {
        BigInt output_value;
        output_value.sign_ = 1;
        output_value.digits_.resize(x.digits_.size() + 1);
        UInt carry = 0;
        for (std::size_t i = 0; i < x.digits_.size(); ++i) {
            UInt digit = UInt(x.digits_[i]) * y + carry;
            output_value.digits_[i] = digit % Base;
            carry = digit / Base;
        }
        output_value.digits_.back() = carry;
        output_value.normalize();
        return output_value;
    }

    static std::pair<BigInt, std::uint32_t> divAbsSmall(const BigInt &x, std::uint32_t y) {
        BigInt output_value;
        output_value.sign_ = 1;
        output_value.digits_.resize(x.digits_.size());
        UInt rem = 0;
        for (Int i = Int(x.digits_.size()) - 1; i >= 0; --i) {
            UInt digit = rem * Base + x.digits_[i];
            output_value.digits_[i] = digit / y;
            rem = digit % y;
        }
        output_value.normalize();
        return {std::move(output_value), std::uint32_t(rem)};
    }

    static BigIntDivision divmodSchoolbook(const BigInt &, const BigInt &);

    struct FixedMultiplier {
        std::unique_ptr<detail::avx2_ntt::FixedConvolution> first, second;
        Int length = 0, fixedLength = 0;
        std::vector<std::uint32_t> c1, c2;
        FixedMultiplier() = default;

        FixedMultiplier(const BigInt &y, Int maxLeftLength) {
            fixedLength = (Int(y.digits_.size()) * 3 + 1) / 2;
            Int leftLength = (maxLeftLength * 3 + 1) / 2;
            if (std::min(fixedLength, leftLength) > LargeBaseMaxDigits ||
                fixedLength + leftLength - 1 > NttMaxLength)
                return;
            length = std::bit_ceil(UInt(fixedLength + leftLength - 1));
            auto digits = splitNttDigits<6>(y);
            first = std::make_unique<detail::avx2_ntt::FixedConvolution>(
                digits.data(), digits.size(), length, Mod1, 11);
            second = std::make_unique<detail::avx2_ntt::FixedConvolution>(
                digits.data(), digits.size(), length, Mod2, 3);
            c1.resize(length);
            c2.resize(length);
        }

        BigInt multiply(const BigInt &x, const BigInt &y) {
            if (!first || Int(x.digits_.size()) <= NttThreshold)
                return x * y;
            auto left = splitNttDigits<6>(x);
            c1.resize(length);
            c2.resize(length);
            first->run(c1.data(), left.data(), left.size());
            second->run(c2.data(), left.data(), left.size());
            c1.resize(left.size() + fixedLength - 1);
            c2.resize(c1.size());
            return restoreNttProduct<6>(c1, c2, x.digits_.size() + y.digits_.size(),
                                        x.sign_ * y.sign_);
        }
    };

    static BigInt shiftLeft(const BigInt &x, Int count) {
        if (x.isZero())
            return {};
        BigInt output_value;
        output_value.sign_ = x.sign_;
        output_value.digits_.resize(x.digits_.size() + count);
        std::copy(x.digits_.begin(), x.digits_.end(), output_value.digits_.begin() + count);
        return output_value;
    }

    static BigInt shiftRight(const BigInt &x, Int count) {
        if (Int(x.digits_.size()) <= count)
            return {};
        BigInt output_value;
        output_value.sign_ = x.sign_;
        output_value.digits_.assign(x.digits_.begin() + count, x.digits_.end());
        return output_value;
    }

    static BigInt basePower(Int n) {
        BigInt output_value;
        output_value.sign_ = 1;
        output_value.digits_.resize(n + 1);
        output_value.digits_[n] = 1;
        return output_value;
    }

    static BigInt reciprocalAbs(const BigInt &);

    static BigInt quotientFromReciprocal(const BigInt &x, const BigInt &inverse, Int divisorLen,
                                         Int precision) {
        Int cut = std::max(Int(0), divisorLen - 1);
        return shiftRight(shiftRight(x, cut) * inverse, divisorLen + precision - cut);
    }

    static BigIntDivision correctDivisionProduct(const BigInt &, const BigInt &, BigInt,
                                                 const BigInt &);
    static BigIntDivision divmodBlocks(const BigInt &, const BigInt &);
    static BigIntDivision divmodNewton(const BigInt &, const BigInt &);

public:
    // 以下 n は対象の10進桁数、二項演算では両入力の最大桁数とする。
    // 零を構築する。計算量 O(1)。
    BigInt() = default;

    // 固定幅整数を変換する。計算量 O(1)。
    template <std::integral T> BigInt(T x) {
        UInt magnitude;
        if constexpr (std::is_signed_v<T>) {
            if (x < 0) {
                sign_ = -1;
                magnitude = UInt(-(x + 1)) + 1;
            } else
                magnitude = UInt(x);
        } else
            magnitude = UInt(x);
        if (magnitude && !sign_)
            sign_ = 1;
        while (magnitude) {
            digits_.push_back(magnitude % Base);
            magnitude /= Base;
        }
    }

    // 符号付き10進文字列を変換し、空文字列や不正な文字にはstd::invalid_argumentを送出する。O(文字数)。
    explicit BigInt(std::string_view s) {
        if (s.empty())
            throw std::invalid_argument("多倍長整数の文字列が空です");
        Int first = 0;
        sign_ = 1;
        if (s[0] == '+' || s[0] == '-') {
            if (s[0] == '-')
                sign_ = -1;
            first = 1;
        }
        if (first == Int(s.size()))
            throw std::invalid_argument("多倍長整数の数字がありません");
        digits_.resize((s.size() - first + 8) / 9);
        Int last = s.size(), index = 0;
        while (last > first) {
            Int start = std::max(first, last - 9);
            std::uint32_t digit = 0;
            for (Int i = start; i < last; ++i) {
                if (s[i] < '0' || s[i] > '9')
                    throw std::invalid_argument("多倍長整数に不正な文字が含まれています");
                digit = digit * 10 + (s[i] - '0');
            }
            digits_[index++] = digit;
            last = start;
        }
        normalize();
    }

    // 零か判定する。計算量 O(1)。
    bool isZero() const {
        return !sign_;
    }

    // 符号を -1, 0, 1 で返す。計算量 O(1)。
    int sgn() const {
        return sign_;
    }

    // 絶対値を返す。計算量 O(n)。
    BigInt abs() const {
        auto output_value = *this;
        if (output_value.sign_ < 0)
            output_value.sign_ = 1;
        return output_value;
    }

    // 符号反転した値を返す。計算量 O(n)。
    BigInt operator-() const {
        auto output_value = *this;
        output_value.sign_ = -output_value.sign_;
        return output_value;
    }

    // 値を複製する。計算量 O(n)。
    BigInt operator+() const {
        return *this;
    }

    // 十進文字列を返す。計算量 O(n)。
    std::string str() const {
        if (!sign_)
            return "0";
        if (digits_.size() <= 2) {
            UInt value = digits_[0];
            if (digits_.size() == 2)
                value += UInt(digits_[1]) * Base;
            return sign_ < 0 ? std::to_string(-Int(value)) : std::to_string(value);
        }
        constexpr std::string_view pairs =
            "00010203040506070809101112131415161718192021222324252627282930313233343536373839404142434445464748495051525354555657585960616263646566676869707172737475767778798081828384858687888990919293949596979899";
        auto top = std::to_string(digits_.back());
        Int signLen = sign_ < 0;
        std::string output_value((digits_.size() - 1) * 9 + top.size() + signLen, '\0');
        if (signLen)
            output_value[0] = '-';
        std::copy(top.begin(), top.end(), output_value.begin() + signLen);
        Int pos = output_value.size();
        for (std::size_t i = 0; i + 1 < digits_.size(); ++i) {
            auto digit = digits_[i];
            for (int j = 0; j < 4; ++j) {
                Int pair = (digit % 100) * 2;
                output_value[pos - 1] = pairs[pair + 1];
                output_value[pos - 2] = pairs[pair];
                pos -= 2;
                digit /= 100;
            }
            output_value[--pos] = char('0' + digit);
        }
        return output_value;
    }

    // 大小を比較し -1, 0, 1 を返す。計算量 O(n)。
    friend int cmp(const BigInt &a, const BigInt &b) {
        return a.sign_ != b.sign_ ? (a.sign_ < b.sign_ ? -1 : 1) : a.sign_ * cmpAbs(a, b);
    }

    // 等しいか判定する。計算量 O(n)。
    friend bool operator==(const BigInt &a, const BigInt &b) {
        return a.sign_ == b.sign_ && a.digits_ == b.digits_;
    }

    // 大小を比較する。計算量 O(n)。
    friend std::strong_ordering operator<=>(const BigInt &a, const BigInt &b) {
        return cmp(a, b) <=> 0;
    }

    // 和を返す。計算量 O(n)。
    friend BigInt operator+(const BigInt &a, const BigInt &b) {
        if (!a.sign_)
            return b;
        if (!b.sign_)
            return a;
        BigInt output_value;
        if (a.sign_ == b.sign_) {
            output_value = addAbs(a, b);
            output_value.sign_ = a.sign_;
        } else if (cmpAbs(a, b) >= 0) {
            output_value = subAbs(a, b);
            output_value.sign_ *= a.sign_;
        } else {
            output_value = subAbs(b, a);
            output_value.sign_ *= b.sign_;
        }
        return output_value;
    }

    // 差を返す。計算量 O(n)。
    friend BigInt operator-(const BigInt &a, const BigInt &b) {
        return a + (-b);
    }

    // 積を返す。NTT使用時 O(n log(n+1))、範囲外では筆算 O(n²)。
    friend BigInt operator*(const BigInt &a, const BigInt &b) {
        if (!a.sign_ || !b.sign_)
            return {};
        Int smaller = std::min(a.digits_.size(), b.digits_.size()),
            total = a.digits_.size() + b.digits_.size();
        if (smaller > NttThreshold) {
            if (smaller <= (LargeBaseMaxDigits * 2) / 3 && total <= (NttMaxLength * 2) / 3)
                return mulNtt<6>(a, b);
            if (total <= (NttMaxLength * 5) / 9)
                return mulNtt<5>(a, b);
        }
        return mulSchoolbook(a, b);
    }

    // 無限長の 2 の補数として論理積を返す。最大桁数を n として O(n^2)。
    friend BigInt operator&(const BigInt &a, const BigInt &b) {
        return bitwise<0>(a, b);
    }

    // 無限長の 2 の補数として論理和を返す。最大桁数を n として O(n^2)。
    friend BigInt operator|(const BigInt &a, const BigInt &b) {
        return bitwise<1>(a, b);
    }

    // 無限長の 2 の補数として排他的論理和を返す。最大桁数を n として O(n^2)。
    friend BigInt operator^(const BigInt &a, const BigInt &b) {
        return bitwise<2>(a, b);
    }

    // 無限長の 2 の補数としてビット反転した -x-1 を返す。桁数を n として O(n)。
    BigInt operator~() const {
        return -*this - BigInt(1);
    }

    // 2^shift倍を返す。負のシフト量はstd::invalid_argument。結果の桁数をnとしてO(n^2)。
    BigInt operator<<(Int shift) const {
        if (shift < 0)
            throw std::invalid_argument("シフト量は非負整数で指定してください");
        if (!shift || !sign_)
            return *this;
        auto digits = binaryDigits();
        Int words = shift / 32, bits = shift % 32;
        std::vector<std::uint32_t> output_value(digits.size() + words + 1);
        UInt carry = 0;
        for (std::size_t i = 0; i < digits.size(); ++i) {
            UInt value = (UInt(digits[i]) << bits) | carry;
            output_value[i + words] = std::uint32_t(value);
            carry = value >> 32;
        }
        output_value[digits.size() + words] = carry;
        auto result = fromBinary(output_value);
        result.sign_ *= sign_;
        return result;
    }

    // 算術右シフト（2^shiftによる床除算）を返す。負のシフト量はstd::invalid_argument。入力n桁でO(n^2)。
    BigInt operator>>(Int shift) const {
        if (shift < 0)
            throw std::invalid_argument("シフト量は非負整数で指定してください");
        if (!shift || !sign_)
            return *this;
        auto digits = binaryDigits();
        Int words = shift / 32, bits = shift % 32;
        if (words >= Int(digits.size()))
            return BigInt(sign_ < 0 ? -1 : 0);
        bool discarded = false;
        for (Int i = 0; i < words; ++i)
            discarded |= digits[i] != 0;
        discarded |= (UInt(digits[words]) & ((UInt(1) << bits) - 1)) != 0;
        std::vector<std::uint32_t> output_value(digits.size() - words);
        for (std::size_t i = 0; i < output_value.size(); ++i) {
            UInt value = UInt(digits[i + words]) >> bits;
            if (bits && i + words + 1 < digits.size())
                value |= UInt(digits[i + words + 1]) << (32 - bits);
            output_value[i] = std::uint32_t(value);
        }
        auto result = fromBinary(output_value);
        if (sign_ < 0) {
            if (discarded)
                result = result + BigInt(1);
            result.sign_ = -result.sign_;
        }
        return result;
    }

    // 0方向に丸めた商と余りを返す。最悪 O(n² log(n+1))。
    static BigIntDivision divmodTrunc(const BigInt &, const BigInt &);
    // 床除算の商と余りを返す。最悪 O(n² log(n+1))。
    static BigIntDivision divmod(const BigInt &, const BigInt &);
    // 0方向に丸めた商を返す。計算量は divmodTrunc と同じ。
    friend BigInt operator/(const BigInt &, const BigInt &);
    // 床除算の余りを返す。計算量は divmod と同じ。
    friend BigInt operator%(const BigInt &, const BigInt &);

    // 加算した結果を代入する。計算量は対応する二項演算と同じ。
    BigInt &operator+=(const BigInt &b) {
        return *this = *this + b;
    }

    // 減算した結果を代入する。計算量は対応する二項演算と同じ。
    BigInt &operator-=(const BigInt &b) {
        return *this = *this - b;
    }

    // 乗算した結果を代入する。計算量は対応する二項演算と同じ。
    BigInt &operator*=(const BigInt &b) {
        return *this = *this * b;
    }

    // 除算した結果を代入する。計算量は対応する二項演算と同じ。
    BigInt &operator/=(const BigInt &b) {
        return *this = *this / b;
    }

    // 床除算の余りした結果を代入する。計算量は対応する二項演算と同じ。
    BigInt &operator%=(const BigInt &b) {
        return *this = *this % b;
    }

    // ビット積した結果を代入する。計算量は対応する二項演算と同じ。
    BigInt &operator&=(const BigInt &b) {
        return *this = *this & b;
    }

    // ビット和した結果を代入する。計算量は対応する二項演算と同じ。
    BigInt &operator|=(const BigInt &b) {
        return *this = *this | b;
    }

    // 排他的ビット和した結果を代入する。計算量は対応する二項演算と同じ。
    BigInt &operator^=(const BigInt &b) {
        return *this = *this ^ b;
    }

    // 左シフトした結果を代入する。計算量は対応する二項演算と同じ。
    BigInt &operator<<=(Int b) {
        return *this = *this << b;
    }

    // 右シフトした結果を代入する。計算量は対応する二項演算と同じ。
    BigInt &operator>>=(Int b) {
        return *this = *this >> b;
    }

    // 非負整数乗を繰り返し二乗法で求め、負の指数にはstd::invalid_argumentを送出する。0^0は1。
    // 結果の桁数を N として最悪 O((N+1)² log(exponent+1))。
    BigInt pow(Int exponent) const {
        if (exponent < 0)
            throw std::invalid_argument("多倍長整数の指数は非負整数で指定してください");
        BigInt output_value = 1, base = *this;
        while (exponent) {
            if (exponent & 1)
                output_value *= base;
            exponent >>= 1;
            if (exponent)
                base *= base;
        }
        return output_value;
    }

    // Intに変換し、範囲外ならstd::overflow_errorを送出する。計算量 O(n)。
    Int toInt() const {
        UInt limit = UInt(std::numeric_limits<Int>::max()) + (sign_ < 0), magnitude = 0;
        for (Int i = Int(digits_.size()) - 1; i >= 0; --i) {
            UInt digit = digits_[i];
            if (magnitude > limit / Base || (magnitude == limit / Base && digit > limit % Base))
                throw std::overflow_error("多倍長整数が int の範囲外です");
            magnitude = magnitude * Base + digit;
        }
        if (sign_ < 0) {
            if (magnitude == UInt(std::numeric_limits<Int>::max()) + 1)
                return std::numeric_limits<Int>::min();
            return -Int(magnitude);
        }
        return Int(magnitude);
    }

    // 格納桁からハッシュ値を返す。計算量 O(n)。
    std::size_t hash() const {
        std::size_t output_value = std::hash<int>{}(sign_);
        for (auto d : digits_)
            output_value ^= std::hash<std::uint32_t>{}(d) + 0x9e3779b97f4a7c15ULL + (output_value << 6) + (output_value >> 2);
        return output_value;
    }

    // 十進文字列として出力する。計算量 O(n)。
    friend std::ostream &operator<<(std::ostream &os, const BigInt &x) {
        return os << x.str();
    }
};

inline BigIntDivision BigInt::divmodSchoolbook(const BigInt &x, const BigInt &y) {
    if (!y.sign_)
        throw std::domain_error("多倍長整数を 0 で割ることはできません");
    BigIntDivision output_value;
    if (cmpAbs(x, y) < 0) {
        output_value.remainder = x;
        return output_value;
    }
    if (x.digits_.size() <= 2) {
        UInt numerator = x.digits_[0], denominator = y.digits_[0];
        if (x.digits_.size() == 2)
            numerator += UInt(x.digits_[1]) * Base;
        if (y.digits_.size() == 2)
            denominator += UInt(y.digits_[1]) * Base;
        output_value.quotient = BigInt(numerator / denominator);
        output_value.remainder = BigInt(numerator % denominator);
        output_value.quotient.sign_ *= x.sign_ * y.sign_;
        output_value.remainder.sign_ *= x.sign_;
        return output_value;
    }
    if (y.digits_.size() == 1) {
        if (y.digits_[0] == 1) {
            output_value.quotient = x;
            output_value.quotient.sign_ *= y.sign_;
            return output_value;
        }
        auto [q, r] = divAbsSmall(x, y.digits_[0]);
        output_value.quotient = std::move(q);
        output_value.quotient.sign_ *= x.sign_ * y.sign_;
        output_value.remainder = BigInt(r);
        output_value.remainder.sign_ *= x.sign_;
        return output_value;
    }
    auto factor = std::uint32_t(Base / (UInt(y.digits_.back()) + 1));
    BigInt dividend = mulAbsSmall(x, factor), divisor = mulAbsSmall(y, factor);
    Int n = divisor.digits_.size(), qlen = Int(x.digits_.size()) - n + 1;
    dividend.digits_.resize(x.digits_.size() + 1);
    output_value.quotient.sign_ = x.sign_ * y.sign_;
    output_value.quotient.digits_.resize(qlen);
    for (Int j = qlen - 1; j >= 0; --j) {
        UInt top = UInt(dividend.digits_[j + n]) * Base + dividend.digits_[j + n - 1],
             estimate = std::min(Base - 1, top / divisor.digits_[n - 1]),
             remainder = top - estimate * divisor.digits_[n - 1];
        while (remainder < Base &&
               estimate * divisor.digits_[n - 2] > remainder * Base + dividend.digits_[j + n - 2]) {
            --estimate;
            remainder += divisor.digits_[n - 1];
        }
        Int borrow = 0;
        for (Int i = 0; i < n; ++i) {
            UInt product = estimate * divisor.digits_[i] + UInt(borrow);
            Int digit = Int(dividend.digits_[j + i]) - Int(product % Base);
            borrow = product / Base;
            if (digit < 0) {
                digit += Base;
                ++borrow;
            }
            dividend.digits_[j + i] = digit;
        }
        Int last = Int(dividend.digits_[j + n]) - borrow;
        if (last < 0) {
            --estimate;
            UInt carry = 0;
            for (Int i = 0; i < n; ++i) {
                UInt digit = UInt(dividend.digits_[j + i]) + divisor.digits_[i] + carry;
                dividend.digits_[j + i] = digit % Base;
                carry = digit / Base;
            }
            last += Int(carry);
        }
        dividend.digits_[j + n] = last;
        output_value.quotient.digits_[j] = estimate;
    }
    output_value.quotient.normalize();
    dividend.digits_.resize(n);
    dividend.normalize();
    output_value.remainder = divAbsSmall(dividend, factor).first;
    output_value.remainder.sign_ *= x.sign_;
    return output_value;
}

// 上位半分にガード2桁を加え、B^(2m)/xの下側・誤差2未満の逆数を求める。
inline BigInt BigInt::reciprocalAbs(const BigInt &x) {
    Int m = x.digits_.size();
    if (m <= 32)
        return divmodSchoolbook(basePower(2 * m), x).quotient;
    Int half = (m + 1) / 2 + 2, low = m - half;
    auto inverse = reciprocalAbs(shiftRight(x, low));
    auto product = x * (inverse * inverse);
    auto correction = shiftRight(product, 2 * half);
    for (Int i = 0; i < std::min(2 * half, Int(product.digits_.size())); ++i)
        if (product.digits_[i]) {
            correction += 1;
            break;
        }
    return shiftLeft(mulAbsSmall(inverse, 2), low) - correction;
}

inline BigIntDivision BigInt::correctDivisionProduct(const BigInt &x, const BigInt &y,
                                                     BigInt quotient, const BigInt &product) {
    BigIntDivision output_value{std::move(quotient), x - product};
    if (output_value.remainder.sign_ < 0) {
        output_value.quotient -= 1;
        output_value.remainder += y;
    }
    for (int i = 0; i < 3; ++i) {
        if (cmpAbs(output_value.remainder, y) < 0)
            break;
        output_value.quotient += 1;
        output_value.remainder -= y;
    }
    return output_value;
}

// 長い被除数を除数長のブロックに分け、逆数と固定側NTTを再利用する。
inline BigIntDivision BigInt::divmodBlocks(const BigInt &x, const BigInt &y) {
    Int n = x.digits_.size(), m = y.digits_.size();
    auto inverse = reciprocalAbs(y);
    FixedMultiplier inverseMultiplier, divisorMultiplier;
    if (n >= 4 * m) {
        inverseMultiplier = FixedMultiplier(inverse, m + 1);
        divisorMultiplier = FixedMultiplier(y, m);
    }
    BigIntDivision output_value;
    output_value.quotient.sign_ = 1;
    output_value.quotient.digits_.resize(n - m + 1);
    for (Int block = (n - 1) / m; block >= 0; --block) {
        Int offset = block * m, width = std::min(m, n - offset);
        BigInt dividend;
        dividend.sign_ = 1;
        dividend.digits_.resize(width + output_value.remainder.digits_.size());
        std::copy_n(x.digits_.begin() + offset, width, dividend.digits_.begin());
        std::copy(output_value.remainder.digits_.begin(), output_value.remainder.digits_.end(),
                  dividend.digits_.begin() + width);
        dividend.normalize();
        if (cmpAbs(dividend, y) < 0) {
            output_value.remainder = std::move(dividend);
            continue;
        }
        auto quotient =
            shiftRight(inverseMultiplier.multiply(shiftRight(dividend, m - 1), inverse), m + 1);
        auto division =
            correctDivisionProduct(dividend, y, quotient, divisorMultiplier.multiply(quotient, y));
        std::copy(division.quotient.digits_.begin(), division.quotient.digits_.end(),
                  output_value.quotient.digits_.begin() + offset);
        output_value.remainder = std::move(division.remainder);
    }
    output_value.quotient.normalize();
    return output_value;
}

inline BigIntDivision BigInt::divmodNewton(const BigInt &x, const BigInt &y) {
    auto dividend = x.abs(), divisor = y.abs();
    Int n = dividend.digits_.size(), m = divisor.digits_.size();
    BigIntDivision output_value;
    if (n > 2 * m)
        output_value = divmodBlocks(dividend, divisor);
    else {
        Int cut = std::max(Int(0), m - (n - m + 2));
        auto high = shiftRight(divisor, cut);
        Int precision = std::max(m - cut, n - m);
        auto inverse = reciprocalAbs(shiftLeft(high, precision - (m - cut)));
        auto quotient =
            quotientFromReciprocal(shiftRight(dividend, cut), inverse, m - cut, precision);
        output_value = correctDivisionProduct(dividend, divisor, quotient, quotient * divisor);
    }
    output_value.quotient.sign_ *= x.sign_ * y.sign_;
    output_value.remainder.sign_ *= x.sign_;
    return output_value;
}

// 0方向に丸めた商と被除数と同符号の余りを返し、0除算にはstd::domain_errorを送出する。
inline BigIntDivision BigInt::divmodTrunc(const BigInt &x, const BigInt &y) {
    Int n = x.digits_.size(), m = y.digits_.size();
    if (std::min(m, n - m) > NewtonThreshold) {
        Int size = n > 2 * m ? 2 * m + 8 : std::max(n + 1, 2 * std::min(m, n - m + 2) + 8);
        if (size <= (NttMaxLength * 5) / 9)
            return divmodNewton(x, y);
    }
    return divmodSchoolbook(x, y);
}

// Pythonと同じ床除算の商と除数と同符号の余りを返す。0除算はstd::domain_error。
inline BigIntDivision BigInt::divmod(const BigInt &x, const BigInt &y) {
    auto output_value = divmodTrunc(x, y);
    if (!output_value.remainder.isZero() && x.sign_ != y.sign_) {
        output_value.quotient -= 1;
        output_value.remainder += y;
    }
    return output_value;
}

inline BigInt operator/(const BigInt &x, const BigInt &y) {
    return BigInt::divmodTrunc(x, y).quotient;
}

// 元の%と同じ床除算の余り。0方向の余りはmod()で取得する。
inline BigInt operator%(const BigInt &x, const BigInt &y) {
    return BigInt::divmod(x, y).remainder;
}

// 床除算の商と余りを返す。計算量は BigInt::divmod と同じ。
inline BigIntDivision divmod(const BigInt &x, const BigInt &y) {
    return BigInt::divmod(x, y);
}

// 0方向に丸めた商と余りを返す。計算量は BigInt::divmodTrunc と同じ。
inline BigIntDivision divmodTrunc(const BigInt &x, const BigInt &y) {
    return BigInt::divmodTrunc(x, y);
}

// 0方向に丸めた商を返す。計算量は x / y と同じ。
inline BigInt div(const BigInt &x, const BigInt &y) {
    return x / y;
}

// 0方向に丸めた余りを返す。計算量は BigInt::divmodTrunc と同じ。
inline BigInt mod(const BigInt &x, const BigInt &y) {
    return BigInt::divmodTrunc(x, y).remainder;
}

// 床除算の商を返す。計算量は BigInt::divmod と同じ。
inline BigInt floor_div(const BigInt &x, const BigInt &y) {
    return BigInt::divmod(x, y).quotient;
}

// 0方向に丸めた商を x に代入する。計算量は x / y と同じ。
inline BigInt &div_assign(BigInt &x, const BigInt &y) {
    return x = x / y;
}

// 0方向に丸めた余りを x に代入する。計算量は mod(x, y) と同じ。
inline BigInt &mod_assign(BigInt &x, const BigInt &y) {
    return x = mod(x, y);
}

// 十進文字列を変換する。計算量 O(文字数)。
inline BigInt parseBigInt(std::string_view s) {
    return BigInt(s);
}

// 指定された値から構築する。計算量は対応する BigInt のコンストラクタと同じ。
template <class T> BigInt initBigInt(const T &x) {
    return BigInt(x);
}

// 固定幅整数を変換する。計算量 O(1)。
template <std::integral T> BigInt toBigInt(T x) {
    return BigInt(x);
}

// 絶対値を返す。桁数 n に対し O(n)。
inline BigInt abs(const BigInt &x) {
    return x.abs();
}

// 符号を返す。計算量 O(1)。
inline int sgn(const BigInt &x) {
    return x.sgn();
}

// 零か判定する。計算量 O(1)。
inline bool isZero(const BigInt &x) {
    return x.isZero();
}

// Int に変換する。桁数 n に対し O(n)。
inline Int toInt(const BigInt &x) {
    return x.toInt();
}

// ハッシュ値を返す。桁数 n に対し O(n)。
inline std::size_t hashValue(const BigInt &x) {
    return x.hash();
}

// 非負整数乗を返す。計算量は x.pow(n) と同じ。
inline BigInt pow(const BigInt &x, Int n) {
    return x.pow(n);
}

// 非負の最大公約数を返す。gcd(0, 0) は 0。最大桁数 n に対し最悪 O(n³ log(n+1))。
inline BigInt gcd(const BigInt &x, const BigInt &y) {
    BigInt a = x.abs(), b = y.abs();
    while (!b.isZero()) {
        auto r = mod(a, b);
        a = std::move(b);
        b = std::move(r);
    }
    return a;
}

// 非負の最小公倍数を返す。いずれかが 0 なら 0。最大桁数 n に対し最悪 O(n³ log(n+1))。
inline BigInt lcm(const BigInt &x, const BigInt &y) {
    if (x.isZero() || y.isZero())
        return {};
    return ((x / gcd(x, y)) * y).abs();
}

// 10進数の_biリテラルを多倍長整数に変換する。桁区切りの_も使用できる。計算量 O(文字数)。
inline BigInt operator""_bi(const char *s) {
    std::string output_value;
    for (; *s; ++s)
        if (*s != '\'' && *s != '_')
            output_value += *s;
    return BigInt(output_value);
}

// 文字列リテラルを多倍長整数に変換する。計算量 O(n)。
inline BigInt operator""_bi(const char *s, std::size_t n) {
    std::string output_value;
    for (std::size_t i = 0; i < n; ++i)
        if (s[i] != '\'' && s[i] != '_')
            output_value += s[i];
    return BigInt(output_value);
}
}

namespace std {
template <> struct hash<cplib::BigInt> {
    // 多倍長整数のハッシュ値を返す。桁数 n に対し O(n)。
    size_t operator()(const cplib::BigInt &x) const {
        return x.hash();
    }
};
}
