#pragma once
#include <cplib/utils/backwards_index.hpp>
#include <stdexcept>
#include <compare>

namespace cplib {
// Size == -1 は動的配列、その他は固定長配列。両版のシフト手順を保持。
template <Int Size = -1> class ScalarBitSet {
    static_assert(Size >= -1);
    // clang-format off
    using Storage = std::conditional_t<Size == -1, std::vector<UInt>,
                                       std::array<UInt, Size < 0 ? 0 : (Size + 63) / 64>>;
    // clang-format on
    Storage bits_{};
    Int size_ = Size < 0 ? 0 : Size;

    void same(const ScalarBitSet &other) const {
        if constexpr (Size < 0)
            if (size_ != other.size_)
                throw std::invalid_argument("BitSet sizes must match");
    }

    void check(Int i) const {
        if (i < 0 || i >= size_)
            throw std::out_of_range("BitSet index out of bounds");
    }

    void trim() {
        Int rem = size_ & 63;
        if (rem)
            bits_.back() &= (UInt(1) << rem) - 1;
    }

    ScalarBitSet empty() const {
        if constexpr (Size < 0)
            return ScalarBitSet(size_);
        else
            return {};
    }

public:
    ScalarBitSet() = default;

    explicit ScalarBitSet(Int n)
        requires(Size < 0)
        : size_(n) {
        if (n < 0)
            throw std::invalid_argument("BitSet size must be non-negative");
        bits_.resize((n + 63) / 64);
    }

    ScalarBitSet(const std::vector<bool> &v, Int n)
        requires(Size < 0)
        : ScalarBitSet(n) {
        if (Int(v.size()) > n)
            throw std::invalid_argument("initial value is longer than BitSet size");
        for (Int i = 0; i < Int(v.size()); ++i)
            if (v[i])
                bits_[i >> 6] |= UInt(1) << (i & 63);
    }

    explicit ScalarBitSet(const std::vector<bool> &v) {
        if constexpr (Size < 0) {
            size_ = v.size();
            bits_.resize((size_ + 63) / 64);
        }
        if (Int(v.size()) > size_)
            throw std::invalid_argument("initial value is longer than BitSet size");
        for (Int i = 0; i < Int(v.size()); ++i)
            if (v[i])
                bits_[i >> 6] |= UInt(1) << (i & 63);
    }

    // ビット数をO(1)で返す。
    Int len() const {
        return size_;
    }

    std::span<const UInt> words() const {
        return bits_;
    }

    bool operator[](Int i) const {
        check(i);
        return (bits_[i >> 6] >> (i & 63)) & 1;
    }

    void set(Int i, bool value) {
        check(i);
        UInt mask = UInt(1) << (i & 63);
        if (value)
            bits_[i >> 6] |= mask;
        else
            bits_[i >> 6] &= ~mask;
    }

    struct Reference {
        ScalarBitSet *owner;
        Int index;

        operator bool() const {
            return std::as_const(*owner)[index];
        }

        Reference &operator=(bool v) {
            owner->set(index, v);
            return *this;
        }

        Reference &operator=(Int v) {
            if (v == 0 || v == 1)
                owner->set(index, v == 1);
            return *this;
        }

        Reference &operator=(int v) {
            return *this = Int(v);
        }

        Reference &operator=(const Reference &v) {
            return *this = bool(v);
        }
    };

    Reference operator[](Int i) {
        return {this, i};
    }

    CPLIB_BACKWARDS_INDEX_OVERLOADS
    ScalarBitSet operator&(const ScalarBitSet &y) const {
        same(y);
        auto out = empty();
        for (std::size_t i = 0; i < bits_.size(); ++i)
            out.bits_[i] = bits_[i] & y.bits_[i];
        return out;
    }

    ScalarBitSet operator|(const ScalarBitSet &y) const {
        same(y);
        auto out = empty();
        for (std::size_t i = 0; i < bits_.size(); ++i)
            out.bits_[i] = bits_[i] | y.bits_[i];
        return out;
    }

    ScalarBitSet operator^(const ScalarBitSet &y) const {
        same(y);
        auto out = empty();
        for (std::size_t i = 0; i < bits_.size(); ++i)
            out.bits_[i] = bits_[i] ^ y.bits_[i];
        return out;
    }

    ScalarBitSet &operator&=(const ScalarBitSet &y) {
        same(y);
        for (std::size_t i = 0; i < bits_.size(); ++i)
            bits_[i] &= y.bits_[i];
        return *this;
    }

    ScalarBitSet &operator|=(const ScalarBitSet &y) {
        same(y);
        for (std::size_t i = 0; i < bits_.size(); ++i)
            bits_[i] |= y.bits_[i];
        return *this;
    }

    ScalarBitSet &operator^=(const ScalarBitSet &y) {
        same(y);
        for (std::size_t i = 0; i < bits_.size(); ++i)
            bits_[i] ^= y.bits_[i];
        return *this;
    }

    ScalarBitSet operator>>(Int x) const {
        if (x < 0)
            throw std::invalid_argument("shift count must be non-negative");
        auto out = empty();
        if (x >= size_)
            return out;
        Int w = x >> 6, b = x & 63;
        if constexpr (Size < 0) {
            for (Int dst = 0; dst < Int(bits_.size()) - w; ++dst) {
                Int src = dst + w;
                out.bits_[dst] = bits_[src] >> b;
                if (b && src + 1 < Int(bits_.size()))
                    out.bits_[dst] |= bits_[src + 1] << (64 - b);
            }
        } else {
            for (Int i = 0; i < Int(bits_.size()); ++i)
                if (i + w < Int(bits_.size()))
                    out.bits_[i] = bits_[i + w];
            UInt tmp = 0;
            if (b)
                for (Int i = bits_.size(); i-- > 0;) {
                    UInt mask = out.bits_[i] & ((UInt(1) << b) - 1);
                    out.bits_[i] >>= b;
                    out.bits_[i] |= tmp << (64 - b);
                    tmp = mask;
                }
        }
        return out;
    }

    // 添字が大きい方向へxビットずらし、範囲外を切り捨てる。
    ScalarBitSet operator<<(Int x) const {
        if (x < 0)
            throw std::invalid_argument("shift count must be non-negative");
        auto out = empty();
        if (x >= size_)
            return out;
        Int w = x >> 6, b = x & 63;
        if constexpr (Size < 0) {
            for (Int dst = bits_.size(); dst-- > w;) {
                Int src = dst - w;
                out.bits_[dst] = bits_[src] << b;
                if (b && src > 0)
                    out.bits_[dst] |= bits_[src - 1] >> (64 - b);
            }
        } else {
            for (Int i = 0; i < Int(bits_.size()); ++i)
                if (i - w >= 0)
                    out.bits_[i] = bits_[i - w];
            UInt tmp = 0;
            if (b)
                for (Int i = 0; i < Int(bits_.size()); ++i) {
                    UInt mask = out.bits_[i];
                    out.bits_[i] <<= b;
                    out.bits_[i] |= tmp >> (64 - b);
                    tmp = mask;
                }
        }
        out.trim();
        return out;
    }

    ScalarBitSet operator~() const {
        auto out = empty();
        if constexpr (Size < 0) {
            for (std::size_t i = 0; i < bits_.size(); ++i)
                out.bits_[i] = ~bits_[i];
            out.trim();
        } else if constexpr (Size > 0) {
            for (std::size_t i = 0; i + 1 < bits_.size(); ++i)
                out.bits_[i] = ~bits_[i];
            if constexpr (Size % 64 == 0)
                out.bits_.back() = ~bits_.back();
            else
                out.bits_.back() = bits_.back() ^ ((UInt(1) << (Size % 64)) - 1);
        }
        return out;
    }

    Int popcount() const {
        Int result = 0;
        for (UInt w : bits_)
            result += std::popcount(w);
        return result;
    }

    Int andpopcount(const ScalarBitSet &y) const {
        same(y);
        Int out = 0;
        for (std::size_t i = 0; i < bits_.size(); ++i)
            out += std::popcount(bits_[i] & y.bits_[i]);
        return out;
    }

    Int orpopcount(const ScalarBitSet &y) const {
        same(y);
        Int out = 0;
        for (std::size_t i = 0; i < bits_.size(); ++i)
            out += std::popcount(bits_[i] | y.bits_[i]);
        return out;
    }

    Int xorpopcount(const ScalarBitSet &y) const {
        same(y);
        Int out = 0;
        for (std::size_t i = 0; i < bits_.size(); ++i)
            out += std::popcount(bits_[i] ^ y.bits_[i]);
        return out;
    }

    Int lowestBit() const {
        for (std::size_t i = 0; i < bits_.size(); ++i)
            if (bits_[i]) {
                Int index = i * 64 + std::countr_zero(bits_[i]);
                if (index < size_)
                    return index;
            }
        return -1;
    }

    // 同じ長さのビット列を添字0からfalse < trueで比較し、-1・0・1を返します。
    // 時間O(1 + N / 64)、追加メモリO(1)。最初の相違で終了します。
    int cmp(const ScalarBitSet &y) const {
        same(y);
        for (std::size_t i = 0; i < bits_.size(); ++i) {
            UInt diff = bits_[i] ^ y.bits_[i];
            if (diff)
                return ((bits_[i] >> std::countr_zero(diff)) & 1) ? 1 : -1;
        }
        return 0;
    }

    // 添字0からfalse < trueの辞書順で小さいかを返します。時間O(1 + N / 64)、追加メモリO(1)。
    bool lexLess(const ScalarBitSet &y) const {
        return cmp(y) < 0;
    }

    auto operator<=>(const ScalarBitSet &y) const {
        return cmp(y) <=> 0;
    }

    bool operator==(const ScalarBitSet &y) const {
        return size_ == y.size_ && bits_ == y.bits_;
    }

    // 有効な全ビットが1かを返します。0を見つけたら終了し、長さ0ではtrueを返します。
    // 最悪時間O(1 + N / 64)、追加メモリO(1)。
    bool all() const {
        Int full = size_ >> 6;
        for (Int i = 0; i < full; ++i)
            if (bits_[i] != ~UInt(0))
                return false;
        Int rem = size_ & 63;
        return !rem || (bits_[full] & ((UInt(1) << rem) - 1)) == ((UInt(1) << rem) - 1);
    }

    // 有効なビットに1があるかを返します。1を見つけたら終了し、長さ0ではfalseを返します。
    // 最悪時間O(1 + N / 64)、追加メモリO(1)。
    bool any() const {
        Int full = size_ >> 6;
        for (Int i = 0; i < full; ++i)
            if (bits_[i])
                return true;
        Int rem = size_ & 63;
        return rem && (bits_[full] & ((UInt(1) << rem) - 1));
    }

    std::string to_string() const {
        std::string out(size_, '0');
        for (Int i = 0; i < size_; ++i)
            if ((*this)[i])
                out[size_ - 1 - i] = '1';
        return out;
    }

    struct Iterator {
        const ScalarBitSet *owner;
        std::size_t word;
        UInt remaining;

        void skip() {
            while (!remaining && word < owner->bits_.size())
                if (++word < owner->bits_.size())
                    remaining = owner->bits_[word];
        }

        Int operator*() const {
            return word * 64 + std::countr_zero(remaining);
        }

        Iterator &operator++() {
            remaining &= remaining - 1;
            skip();
            return *this;
        }

        bool operator==(const Iterator &other) const {
            return owner == other.owner && word == other.word && remaining == other.remaining;
        }
    };

    Iterator begin() const {
        Iterator it{this, 0, bits_.empty() ? 0 : bits_[0]};
        it.skip();
        return it;
    }

    Iterator end() const {
        return {this, bits_.size(), 0};
    }

    const ScalarBitSet &items() const {
        return *this;
    }
};

template <Int N> Int len(const ScalarBitSet<N> &s) {
    return s.len();
}

template <Int N> Int popcount(const ScalarBitSet<N> &s) {
    return s.popcount();
}

template <Int N> Int lowestBit(const ScalarBitSet<N> &s) {
    return s.lowestBit();
}

template <Int N> bool all(const ScalarBitSet<N> &s) {
    return s.all();
}

template <Int N> bool any(const ScalarBitSet<N> &s) {
    return s.any();
}

template <Int N> int cmp(const ScalarBitSet<N> &s, const ScalarBitSet<N> &t) {
    return s.cmp(t);
}

template <Int N> bool lexLess(const ScalarBitSet<N> &s, const ScalarBitSet<N> &t) {
    return s.lexLess(t);
}

template <Int N> Int andpopcount(const ScalarBitSet<N> &s, const ScalarBitSet<N> &t) {
    return s.andpopcount(t);
}

template <Int N> Int orpopcount(const ScalarBitSet<N> &s, const ScalarBitSet<N> &t) {
    return s.orpopcount(t);
}

template <Int N> Int xorpopcount(const ScalarBitSet<N> &s, const ScalarBitSet<N> &t) {
    return s.xorpopcount(t);
}

template <Int N> std::string to_string(const ScalarBitSet<N> &s) {
    return s.to_string();
}
}
