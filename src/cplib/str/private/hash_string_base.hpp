#pragma once
#include <cplib/utils/backwards_index.hpp>
#include <memory>
#include <random>
#include <string_view>

namespace cplib::detail {
inline constexpr UInt hash_string_mod = (UInt(1) << 61) - 1;

inline UInt hash_string_modulo(UInt x) {
    x = (x >> 61) + (x & hash_string_mod);
    if (x >= hash_string_mod)
        x -= hash_string_mod;
    return x;
}

// 31/30 bit 分割による元の乗算。128bit乗算への置き換えはしない。
inline UInt hash_string_mul(UInt a, UInt b) {
    UInt au = a >> 31, al = a & ((UInt(1) << 31) - 1), bu = b >> 31, bl = b & ((UInt(1) << 31) - 1),
         mid = al * bu + au * bl;
    return au * bu * 2 + (mid >> 30) + ((mid & ((UInt(1) << 30) - 1)) << 31) + al * bl;
}

inline UInt hash_string_pow(UInt a, Int n) {
    UInt result = 1;
    while (n > 0) {
        if (n & 1)
            result = hash_string_modulo(hash_string_mul(result, a));
        a = hash_string_modulo(hash_string_mul(a, a));
        n >>= 1;
    }
    return result;
}

template <bool Reverse> struct HashStringConfig {
    UInt base, inverse;
    std::vector<UInt> pows, invpows;

    HashStringConfig() : pows(500001), invpows(500001) {
        std::mt19937_64 rng(std::random_device{}());
        base = std::uniform_int_distribution<UInt>(129, UInt(1) << 30)(rng);
        inverse = hash_string_pow(base, hash_string_mod - 2);
        pows[0] = invpows[0] = 1;
        for (std::size_t i = 1; i < pows.size(); ++i) {
            pows[i] = hash_string_modulo(hash_string_mul(pows[i - 1], base));
            invpows[i] = hash_string_modulo(hash_string_mul(invpows[i - 1], inverse));
        }
    }

    UInt power(Int n) const {
        assert(n >= 0);
        return n < Int(pows.size()) ? pows[n] : hash_string_pow(base, n);
    }

    UInt invpower(Int n) const {
        assert(n >= 0);
        return n < Int(invpows.size()) ? invpows[n] : hash_string_pow(inverse, n);
    }
};

template <bool R> inline HashStringConfig<R> hash_string_config;

template <bool Reverse> struct ReverseHashField {};

template <> struct ReverseHashField<true> {
    UInt rhash = 0;
};
}

namespace cplib {
// 文字列の長さとハッシュを保持する。ハッシュの衝突により異なる文字列を等しいと判定する可能性がある。
template <bool Reverse> struct BasicHashString : detail::ReverseHashField<Reverse> {
    UInt hash = 0, bpow = 1;
    Int length_ = 0;

    // 1要素からハッシュを作る。計算量 O(1)。
    static BasicHashString character(Int s) {
        if constexpr (Reverse)
            assert(0 <= s && UInt(s) < detail::hash_string_mod);
        BasicHashString out;
        out.hash = UInt(s) % detail::hash_string_mod;
        if constexpr (Reverse)
            out.rhash = out.hash;
        out.bpow = detail::hash_string_config<Reverse>.base;
        out.length_ = 1;
        return out;
    }

    // 長さNの列からハッシュを作る。計算量 O(N+1)。
    template <class Range> static BasicHashString sequence(const Range &s) {
        BasicHashString out;
        UInt tmp = 1;
        for (Int i = Int(s.size()); i-- > 0;) {
            auto value = [](auto x) {
                if constexpr (std::is_same_v<decltype(x), char>)
                    return Int(static_cast<unsigned char>(x));
                else
                    return Int(x);
            };
            out.hash = detail::hash_string_modulo(
                out.hash + detail::hash_string_mul(character(value(s[i])).hash, tmp));
            if constexpr (Reverse)
                out.rhash = detail::hash_string_modulo(
                    out.rhash +
                    detail::hash_string_mul(character(value(s[s.size() - 1 - i])).hash, tmp));
            tmp = detail::hash_string_modulo(
                detail::hash_string_mul(tmp, detail::hash_string_config<Reverse>.base));
        }
        out.bpow = detail::hash_string_config<Reverse>.power(s.size());
        out.length_ = s.size();
        return out;
    }

    // 列の長さを返す。計算量 O(1)。
    Int size() const {
        return length_;
    }

    // 2つのハッシュを結合する。計算量 O(1)。
    friend BasicHashString operator&(const BasicHashString &l, const BasicHashString &r) {
        BasicHashString out;
        out.hash = detail::hash_string_modulo(
            detail::hash_string_modulo(detail::hash_string_mul(l.hash, r.bpow)) + r.hash);
        if constexpr (Reverse)
            out.rhash = detail::hash_string_modulo(
                detail::hash_string_modulo(detail::hash_string_mul(r.rhash, l.bpow)) + l.rhash);
        out.bpow = detail::hash_string_modulo(detail::hash_string_mul(l.bpow, r.bpow));
        out.length_ = l.length_ + r.length_;
        return out;
    }

    // 長さとハッシュが等しいか判定する。計算量 O(1)。
    bool operator==(const BasicHashString &r) const {
        return length_ == r.length_ && hash == r.hash;
    }

    // 同じ列をx回繰り返したハッシュを返す。計算量 O(log(x+1))。
    BasicHashString operator*(Int x) const {
        BasicHashString out, tmp = *this;
        Int total = length_ * x;
        while (x > 0) {
            if (x & 1)
                out = out & tmp;
            if (x > 1)
                tmp = tmp & tmp;
            x >>= 1;
        }
        out.length_ = total;
        return out;
    }

    // 指定した接頭辞のハッシュを取り除く。残り長が表の範囲内なら O(1)、範囲外は O(log(残り長+1))。
    BasicHashString removePrefix(const BasicHashString &prefix) const
        requires(!Reverse)
    {
        Int n = length_ - prefix.length_;
        BasicHashString out;
        out.hash = detail::hash_string_modulo(
            hash + detail::hash_string_mod -
            detail::hash_string_modulo(detail::hash_string_mul(
                prefix.hash, detail::hash_string_config<Reverse>.power(n))));
        out.bpow = detail::hash_string_config<Reverse>.power(n);
        out.length_ = n;
        return out;
    }

    // 正逆のハッシュを比べて回文か判定する。計算量 O(1)。
    bool isPalindrome() const
        requires Reverse
    {
        return hash == this->rhash;
    }

    // 反転した列のハッシュを返す。計算量 O(1)。
    BasicHashString reversed() const
        requires Reverse
    {
        auto out = *this;
        std::swap(out.hash, out.rhash);
        return out;
    }
};

template <bool Reverse> struct BasicRollingHash {
    using Hash = BasicHashString<Reverse>;

    struct Base {
        std::string S;
        std::vector<UInt> prefixs, rprefixs;

        // 基底文字列の長さを返す。計算量 O(1)。
        Int size() const {
            return S.size();
        }

        // 基底文字列の長さを返す。計算量 O(1)。
    };

    std::shared_ptr<Base> R;
    Int l = 0, r = 0;
    // 空のハンドルを作る。計算量 O(1)。
    BasicRollingHash() = default;

    // 長さNの文字列をコピーして累積ハッシュを作る。計算量 O(N+1)。
    explicit BasicRollingHash(std::string_view s) : R(std::make_shared<Base>()), r(s.size()) {
        R->S = s;
        R->prefixs.resize(s.size() + 1);
        if constexpr (Reverse)
            R->rprefixs.resize(s.size() + 1);
        UInt tmp = 1;
        for (std::size_t i = 1; i <= s.size(); ++i) {
            UInt ch = static_cast<unsigned char>(s[i - 1]);
            R->prefixs[i] = detail::hash_string_modulo(
                detail::hash_string_mul(R->prefixs[i - 1],
                                        detail::hash_string_config<Reverse>.base) +
                ch);
            if constexpr (Reverse) {
                R->rprefixs[i] = detail::hash_string_modulo(detail::hash_string_mul(ch, tmp) +
                                                            R->rprefixs[i - 1]);
                tmp = detail::hash_string_modulo(
                    detail::hash_string_mul(tmp, detail::hash_string_config<Reverse>.base));
            }
        }
    }

    // ビューの長さを返す。計算量 O(1)。
    Int size() const {
        return r - l;
    }

    // i番目の文字を返す。計算量 O(1)。
    char operator[](Int i) const {
        assert(0 <= i && i < size());
        return R->S[l + i];
    }

    // 末尾基準の添字で文字を返す。計算量 O(1)。
    char operator[](BackwardsIndex i) const {
        return (*this)[size() - i.value];
    }

    // 半開区間[a,b)の部分文字列を返す。空区間は返り値のl = r = 0へ正規化する。計算量 O(1)。
    BasicRollingHash substr(Int a, Int b) const {
        assert(0 <= a && a <= b && b <= size());
        auto out = *this;
        if (a == b)
            out.l = out.r = 0;
        else {
            out.l = l + a;
            out.r = l + b;
        }
        return out;
    }

    // 両端を含む区間のビューを返す。計算量 O(1)。
    template <class A, class B> auto operator[](ClosedSlice<A, B> s) const {
        return substr(resolve_index(size(), s.a), resolve_index(size(), s.b) + 1);
    }

    // 半開区間のハッシュを返す。長さが累乗表の範囲内なら O(1)、範囲外は O(log(長さ+1))。
    UInt gethash(Int a, Int b) const {
        return detail::hash_string_modulo(
            R->prefixs[l + b] + detail::hash_string_mod -
            detail::hash_string_modulo(detail::hash_string_mul(
                R->prefixs[l + a], detail::hash_string_config<Reverse>.power(b - a))));
    }

    // ビューのハッシュを作る。長さと開始位置が累乗表の範囲内なら O(1)、範囲外は O(log(基底長+1))。
    operator Hash() const {
        Hash out;
        out.hash = gethash(0, size());
        out.bpow = detail::hash_string_config<Reverse>.power(size());
        out.length_ = size();
        if constexpr (Reverse)
            out.rhash = detail::hash_string_modulo(detail::hash_string_mul(
                detail::hash_string_modulo(R->rprefixs[r] + detail::hash_string_mod -
                                           R->rprefixs[l]),
                detail::hash_string_config<Reverse>.invpower(l)));
        return out;
    }

    // 2つのビューの長さとハッシュが等しいか判定する。計算量はgethashの2回分。
    bool operator==(const BasicRollingHash &t) const {
        return size() == t.size() && gethash(0, size()) == t.gethash(0, t.size());
    }

    // ビューの文字列をコピーする。長さNとして O(N+1)。
    std::string to_string() const {
        return R->S.substr(l, size());
    }

    // 2つのビューのハッシュを結合する。計算量は各ビューのHash変換と O(1)。
    friend Hash operator&(const BasicRollingHash &a, const BasicRollingHash &b) {
        return Hash(a) & Hash(b);
    }

    // 正逆のハッシュを比べて回文か判定する。計算量はHash変換と同じ。
    bool isPalindrome() const
        requires Reverse
    {
        return Hash(*this).isPalindrome();
    }
};

// ビューのハッシュを作る。計算量はHash変換と同じ。
template <bool R> auto toHashString(const BasicRollingHash<R> &s) {
    return BasicHashString<R>(s);
}

// ハッシュを二分探索して最長共通接頭辞長を返す。短い方の長さNとして表内は O(log(N+1))、表外は O(log²(N+1))。
template <bool R> Int LCP(const BasicRollingHash<R> &s, const BasicRollingHash<R> &t) {
    Int ok = 0, ng = std::min(s.size(), t.size()) + 1;
    while (ng - ok > 1) {
        Int mid = (ok + ng) / 2;
        if (s.gethash(0, mid) == t.gethash(0, mid))
            ok = mid;
        else
            ng = mid;
    }
    return ok;
}

// ハッシュを用いて辞書順に比較する。計算量はLCPと同じ。
template <bool R> int cmp(const BasicRollingHash<R> &s, const BasicRollingHash<R> &t) {
    Int n = LCP(s, t);
    if (n == std::min(s.size(), t.size()))
        return (s.size() > t.size()) - (s.size() < t.size());
    return static_cast<unsigned char>(s[n]) < static_cast<unsigned char>(t[n]) ? -1 : 1;
}

// ハッシュを用いて辞書順に比較する。計算量はLCPと同じ。
template <bool R> bool operator<(const BasicRollingHash<R> &s, const BasicRollingHash<R> &t) {
    return cmp(s, t) < 0;
}

// ビューの文字列をコピーする。長さNとして O(N+1)。
template <bool R> std::string to_string(const BasicRollingHash<R> &s) {
    return s.to_string();
}

// 接頭辞のハッシュを取り除く。計算量はBasicHashString::removePrefixと同じ。
inline auto removePrefix(const BasicHashString<false> &s, const BasicHashString<false> &p) {
    return s.removePrefix(p);
}

// 正逆ハッシュを比べて回文か判定する。計算量 O(1)。
inline bool isPalindrome(const BasicHashString<true> &s) {
    return s.isPalindrome();
}

// 反転した列のハッシュを返す。計算量 O(1)。
inline auto reversed(const BasicHashString<true> &s) {
    return s.reversed();
}

// 正逆ハッシュを比べて回文か判定する。計算量はHash変換と同じ。
inline bool isPalindrome(const BasicRollingHash<true> &s) {
    return s.isPalindrome();
}

// 反転した列のハッシュを返す。計算量はHash変換と同じ。
inline auto reversed(const BasicRollingHash<true> &s) {
    return BasicHashString<true>(s).reversed();
}
}
