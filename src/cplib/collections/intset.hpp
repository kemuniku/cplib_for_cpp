#pragma once
#include <cplib/common.hpp>
#include <cplib/utils/backwards_index.hpp>
#include <stdexcept>
#include <iterator>

namespace cplib {
// Bit blocks linked through the nonempty blocks: worst-case O(1) updates.
class IntSet {
    struct Block {
        UInt bits = 0;
        Int prev = 0, next = 0;
    };

    std::vector<Block> blocks;
    Int lower = 0, size_ = 0, count_ = 0, head = 0;

    Int offset(Int x) const {
        if (x < lower)
            return -1;
        UInt d = UInt(x) - UInt(lower);
        return d < UInt(size_) ? Int(d) : -1;
    }

public:
    explicit IntSet(Int n = 0) : size_(n) {
        if (n < 0)
            throw std::invalid_argument("negative IntSet size");
        blocks.resize((UInt(n) + 63) / 64);
    }

    template <class L, class R> explicit IntSet(ClosedSlice<L, R> b) {
        lower = b.a;
        if (b.a <= b.b) {
            UInt d = UInt(b.b) - UInt(b.a);
            if (d >= UInt(std::numeric_limits<Int>::max()))
                throw std::invalid_argument("IntSet range too wide");
            size_ = Int(d) + 1;
            blocks.resize((UInt(size_) + 63) / 64);
        }
    }

    Int len() const {
        return count_;
    }

    Int size() const {
        return len();
    }

    // xが含まれるかを返します。範囲外はfalseです。O(1)。
    bool contains(Int x) const {
        Int p = offset(x);
        return p >= 0 && (blocks[p >> 6].bits & (UInt(1) << (p & 63)));
    }

    // xを挿入し、既に含まれていたかを返す。範囲外はstd::out_of_range。O(1)。
    bool containsOrIncl(Int x) {
        Int p = offset(x);
        if (p < 0)
            throw std::out_of_range("IntSet value");
        Int i = p >> 6;
        UInt m = UInt(1) << (p & 63);
        auto &b = blocks[i];
        if (b.bits & m)
            return true;
        if (!b.bits) {
            b.prev = 0;
            b.next = head;
            if (head)
                blocks[head - 1].prev = i + 1;
            head = i + 1;
        }
        b.bits |= m;
        ++count_;
        return false;
    }

    void incl(Int x) {
        containsOrIncl(x);
    }

    // xを削除し、元々含まれていなかったかを返します。範囲外もtrueです。O(1)。
    bool missingOrExcl(Int x) {
        Int p = offset(x);
        if (p < 0)
            return true;
        auto &b = blocks[p >> 6];
        UInt m = UInt(1) << (p & 63);
        if (!(b.bits & m))
            return true;
        b.bits &= ~m;
        --count_;
        if (!b.bits) {
            if (!b.prev)
                head = b.next;
            else
                blocks[b.prev - 1].next = b.next;
            if (b.next)
                blocks[b.next - 1].prev = b.prev;
        }
        return false;
    }

    void excl(Int x) {
        missingOrExcl(x);
    }

    // 任意の要素を一つ取り出して削除する。空集合はstd::out_of_range。O(1)。
    Int pop() {
        if (!count_)
            throw std::out_of_range("empty IntSet");
        Int x = lower + (head - 1) * 64 + std::countr_zero(blocks[head - 1].bits);
        excl(x);
        return x;
    }

    // 範囲と確保済み領域を保って空にします。非空ブロック数をBとしてO(B + 1)。
    void clear() {
        while (head) {
            Int i = head - 1;
            head = blocks[i].next;
            blocks[i] = {};
        }
        count_ = 0;
    }

    struct Iterator {
        const IntSet *set;
        Int node;
        UInt bits;

        Int operator*() const {
            return set->lower + (node - 1) * 64 + std::countr_zero(bits);
        }

        Iterator &operator++() {
            bits &= bits - 1;
            if (!bits) {
                node = set->blocks[node - 1].next;
                if (node)
                    bits = set->blocks[node - 1].bits;
            }
            return *this;
        }

        bool operator==(std::default_sentinel_t) const {
            return node == 0;
        }
    };

    // 要素を任意順に列挙する。列挙中の変更は不可。要素数kとして全体O(k + 1)、一要素あたりO(1)。
    Iterator begin() const {
        return {this, head, head ? blocks[head - 1].bits : 0};
    }

    std::default_sentinel_t end() const {
        return {};
    }

    // 範囲や挿入順によらず要素が等しいかを返す。O(a.len() + 1)。
    friend bool operator==(const IntSet &a, const IntSet &b) {
        if (a.len() != b.len())
            return false;
        for (Int x : a)
            if (!b.contains(x))
                return false;
        return true;
    }

    std::string str() const {
        std::string s = "{";
        for (Int x : *this) {
            if (s.size() > 1)
                s += ", ";
            s += std::to_string(x);
        }
        return s + "}";
    }
};

inline IntSet initIntSet(Int n) {
    return IntSet(n);
}

template <class L, class R> IntSet initIntSet(ClosedSlice<L, R> b) {
    return IntSet(b);
}

template <class R, class B> IntSet toIntSet(const R &v, B bounds) {
    IntSet s(bounds);
    for (Int x : v)
        s.incl(x);
    return s;
}
}
