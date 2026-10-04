#pragma once
#include <cplib/common.hpp>
#include <cplib/utils/backwards_index.hpp>

namespace cplib {
template <class T> class FenwickTree {
    Int size_;
    std::vector<T> data;

    static Int slot(Int i) {
        return i + (i >> 10);
    }

public:
    // 1024要素ごとに余白を設けて構築する。O(n)。
    // 長さnの零配列から構築します。O(n)時間・領域です。
    // T{}を零とし、Tには可換加法群の演算+=と減算が必要です。
    explicit FenwickTree(Int n) : size_(n), data(slot(n) + 1) {
        assert(n >= 0);
    }

    explicit FenwickTree(std::span<const T> values) : FenwickTree(values.size()) {
        for (Int i = 1; i <= size_; ++i)
            data[slot(i)] = values[i - 1];
        for (Int i = 1; i <= size_; ++i) {
            Int parent = i + (i & -i);
            if (parent <= size_)
                data[slot(parent)] += data[slot(i)];
        }
    }

    Int len() const {
        return size_;
    }

    Int size() const {
        return len();
    }

    // 一点加算。O(log n)。
    void add(Int p, const T &delta) {
        assert(0 <= p && p < size_);
        for (Int i = p + 1; i <= size_; i += i & -i)
            data[slot(i)] += delta;
    }

    // [0, r)の和をO(log n)で返します。
    T prefix(Int r) const {
        assert(0 <= r && r <= size_);
        T result{};
        while (r > 0) {
            result += data[slot(r)];
            r &= r - 1;
        }
        return result;
    }

    // 共通祖先の走査を省いて半開区間の和を求める。O(log n)。
    T get(Int l, Int r) const {
        assert(0 <= l && l <= r && r <= size_);
        T result{}, left{};
        while (r > l) {
            result += data[slot(r)];
            r &= r - 1;
        }
        while (l > r) {
            left += data[slot(l)];
            l &= l - 1;
        }
        return result - left;
    }

    T operator[](Int p) const {
        return get(p, p + 1);
    }

    void set(Int p, const T &value) {
        add(p, value - get(p, p + 1));
    }

    struct Reference {
        FenwickTree *tree;
        Int p;

        operator T() const {
            return tree->get(p, p + 1);
        }

        Reference &operator=(const T &value) {
            tree->set(p, value);
            return *this;
        }

        Reference &operator=(const Reference &other) {
            return *this = T(other);
        }

        Reference &operator+=(const T &delta) {
            tree->add(p, delta);
            return *this;
        }

        Reference &operator-=(const T &delta) {
            tree->add(p, -delta);
            return *this;
        }
    };

    Reference operator[](Int p) {
        return {this, p};
    }

    CPLIB_BACKWARDS_INDEX_OVERLOADS
    template <class L, class R> T operator[](ClosedSlice<L, R> segment) const {
        return get(resolve_index(len(), segment.a), resolve_index(len(), segment.b) + 1);
    }
};

template <class T> FenwickTree<T> initFenwickTree(Int n) {
    return FenwickTree<T>(n);
}

template <class T> FenwickTree<T> initFenwickTree(std::span<const T> a) {
    return FenwickTree<T>(a);
}

template <class T> FenwickTree<T> initFenwickTree(const std::vector<T> &a) {
    return FenwickTree<T>(a);
}

template <class T> Int len(const FenwickTree<T> &f) {
    return f.len();
}

template <class T> void add(FenwickTree<T> &f, Int p, const T &delta) {
    f.add(p, delta);
}

template <class T> T prefix(const FenwickTree<T> &f, Int r) {
    return f.prefix(r);
}

template <class T> T get(const FenwickTree<T> &f, Int l, Int r) {
    return f.get(l, r);
}
}
