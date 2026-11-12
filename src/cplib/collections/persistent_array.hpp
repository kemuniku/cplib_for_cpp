#pragma once
#include <cplib/common.hpp>
#include <bit>
#include <memory>

namespace cplib {
template <int Shift, class T> struct PersistentArrayNode {
    std::vector<std::shared_ptr<PersistentArrayNode>> arr;
    T value{};
};

// 2^Shift分木による永続配列。取得O(log N)、変更O(2^Shift log N)で経路をコピー。
template <int Shift, class T> class PersistentArray {
    static_assert(Shift > 0 && Shift < 31);
    using Node = PersistentArrayNode<Shift, T>;
    using Ptr = std::shared_ptr<Node>;
    Int size_ = 0, h = 0;
    Ptr root;

    PersistentArray(Int n, Int height, Ptr r) : size_(n), h(height), root(std::move(r)) {
    }

public:
    // 配列から永続配列を構築する。要素数をn、木の高さをhとしてO(n+h+1)。
    explicit PersistentArray(std::span<const T> values)
        : size_(values.size()), h((std::bit_width(UInt(size_)) + Shift - 1) / Shift),
          root(std::make_shared<Node>()) {
        auto dfs = [&](auto &&self, Ptr node, Int now, Int depth) -> bool {
            if (depth == h) {
                if (now < size_)
                    node->value = values[now];
                else
                    return false;
            } else
                for (Int i = 0; i < (Int(1) << Shift); ++i) {
                    node->arr.push_back(std::make_shared<Node>());
                    if (!self(self, node->arr[i], (now << Shift) | i, depth + 1))
                        return false;
                }
            return true;
        };
        dfs(dfs, root, 0, 0);
    }

    // 要素数を返す。O(1)。
    Int size() const {
        return size_;
    }

    // この版の全要素を配列にして返す。要素数をn、木の高さをhとしてO(n+h+1)。
    std::vector<T> toseq() const {
        std::vector<T> result;
        result.reserve(size_);
        auto dfs = [&](auto &&self, Ptr node, Int now, Int depth) -> bool {
            if (depth == h) {
                if (now < size_)
                    result.push_back(node->value);
                else
                    return false;
            } else
                for (Int i = 0; i < (Int(1) << Shift); ++i)
                    if (!self(self, node->arr[i], (now << Shift) | i, depth + 1))
                        return false;
            return true;
        };
        dfs(dfs, root, 0, 0);
        return result;
    }

    // この版の全要素を配列にして返す。要素数をn、木の高さをhとしてO(n+h+1)。
    std::vector<T> toSeq() const {
        return toseq();
    }

    // 指定位置の値を返す。木の高さをhとしてO(h+1)。
    T operator[](Int index) const {
        assert(0 <= index && index < size_);
        std::vector<Int> indices(h);
        for (Int i = h - 1; i >= 0; --i) {
            indices[i] = index & ((Int(1) << Shift) - 1);
            index >>= Shift;
        }
        auto now = root;
        for (Int i : indices)
            now = now->arr[i];
        return now->value;
    }

    // 指定位置を変更した新版を返す。木の高さをh、分岐数をBとしてO(Bh+1)。
    PersistentArray change_value(Int index, const T &value) const {
        assert(0 <= index && index < size_);
        std::vector<Int> indices(h);
        for (Int i = h - 1; i >= 0; --i) {
            indices[i] = index & ((Int(1) << Shift) - 1);
            index >>= Shift;
        }
        auto now = root;
        std::vector<Ptr> stack;
        stack.reserve(h);
        for (Int i : indices) {
            stack.push_back(now);
            now = now->arr[i];
        }
        now = std::make_shared<Node>();
        now->value = value;
        for (Int i = h - 1; i >= 0; --i) {
            auto tmp = std::make_shared<Node>();
            tmp->arr = stack[i]->arr;
            tmp->arr[indices[i]] = now;
            now = tmp;
        }
        return PersistentArray(size_, h, now);
    }
};

// 配列から永続配列を構築する。要素数をn、木の高さをhとしてO(n+h+1)。
template <class T, int Shift = 5> auto initPersistentArray(std::span<const T> v) {
    return PersistentArray<Shift, T>(v);
}

// 配列から永続配列を構築する。要素数をn、木の高さをhとしてO(n+h+1)。
template <class T, int Shift = 5> auto initPersistentArray(const std::vector<T> &v) {
    return PersistentArray<Shift, T>(v);
}

// この版の全要素を配列にして返す。要素数をn、木の高さをhとしてO(n+h+1)。
template <int Shift, class T> auto toseq(const PersistentArray<Shift, T> &p) {
    return p.toseq();
}

// この版の全要素を配列にして返す。要素数をn、木の高さをhとしてO(n+h+1)。
template <int Shift, class T> auto toSeq(const PersistentArray<Shift, T> &p) {
    return p.toseq();
}

// 位置iをvに変更した新版を返す。木の高さをh、分岐数をBとしてO(Bh+1)。
template <int Shift, class T>
auto change_value(const PersistentArray<Shift, T> &p, Int i, const T &v) {
    return p.change_value(i, v);
}
}
