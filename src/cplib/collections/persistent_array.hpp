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
    Int size = 0, h = 0;
    Ptr root;

    PersistentArray(Int n, Int height, Ptr r) : size(n), h(height), root(std::move(r)) {
    }

public:
    explicit PersistentArray(std::span<const T> values)
        : size(values.size()), h((std::bit_width(UInt(size)) + Shift - 1) / Shift),
          root(std::make_shared<Node>()) {
        auto dfs = [&](auto &&self, Ptr node, Int now, Int depth) -> bool {
            if (depth == h) {
                if (now < size)
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

    std::vector<T> toseq() const {
        std::vector<T> result;
        result.reserve(size);
        auto dfs = [&](auto &&self, Ptr node, Int now, Int depth) -> bool {
            if (depth == h) {
                if (now < size)
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

    std::vector<T> toSeq() const {
        return toseq();
    }

    T operator[](Int index) const {
        assert(0 <= index && index < size);
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

    PersistentArray change_value(Int index, const T &value) const {
        assert(0 <= index && index < size);
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
        return PersistentArray(size, h, now);
    }
};

template <class T, int Shift = 5> auto initPersistentArray(std::span<const T> v) {
    return PersistentArray<Shift, T>(v);
}

template <class T, int Shift = 5> auto initPersistentArray(const std::vector<T> &v) {
    return PersistentArray<Shift, T>(v);
}

template <int Shift, class T> auto toseq(const PersistentArray<Shift, T> &p) {
    return p.toseq();
}

template <int Shift, class T> auto toSeq(const PersistentArray<Shift, T> &p) {
    return p.toseq();
}

template <int Shift, class T>
auto change_value(const PersistentArray<Shift, T> &p, Int i, const T &v) {
    return p.change_value(i, v);
}
}
