#pragma once
#include <cplib/common.hpp>
#include <functional>

namespace cplib {
namespace detail {
enum class OpenHashState { empty, active, inactive };

inline std::size_t open_hash_capacity(std::size_t n) {
    return n ? std::size_t(1) << (std::bit_width(n) + 1) : 4;
}
}

// 線形探索と削除済み印による集合。再挿入時の fill 加算も元実装と同じ。
template <class T, class Hash = std::hash<T>> class HashSet {
    std::size_t length_ = 0, fill_ = 0, mask_ = 3;

    std::size_t find(const T &x) const {
        auto p = Hash{}(x)&mask_;
        while (values[p].state != detail::OpenHashState::empty && !(values[p].value == x))
            p = (p + 1) & mask_;
        return p;
    }

    void add_item(const T &x) {
        auto p = find(x);
        if (values[p].state == detail::OpenHashState::active)
            return;
        ++length_;
        ++fill_;
        values[p] = {detail::OpenHashState::active, x};
    }

    void resize() {
        auto old = std::move(values);
        values.resize(detail::open_hash_capacity(length_));
        mask_ = values.size() - 1;
        length_ = fill_ = 0;
        for (const auto &n : old)
            if (n.state == detail::OpenHashState::active)
                add_item(n.value);
    }

public:
    struct Node {
        detail::OpenHashState state = detail::OpenHashState::empty;
        T value{};
    };

    std::vector<Node> values = std::vector<Node>(4);

    // 登録済みの要素数を返す。O(1)。
    Int size() const {
        return length_;
    }

    // xを追加する。一様なハッシュを仮定して期待償却O(1)。容量Cで再構築を含む最悪O(C²)。
    void incl(const T &x) {
        add_item(x);
        if (detail::open_hash_capacity(fill_) > values.size())
            resize();
    }

    // xが登録済みか返す。一様なハッシュを仮定して期待O(1)、容量Cで最悪O(C)。
    bool contains(const T &x) const {
        return values[find(x)].state == detail::OpenHashState::active;
    }

    // xを削除する。一様なハッシュを仮定して期待O(1)、容量Cで最悪O(C)。
    void excl(const T &x) {
        auto p = find(x);
        if (values[p].state != detail::OpenHashState::active)
            return;
        --length_;
        values[p].state = detail::OpenHashState::inactive;
    }

    struct Iterator {
        const HashSet *owner;
        std::size_t pos;

        // 次の有効スロットまで進む。容量Cとして最悪O(C)、全列挙で合計O(C)。
        void skip() {
            while (pos < owner->values.size() &&
                   owner->values[pos].state != detail::OpenHashState::active)
                ++pos;
        }

        // 現在の要素への参照を返す。O(1)。
        const T &operator*() const {
            return owner->values[pos].value;
        }

        // 次の要素へ進む。容量Cとして最悪O(C)、全列挙で合計O(C)。
        Iterator &operator++() {
            ++pos;
            skip();
            return *this;
        }

        // 同じ位置を指すか比較する。O(1)。
        bool operator==(const Iterator &) const = default;
    };

    // 先頭イテレータを返す。容量Cとして最悪O(C)。
    Iterator begin() const {
        Iterator it{this, 0};
        it.skip();
        return it;
    }

    // 末尾の次のイテレータを返す。O(1)。
    Iterator end() const {
        return {this, values.size()};
    }

    // 列挙用に自身への参照を返す。O(1)。
    const HashSet &items() const {
        return *this;
    }
};

// 空のハッシュ集合を構築する。O(1)。
template <class T> HashSet<T> initHashSet() {
    return {};
}

// 配列からハッシュ集合を構築する。入力長nとして期待O(n)、最悪O(n²)。
template <class T> HashSet<T> toHashSet(const std::vector<T> &v) {
    HashSet<T> out;
    for (const T &x : v)
        out.incl(x);
    return out;
}
}
