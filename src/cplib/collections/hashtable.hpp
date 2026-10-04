#pragma once
#include <cplib/collections/hashset.hpp>
#include <cplib/iterator_views.hpp>

namespace cplib {
// 線形探索ハッシュ表。参照では存在を要求し、代入では新規キーを挿入する。
template <class K, class V, class Hash = std::hash<K>> class HashTable {
    std::size_t length_ = 0, fill_ = 0, mask_;

    std::size_t find(const K &x) const {
        auto p = Hash{}(x)&mask_;
        while (storage[p].state != detail::OpenHashState::empty && !(storage[p].value.first == x))
            p = (p + 1) & mask_;
        return p;
    }

    void add_item(const K &k, const V &v) {
        auto p = find(k);
        if (storage[p].state == detail::OpenHashState::active) {
            storage[p].value.second = v;
            return;
        }
        ++length_;
        ++fill_;
        storage[p] = {detail::OpenHashState::active, {k, v}};
    }

    void resize() {
        auto old = std::move(storage);
        storage.resize(detail::open_hash_capacity(length_));
        mask_ = storage.size() - 1;
        length_ = fill_ = 0;
        for (const auto &n : old)
            if (n.state == detail::OpenHashState::active)
                add_item(n.value.first, n.value.second);
    }

public:
    struct Node {
        detail::OpenHashState state = detail::OpenHashState::empty;
        std::pair<K, V> value{};
    };

    std::vector<Node> storage; // Nim の values フィールド。values() と衝突するため改名。

    explicit HashTable(Int capacity = 0)
        : mask_(detail::open_hash_capacity(capacity) - 1), storage(mask_ + 1) {
        assert(capacity >= 0);
    }

    Int len() const {
        return length_;
    }

    bool contains(const K &k) const {
        return storage[find(k)].state == detail::OpenHashState::active;
    }

    bool hasKey(const K &k) const {
        return contains(k);
    }

    void incl(const std::pair<K, V> &p) {
        set(p.first, p.second);
    }

    void set(const K &k, const V &v) {
        add_item(k, v);
        if (detail::open_hash_capacity(fill_) > storage.size())
            resize();
    }

    const V &at(const K &k) const {
        auto p = find(k);
        assert(storage[p].state == detail::OpenHashState::active);
        return storage[p].value.second;
    }

    V &at(const K &k) {
        auto p = find(k);
        assert(storage[p].state == detail::OpenHashState::active);
        return storage[p].value.second;
    }

    struct Reference {
        HashTable *owner;
        K key;

        operator V &() const {
            return owner->at(key);
        }

        Reference &operator=(const V &v) {
            owner->set(key, v);
            return *this;
        }

        Reference &operator=(const Reference &other) {
            V v = other.owner->at(other.key);
            return *this = v;
        }

        template <class U> Reference &operator+=(const U &v) {
            owner->at(key) += v;
            return *this;
        }

        template <class U> Reference &operator-=(const U &v) {
            owner->at(key) -= v;
            return *this;
        }

        template <class U> Reference &operator*=(const U &v) {
            owner->at(key) *= v;
            return *this;
        }

        template <class U> Reference &operator/=(const U &v) {
            owner->at(key) /= v;
            return *this;
        }
    };

    const V &operator[](const K &k) const {
        return at(k);
    }

    Reference operator[](const K &k) {
        return {this, k};
    }

    void clear() {
        *this = HashTable();
    }

    void del(const K &k) {
        auto p = find(k);
        if (storage[p].state != detail::OpenHashState::active)
            return;
        --length_;
        storage[p].state = detail::OpenHashState::inactive;
    }

    void excl(const K &k) {
        del(k);
    }

    struct Iterator {
        const HashTable *owner;
        std::size_t pos;

        void skip() {
            while (pos < owner->storage.size() &&
                   owner->storage[pos].state != detail::OpenHashState::active)
                ++pos;
        }

        const std::pair<K, V> &operator*() const {
            return owner->storage[pos].value;
        }

        Iterator &operator++() {
            ++pos;
            skip();
            return *this;
        }

        bool operator==(const Iterator &) const = default;
    };

    struct PairRange {
        const HashTable *owner;

        Iterator begin() const {
            Iterator it{owner, 0};
            it.skip();
            return it;
        }

        Iterator end() const {
            return {owner, owner->storage.size()};
        }
    };

    PairRange pairs() const {
        return {this};
    }

    auto keys() const {
        return detail::transform_range(pairs(), [](const auto &p) { return p.first; });
    }

    auto values() const {
        return detail::transform_range(pairs(), [](const auto &p) { return p.second; });
    }

    std::size_t hash() const {
        std::size_t result = 0;
        for (const auto &[k, v] : pairs()) {
            auto a = Hash{}(k), b = std::hash<V>{}(v);
            auto h = a ^ (b + 0x9e3779b97f4a7c15ULL + (a << 6) + (a >> 2));
            result += h;
            result += result << 10;
            result ^= result >> 6;
        }
        return result;
    }
};

template <class K, class V> auto initHashTable(Int capacity = 0) {
    return HashTable<K, V>(capacity);
}

template <class K, class V, class H> auto hash(const HashTable<K, V, H> &t) {
    return t.hash();
}
}
