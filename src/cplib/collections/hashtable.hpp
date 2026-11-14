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

    // 指定容量の空ハッシュ表を構築する。O(capacity+1)。
    explicit HashTable(Int capacity = 0)
        : mask_(detail::open_hash_capacity(capacity) - 1), storage(mask_ + 1) {
        assert(capacity >= 0);
    }

    // 登録済みの要素数を返す。O(1)。
    Int size() const {
        return length_;
    }

    // kが登録済みか返す。一様なハッシュを仮定して期待O(1)、容量Cで最悪O(C)。
    bool contains(const K &k) const {
        return storage[find(k)].state == detail::OpenHashState::active;
    }

    // kが登録済みか返す。一様なハッシュを仮定して期待O(1)、容量Cで最悪O(C)。
    bool hasKey(const K &k) const {
        return contains(k);
    }

    // キーと値の組を登録または更新する。一様なハッシュを仮定して期待償却O(1)、容量Cで再構築を含む最悪O(C²)。
    void incl(const std::pair<K, V> &p) {
        set(p.first, p.second);
    }

    // kの値を登録または更新する。一様なハッシュを仮定して期待償却O(1)、容量Cで再構築を含む最悪O(C²)。
    void set(const K &k, const V &v) {
        add_item(k, v);
        if (detail::open_hash_capacity(fill_) > storage.size())
            resize();
    }

    // 登録済みkの値への読み取り専用参照を返す。一様なハッシュを仮定して期待O(1)、容量Cで最悪O(C)。
    const V &at(const K &k) const {
        auto p = find(k);
        assert(storage[p].state == detail::OpenHashState::active);
        return storage[p].value.second;
    }

    // 登録済みkの値への変更可能な参照を返す。一様なハッシュを仮定して期待O(1)、容量Cで最悪O(C)。
    V &at(const K &k) {
        auto p = find(k);
        assert(storage[p].state == detail::OpenHashState::active);
        return storage[p].value.second;
    }

    struct Reference {
        HashTable *owner;
        K key;

        // 登録済みキーの値への参照を返す。一様なハッシュを仮定して期待O(1)、容量Cで最悪O(C)。
        operator V &() const {
            return owner->at(key);
        }

        // 対象キーの値を登録または更新する。一様なハッシュを仮定して期待償却O(1)、容量Cで再構築を含む最悪O(C²)。
        Reference &operator=(const V &v) {
            owner->set(key, v);
            return *this;
        }

        // 別の登録済みキーから値を読み取り代入する。一様なハッシュを仮定して期待償却O(1)、容量Cで再構築を含む最悪O(C²)。
        Reference &operator=(const Reference &other) {
            V v = other.owner->at(other.key);
            return *this = v;
        }

        // 登録済みキーの値にvを加算する。一様なハッシュを仮定して期待O(1)、容量Cで最悪O(C)。
        template <class U> Reference &operator+=(const U &v) {
            owner->at(key) += v;
            return *this;
        }

        // 登録済みキーの値からvを減算する。一様なハッシュを仮定して期待O(1)、容量Cで最悪O(C)。
        template <class U> Reference &operator-=(const U &v) {
            owner->at(key) -= v;
            return *this;
        }

        // 登録済みキーの値にvを乗算する。一様なハッシュを仮定して期待O(1)、容量Cで最悪O(C)。
        template <class U> Reference &operator*=(const U &v) {
            owner->at(key) *= v;
            return *this;
        }

        // 登録済みキーの値をvで除算する。一様なハッシュを仮定して期待O(1)、容量Cで最悪O(C)。
        template <class U> Reference &operator/=(const U &v) {
            owner->at(key) /= v;
            return *this;
        }
    };

    // 登録済みkの値への読み取り専用参照を返す。一様なハッシュを仮定して期待O(1)、容量Cで最悪O(C)。
    const V &operator[](const K &k) const {
        return at(k);
    }

    // kに対する読み書き参照を返す。O(1)とキーのコピーに要する時間。
    Reference operator[](const K &k) {
        return {this, k};
    }

    // 全要素を削除して初期容量に戻す。容量CとしてO(C)。
    void clear() {
        *this = HashTable();
    }

    // kを削除する。一様なハッシュを仮定して期待O(1)、容量Cで最悪O(C)。
    void del(const K &k) {
        auto p = find(k);
        if (storage[p].state != detail::OpenHashState::active)
            return;
        --length_;
        storage[p].state = detail::OpenHashState::inactive;
    }

    // kを削除する。一様なハッシュを仮定して期待O(1)、容量Cで最悪O(C)。
    void excl(const K &k) {
        del(k);
    }

    struct Iterator {
        const HashTable *owner;
        std::size_t pos;

        // 次の有効スロットまで進む。容量Cとして最悪O(C)、全列挙で合計O(C)。
        void skip() {
            while (pos < owner->storage.size() &&
                   owner->storage[pos].state != detail::OpenHashState::active)
                ++pos;
        }

        // 現在のキーと値への参照を返す。O(1)。
        const std::pair<K, V> &operator*() const {
            return owner->storage[pos].value;
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

    struct PairRange {
        const HashTable *owner;

        // 先頭イテレータを返す。容量Cとして最悪O(C)。
        Iterator begin() const {
            Iterator it{owner, 0};
            it.skip();
            return it;
        }

        // 末尾の次のイテレータを返す。O(1)。
        Iterator end() const {
            return {owner, owner->storage.size()};
        }
    };

    // 登録済みの組を列挙するビューを返す。生成はO(1)、全列挙は容量CとしてO(C)。
    PairRange pairs() const {
        return {this};
    }

    // 登録済みキーを列挙するビューを返す。生成はO(1)、全列挙は容量CとしてO(C)。
    auto keys() const {
        return detail::transform_range(pairs(), [](const auto &p) { return p.first; });
    }

    // 登録済み値を列挙するビューを返す。生成はO(1)、全列挙は容量CとしてO(C)。
    auto values() const {
        return detail::transform_range(pairs(), [](const auto &p) { return p.second; });
    }

    // 登録済みの組からハッシュを計算する。容量CとしてO(C)と各組のハッシュ計算に要する時間。
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

// 指定容量の空ハッシュ表を構築する。O(capacity+1)。
template <class K, class V> auto initHashTable(Int capacity = 0) {
    return HashTable<K, V>(capacity);
}

// 登録済みの組からハッシュを計算する。容量CとしてO(C)と各組のハッシュ計算に要する時間。
template <class K, class V, class H> auto hashValue(const HashTable<K, V, H> &t) {
    return t.hash();
}
}
