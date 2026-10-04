#pragma once
#include <cplib/common.hpp>
#include <unordered_map>
#include <cplib/iterator_views.hpp>
#include <sstream>

namespace cplib {
template <class K, class V> class DefaultDict {
    std::unordered_map<K, V> table;
    V defaultValue;

    struct TableRange {
        const std::unordered_map<K, V> *table;

        auto begin() const {
            return table->begin();
        }

        auto end() const {
            return table->end();
        }
    };

public:
    explicit DefaultDict(V value) : defaultValue(std::move(value)) {
    }

    DefaultDict(std::unordered_map<K, V> table, V value)
        : table(std::move(table)), defaultValue(std::move(value)) {
    }

    // constでの参照は挿入せず、非constでは既定値を挿入する。期待O(1)。
    V operator[](const K &key) const {
        auto it = table.find(key);
        return it == table.end() ? defaultValue : it->second;
    }

    V &operator[](const K &key) {
        return table.try_emplace(key, defaultValue).first->second;
    }

    bool operator==(const DefaultDict &other) const {
        return table == other.table;
    }

    void clear() {
        table = std::unordered_map<K, V>{};
    }

    bool contains(const K &key) const {
        return table.contains(key);
    }

    bool hasKey(const K &key) const {
        return contains(key);
    }

    void del(const K &key) {
        table.erase(key);
    }

    Int len() const {
        return table.size();
    }

    Int size() const {
        return len();
    }

    bool pop(const K &key, V &value) {
        auto it = table.find(key);
        if (it == table.end())
            return false;
        value = it->second;
        table.erase(it);
        return true;
    }

    bool take(const K &key, V &value) {
        return pop(key, value);
    }

    const auto &pairs() const {
        return table;
    }

    auto &mpairs() {
        return table;
    }

    auto keys() const {
        return detail::transform_range(TableRange{&table}, [](const auto &p) { return p.first; });
    }

    auto values() const {
        return detail::transform_range(TableRange{&table}, [](const auto &p) { return p.second; });
    }

    std::size_t hash() const {
        std::size_t out = 0;
        for (const auto &[k, v] : table) {
            std::size_t a = std::hash<K>{}(k), b = std::hash<V>{}(v);
            out += a ^ (b + 0x9e3779b97f4a7c15ULL + (a << 6) + (a >> 2));
        }
        return out;
    }

    friend std::ostream &operator<<(std::ostream &out, const DefaultDict &d) {
        out << '{';
        bool first = true;
        for (const auto &[k, v] : d.table) {
            if (!first)
                out << ", ";
            first = false;
            out << k << ": " << v;
        }
        return out << '}';
    }
};

template <class K, class V> auto initDefaultDict(V value) {
    return DefaultDict<K, V>(std::move(value));
}

template <class K, class V> auto toDefaultDict(std::span<const std::pair<K, V>> pairs, V value) {
    std::unordered_map<K, V> table;
    for (const auto &[k, v] : pairs)
        table[k] = v;
    return DefaultDict<K, V>(std::move(table), std::move(value));
}

template <class K, class V> auto toDefaultDict(const std::vector<std::pair<K, V>> &pairs, V value) {
    return toDefaultDict<K, V>(std::span<const std::pair<K, V>>(pairs), std::move(value));
}

template <class K, class V> auto toDefaultDict(std::unordered_map<K, V> table, V value) {
    return DefaultDict<K, V>(std::move(table), std::move(value));
}

template <class K, class V> Int len(const DefaultDict<K, V> &d) {
    return d.len();
}

template <class K, class V> bool hasKey(const DefaultDict<K, V> &d, const K &k) {
    return d.hasKey(k);
}

template <class K, class V> bool contains(const DefaultDict<K, V> &d, const K &k) {
    return d.contains(k);
}

template <class K, class V> void del(DefaultDict<K, V> &d, const K &k) {
    d.del(k);
}

template <class K, class V> void clear(DefaultDict<K, V> &d) {
    d.clear();
}

template <class K, class V> auto hash(const DefaultDict<K, V> &d) {
    return d.hash();
}

template <class K, class V> bool pop(DefaultDict<K, V> &d, const K &k, V &v) {
    return d.pop(k, v);
}

template <class K, class V> bool take(DefaultDict<K, V> &d, const K &k, V &v) {
    return d.take(k, v);
}

template <class K, class V> const auto &pairs(const DefaultDict<K, V> &d) {
    return d.pairs();
}

template <class K, class V> auto &mpairs(DefaultDict<K, V> &d) {
    return d.mpairs();
}

template <class K, class V> auto keys(const DefaultDict<K, V> &d) {
    return d.keys();
}

template <class K, class V> auto values(const DefaultDict<K, V> &d) {
    return d.values();
}

template <class K, class V> auto to_string(const DefaultDict<K, V> &d) {
    std::ostringstream out;
    out << d;
    return out.str();
}
}

namespace std {
template <class K, class V> struct hash<cplib::DefaultDict<K, V>> {
    size_t operator()(const cplib::DefaultDict<K, V> &d) const {
        return d.hash();
    }
};
}
