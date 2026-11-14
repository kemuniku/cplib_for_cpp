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
    // 既定値を指定して空の辞書を構築する。O(1)と値のムーブに要する時間。
    explicit DefaultDict(V value) : defaultValue(std::move(value)) {
    }

    // 辞書と既定値を受け取る。引数の受け渡し後はO(1)と値のムーブに要する時間。
    DefaultDict(std::unordered_map<K, V> table, V value)
        : table(std::move(table)), defaultValue(std::move(value)) {
    }

    // キーを挿入せず登録値または既定値を返す。期待 O(1)、最悪 O(N) と値のコピー費用。
    V get(const K &key) const {
        auto it = table.find(key);
        return it == table.end() ? defaultValue : it->second;
    }

    // constでの参照は挿入せず、非constでは既定値を挿入する。期待O(1)。
    V operator[](const K &key) const {
        auto it = table.find(key);
        return it == table.end() ? defaultValue : it->second;
    }

    // keyの値への参照を返し、存在しなければ既定値を挿入する。期待償却O(1)、要素数nで最悪O(n)。
    V &operator[](const K &key) {
        return table.try_emplace(key, defaultValue).first->second;
    }

    // 保存されたキーと値の組が等しいか返す。要素数nとして期待O(n)、最悪O(n²)。既定値は比較しない。
    bool operator==(const DefaultDict &other) const {
        return table == other.table;
    }

    // 全要素を削除する。要素数nとバケット数bに対しO(n+b)。
    void clear() {
        table = std::unordered_map<K, V>{};
    }

    // keyが登録済みか返す。要素数nとして期待O(1)、最悪O(n)。
    bool contains(const K &key) const {
        return table.contains(key);
    }

    // keyが登録済みか返す。要素数nとして期待O(1)、最悪O(n)。
    bool hasKey(const K &key) const {
        return contains(key);
    }

    // keyを削除する。要素数nとして期待O(1)、最悪O(n)。
    void del(const K &key) {
        table.erase(key);
    }

    // 登録済みの要素数を返す。O(1)。
    Int size() const {
        return table.size();
    }

    // keyを削除して値をvalueに格納し、存在したか返す。要素数nとして期待O(1)、最悪O(n)。
    bool pop(const K &key, V &value) {
        auto it = table.find(key);
        if (it == table.end())
            return false;
        value = it->second;
        table.erase(it);
        return true;
    }

    // keyを削除して値をvalueに格納し、存在したか返す。要素数nとして期待O(1)、最悪O(n)。
    bool take(const K &key, V &value) {
        return pop(key, value);
    }

    // 登録済みのキーと値を参照する。O(1)。
    const auto &pairs() const {
        return table;
    }

    // 登録済みのキーと値を変更可能な参照で返す。O(1)。
    auto &mpairs() {
        return table;
    }

    // キーを列挙するビューを返す。生成はO(1)、全列挙は要素数nとしてO(n)。
    auto keys() const {
        return detail::transform_range(TableRange{&table}, [](const auto &p) { return p.first; });
    }

    // 値を列挙するビューを返す。生成はO(1)、全列挙は要素数nとしてO(n)。
    auto values() const {
        return detail::transform_range(TableRange{&table}, [](const auto &p) { return p.second; });
    }

    // 登録済みの組から順序によらないハッシュを計算する。要素数nとしてO(n)回のハッシュ計算。
    std::size_t hash() const {
        std::size_t out = 0;
        for (const auto &[k, v] : table) {
            std::size_t a = std::hash<K>{}(k), b = std::hash<V>{}(v);
            out += a ^ (b + 0x9e3779b97f4a7c15ULL + (a << 6) + (a >> 2));
        }
        return out;
    }

    // 登録済みの組をストリームに出力する。要素数nとしてO(n)と各組の出力に要する時間。
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

// 既定値を指定して空の辞書を生成する。O(1)と値のムーブに要する時間。
template <class K, class V> auto initDefaultDict(V value) {
    return DefaultDict<K, V>(std::move(value));
}

// キーと値の配列から辞書を構築する。入力長nとして期待O(n)、最悪O(n²)。
template <class K, class V> auto toDefaultDict(std::span<const std::pair<K, V>> pairs, V value) {
    std::unordered_map<K, V> table;
    for (const auto &[k, v] : pairs)
        table[k] = v;
    return DefaultDict<K, V>(std::move(table), std::move(value));
}

// キーと値の配列から辞書を構築する。入力長nとして期待O(n)、最悪O(n²)。
template <class K, class V> auto toDefaultDict(const std::vector<std::pair<K, V>> &pairs, V value) {
    return toDefaultDict<K, V>(std::span<const std::pair<K, V>>(pairs), std::move(value));
}

// 連想配列と既定値から辞書を構築する。引数の受け渡し後はO(1)と値のムーブに要する時間。
template <class K, class V> auto toDefaultDict(std::unordered_map<K, V> table, V value) {
    return DefaultDict<K, V>(std::move(table), std::move(value));
}

// kが登録済みか返す。要素数nとして期待O(1)、最悪O(n)。
template <class K, class V> bool hasKey(const DefaultDict<K, V> &d, const K &k) {
    return d.hasKey(k);
}

// kが登録済みか返す。要素数nとして期待O(1)、最悪O(n)。
template <class K, class V> bool contains(const DefaultDict<K, V> &d, const K &k) {
    return d.contains(k);
}

// kを削除する。要素数nとして期待O(1)、最悪O(n)。
template <class K, class V> void del(DefaultDict<K, V> &d, const K &k) {
    d.del(k);
}

// 全要素を削除する。要素数nとバケット数bに対しO(n+b)。
template <class K, class V> void clear(DefaultDict<K, V> &d) {
    d.clear();
}

// 登録済みの組のハッシュを返す。要素数nとしてO(n)回のハッシュ計算。
template <class K, class V> auto hashValue(const DefaultDict<K, V> &d) {
    return d.hash();
}

// kを削除し、値をvに格納して存在したか返す。要素数nとして期待O(1)、最悪O(n)。
template <class K, class V> bool pop(DefaultDict<K, V> &d, const K &k, V &v) {
    return d.pop(k, v);
}

// kを削除し、値をvに格納して存在したか返す。要素数nとして期待O(1)、最悪O(n)。
template <class K, class V> bool take(DefaultDict<K, V> &d, const K &k, V &v) {
    return d.take(k, v);
}

// 登録済みのキーと値を参照する。O(1)。
template <class K, class V> const auto &pairs(const DefaultDict<K, V> &d) {
    return d.pairs();
}

// 登録済みのキーと値を変更可能な参照で返す。O(1)。
template <class K, class V> auto &mpairs(DefaultDict<K, V> &d) {
    return d.mpairs();
}

// キーを列挙するビューを返す。生成O(1)、全列挙O(n)。
template <class K, class V> auto keys(const DefaultDict<K, V> &d) {
    return d.keys();
}

// 値を列挙するビューを返す。生成O(1)、全列挙O(n)。
template <class K, class V> auto values(const DefaultDict<K, V> &d) {
    return d.values();
}

// 登録済みの組を文字列化する。要素数nとしてO(n)と各組の出力に要する時間。
template <class K, class V> auto to_string(const DefaultDict<K, V> &d) {
    std::ostringstream out;
    out << d;
    return out.str();
}
}

namespace std {
template <class K, class V> struct hash<cplib::DefaultDict<K, V>> {
    // 登録済みの組のハッシュを返す。要素数nとしてO(n)回のハッシュ計算。
    size_t operator()(const cplib::DefaultDict<K, V> &d) const {
        return d.hash();
    }
};
}
