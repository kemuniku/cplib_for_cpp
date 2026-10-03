#pragma once
#include <cplib/common.hpp>
#include <type_traits>

namespace cplib {
// Monotone integer heap. O(B) amortized per removed item, B = key width.
template <class K, class V> class RadixHeap {
    static_assert(std::is_integral_v<K> && !std::is_same_v<K, bool>);
    using U = std::make_unsigned_t<K>;
    std::array<std::vector<std::pair<K, V>>, sizeof(K) * 8 + 1> buckets;
    std::array<K, sizeof(K) * 8 + 1> minima{};
    K last;
    Int size_ = 0;

    static U key(K x) {
        U u = U(x);
        if constexpr (std::is_signed_v<K>)
            u ^= U(1) << (sizeof(K) * 8 - 1);
        return u;
    }

    static int bucket(K a, K b) {
        return std::bit_width(U(key(a) ^ key(b)));
    }

    void prepare() {
        assert(size_);
        if (!buckets[0].empty())
            return;
        int b = 1;
        while (buckets[b].empty())
            ++b;
        last = minima[b];
        for (auto &item : buckets[b]) {
            int d = bucket(item.first, last);
            if (buckets[d].empty() || item.first < minima[d])
                minima[d] = item.first;
            buckets[d].push_back(std::move(item));
        }
        buckets[b].clear();
    }

public:
    // minKey以上のキーを扱う単調最小ヒープを作ります。O(B)、Bはキーのビット数。
    explicit RadixHeap(K minKey = std::numeric_limits<K>::lowest()) : last(minKey) {
    }

    Int len() const {
        return size_;
    }

    bool isEmpty() const {
        return size_ == 0;
    }

    // 最後にtop/popで参照したキー以上のkを追加します。償却O(1)。同値の順序は不定です。
    void push(K k, V v) {
        assert(k >= last);
        int b = bucket(k, last);
        if (buckets[b].empty() || k < minima[b])
            minima[b] = k;
        buckets[b].emplace_back(k, std::move(v));
        ++size_;
    }

    void push(std::pair<K, V> v) {
        push(v.first, std::move(v.second));
    }

    // 最小要素を返し、以後追加できるキーの下限を更新する。push と合わせて一要素あたり償却 O(B)。
    const std::pair<K, V> &top() {
        prepare();
        return buckets[0].back();
    }

    const std::pair<K, V> &operator[](Int i) {
        assert(i == 0);
        return top();
    }

    std::pair<K, V> pop() {
        prepare();
        auto v = std::move(buckets[0].back());
        buckets[0].pop_back();
        --size_;
        return v;
    }

    // 全要素を削除してキーの下限をリセットする。O(N + B)。
    void clear(K minKey = std::numeric_limits<K>::lowest()) {
        for (auto &b : buckets)
            b.clear();
        last = minKey;
        size_ = 0;
    }
};

template <class K, class V> auto initRadixHeap(K minKey = std::numeric_limits<K>::lowest()) {
    return RadixHeap<K, V>(minKey);
}
}
