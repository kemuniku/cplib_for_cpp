#pragma once
#include <cplib/common.hpp>
#include <queue>

namespace cplib {
template <class T> class Deletable_HeapQueue {
    struct Greater {
        bool operator()(const T &a, const T &b) const {
            return b < a;
        }
    };

    std::priority_queue<T, std::vector<T>, Greater> hq, dlhq;

    void clean() {
        while (!dlhq.empty() && !hq.empty() && dlhq.top() == hq.top()) {
            dlhq.pop();
            hq.pop();
        }
    }

public:
    Deletable_HeapQueue() = default;

    explicit Deletable_HeapQueue(std::span<const T> v)
        : hq(Greater{}, std::vector<T>(v.begin(), v.end())) {
    }

    T operator[](Int i) const {
        assert(i == 0 && !hq.empty());
        return hq.top();
    }

    // 削除予定用の第2ヒープを用いる。存在する要素のみ削除可。各操作償却O(log N)。
    void erase(T x) {
        dlhq.push(x);
        clean();
    }

    void delete_value(T x) {
        erase(x);
    }

    void push(T x) {
        hq.push(x);
    }

    T pop() {
        assert(!hq.empty());
        T out = hq.top();
        hq.pop();
        clean();
        return out;
    }

    Int len() const {
        return Int(hq.size()) - Int(dlhq.size());
    }
};

template <class T> auto initDeletableHeapQueue() {
    return Deletable_HeapQueue<T>();
}

template <class T> auto toDeletableHeapQueue(std::span<const T> v) {
    return Deletable_HeapQueue<T>(v);
}

template <class T> auto toDeletableHeapQueue(const std::vector<T> &v) {
    return Deletable_HeapQueue<T>(v);
}

template <class T> Int len(const Deletable_HeapQueue<T> &q) {
    return q.len();
}

template <class T> void delete_value(Deletable_HeapQueue<T> &q, T x) {
    q.erase(x);
}

template <class T> void push(Deletable_HeapQueue<T> &q, T x) {
    q.push(x);
}

template <class T> T pop(Deletable_HeapQueue<T> &q) {
    return q.pop();
}
}
