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
    // 空のヒープを構築する。O(1)。
    Deletable_HeapQueue() = default;

    // 配列の要素から最小ヒープを構築する。長さをnとしてO(n)。
    explicit Deletable_HeapQueue(std::span<const T> v)
        : hq(Greater{}, std::vector<T>(v.begin(), v.end())) {
    }

    // i=0のとき最小要素を返す。O(1)。
    T operator[](Int i) const {
        assert(i == 0 && !hq.empty());
        return hq.top();
    }

    // 存在するxを一つ削除する。内部ヒープの大きさをHとして償却O(log(H+1))。
    void erase(T x) {
        dlhq.push(x);
        clean();
    }

    // 存在する要素xを一つ削除する。内部ヒープの大きさをHとして償却O(log(H+1))。
    void delete_value(T x) {
        erase(x);
    }

    // 要素を追加する。内部ヒープの大きさをHとして償却O(log(H+1))。
    void push(T x) {
        hq.push(x);
    }

    // 最小要素を取り除いて返す。内部ヒープの大きさをHとして償却O(log(H+1))。
    T pop() {
        assert(!hq.empty());
        T out = hq.top();
        hq.pop();
        clean();
        return out;
    }

    // 要素数を返す。O(1)。
    Int len() const {
        return Int(hq.size()) - Int(dlhq.size());
    }

    // 要素数を返す。O(1)。
    Int size() const {
        return len();
    }
};

// 空のヒープを生成する。O(1)。
template <class T> auto initDeletableHeapQueue() {
    return Deletable_HeapQueue<T>();
}

// 配列の要素からヒープを生成する。長さをnとしてO(n)。
template <class T> auto toDeletableHeapQueue(std::span<const T> v) {
    return Deletable_HeapQueue<T>(v);
}

// 配列の要素からヒープを生成する。長さをnとしてO(n)。
template <class T> auto toDeletableHeapQueue(const std::vector<T> &v) {
    return Deletable_HeapQueue<T>(v);
}

// 要素数を返す。O(1)。
template <class T> Int len(const Deletable_HeapQueue<T> &q) {
    return q.len();
}

// 存在する要素xを一つ削除する。内部ヒープの大きさをHとして償却O(log(H+1))。
template <class T> void delete_value(Deletable_HeapQueue<T> &q, T x) {
    q.erase(x);
}

// 要素を追加する。内部ヒープの大きさをHとして償却O(log(H+1))。
template <class T> void push(Deletable_HeapQueue<T> &q, T x) {
    q.push(x);
}

// 最小要素を取り除いて返す。内部ヒープの大きさをHとして償却O(log(H+1))。
template <class T> T pop(Deletable_HeapQueue<T> &q) {
    return q.pop();
}
}
