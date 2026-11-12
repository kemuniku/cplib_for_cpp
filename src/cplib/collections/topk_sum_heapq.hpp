#pragma once
#include <cplib/collections/deletable_heapqueue.hpp>

namespace cplib {
class TopK_sum_heapq {
    Deletable_HeapQueue<Int> G, L;

public:
    Int sm = 0, topk = 0, k;

    // 配列から構築し、大きいk個の和を管理する。要素数をnとしてO(n log(n+1))。
    TopK_sum_heapq(std::span<const Int> input, Int k) : k(k) {
        assert(0 <= k && k <= Int(input.size()));
        std::vector<Int> v(input.begin(), input.end());
        std::sort(v.begin(), v.end(), std::greater<Int>());
        for (Int i = 0; i < k; ++i) {
            G.push(v[i]);
            sm += v[i];
            topk += v[i];
        }
        for (Int i = k; i < Int(v.size()); ++i) {
            L.push(-v[i]);
            sm += v[i];
        }
    }

    // xを追加する。内部ヒープの大きさをHとして償却O(log(H+1))。
    void incl(Int x) {
        sm += x;
        if (G.len() < k) {
            topk += x;
            G.push(x);
        } else if (k == 0)
            L.push(-x);
        else if (G[0] <= x) {
            topk += x;
            G.push(x);
            Int tmp = G.pop();
            L.push(-tmp);
            topk -= tmp;
        } else
            L.push(-x);
    }

    // 存在するxを一つ削除する。内部ヒープの大きさをHとして償却O(log(H+1))。
    void excl(Int x) {
        sm -= x;
        if (k == 0)
            L.erase(-x);
        else if (G[0] <= x) {
            topk -= x;
            G.erase(x);
            if (L.len() > 0) {
                Int tmp = -L.pop();
                G.push(tmp);
                topk += tmp;
            }
        } else
            L.erase(-x);
    }

    // kを1増やす。内部ヒープの大きさをHとして償却O(log(H+1))。
    void addK() {
        Int tmp = -L.pop();
        topk += tmp;
        ++k;
        G.push(tmp);
    }

    // kを1減らす。内部ヒープの大きさをHとして償却O(log(H+1))。
    void minusK() {
        Int tmp = G.pop();
        topk -= tmp;
        --k;
        L.push(-tmp);
    }

    // kをvalueに変更する。内部ヒープの大きさをHとして償却O(|k-value| log(H+1))。
    // 変更幅に比例する回数のヒープ操作。
    void setK(Int value) {
        while (k > value)
            minusK();
        while (k < value)
            addK();
    }
};

// 配列から上位k個の和を管理する構造を生成する。要素数をnとしてO(n log(n+1))。
inline auto initTopKHeapq(std::span<const Int> v, Int k) {
    return TopK_sum_heapq(v, k);
}

// xを追加する。内部ヒープの大きさをHとして償却O(log(H+1))。
inline void incl(TopK_sum_heapq &q, Int x) {
    q.incl(x);
}

// 存在するxを一つ削除する。内部ヒープの大きさをHとして償却O(log(H+1))。
inline void excl(TopK_sum_heapq &q, Int x) {
    q.excl(x);
}

// kを1増やす。内部ヒープの大きさをHとして償却O(log(H+1))。
inline void addK(TopK_sum_heapq &q) {
    q.addK();
}

// kを1減らす。内部ヒープの大きさをHとして償却O(log(H+1))。
inline void minusK(TopK_sum_heapq &q) {
    q.minusK();
}

// kを変更する。変更幅をd、内部ヒープの大きさをHとして償却O(d log(H+1))。
inline void setK(TopK_sum_heapq &q, Int k) {
    q.setK(k);
}
}
