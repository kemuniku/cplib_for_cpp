#pragma once
#include <cplib/graph/graph.hpp>
#include <cplib/collections/lazy_leftist_heap.hpp>
#include <optional>
#include <stdexcept>

namespace cplib {
struct DirectedMSTResult {
    Int cost = 0;
    std::vector<Int> inEdge;
};
// 根から全頂点へ届く最小有向全域木の費用と各頂点への元辺IDを返す。不可能ならnullopt。
// 時間 O((V+E) log(V+E))、領域 O(V+E)。費用の桁あふれはoverflow_error。
template <class T, bool Static>
    requires std::is_signed_v<T>
std::optional<DirectedMSTResult> directedMST(const BasicGraph<T, Static, true> &g, Int root) {
    Int n = g.size();
    if (n <= 0 || root < 0 || root >= n)
        throw std::invalid_argument("非空グラフと有効な根が必要です");
    LazyLeftistHeapPool<__int128, Int> pool(g.edge_count());
    std::vector<Int> heap(2 * n, -1), parent(2 * n, -1), contractionParent(2 * n, -1);
    std::vector<Int> chosen(2 * n, -1), seen(2 * n);
    auto find = [&](Int v) {
        Int result = v;
        while (parent[result] >= 0)
            result = parent[result];
        while (parent[v] >= 0) {
            Int next = parent[v];
            parent[v] = result;
            v = next;
        }
        return result;
    };
    for (Int id = 0; id < g.edge_count(); ++id) {
        const auto &e = g.edge_info[id];
        if (e.src < 0 || e.src >= n || e.dst < 0 || e.dst >= n)
            throw std::invalid_argument("辺の端点が範囲外です");
        if (e.dst != root && e.src != e.dst)
            heap[e.dst] = pool.meld(heap[e.dst], pool.singleton(e.cost, id));
    }
    seen[root] = -1;
    Int count = n;
    for (Int start = 0; start < n; ++start) {
        Int v = find(start), stamp = start + 1;
        while (seen[v] == 0 || seen[v] == stamp) {
            if (seen[v] == stamp) {
                std::vector<Int> cycle;
                Int u = v;
                do {
                    cycle.push_back(u);
                    u = find(g.edge_info[chosen[u]].src);
                } while (u != v);
                Int contracted = count++;
                for (Int x : cycle) {
                    heap[contracted] = pool.meld(heap[contracted], heap[x]);
                    parent[x] = contractionParent[x] = contracted;
                }
                v = contracted;
            }
            seen[v] = stamp;
            while (heap[v] >= 0 && find(g.edge_info[pool.top(heap[v]).second].src) == v)
                heap[v] = pool.pop(heap[v]);
            if (heap[v] < 0)
                return std::nullopt;
            auto [weight, id] = pool.top(heap[v]);
            chosen[v] = id;
            heap[v] = pool.pop(heap[v]);
            pool.addAll(heap[v], -weight);
            v = find(g.edge_info[id].src);
        }
    }
    DirectedMSTResult result{0, std::vector<Int>(n, -1)};
    std::vector<bool> expanded(count);
    __int128 total = 0;
    for (Int v = count; v-- > 0;) {
        if (v == root || expanded[v])
            continue;
        Int id = chosen[v];
        const auto &e = g.edge_info[id];
        result.inEdge[e.dst] = id;
        total += e.cost;
        Int u = e.dst;
        while (u >= 0 && !expanded[u]) {
            expanded[u] = true;
            u = contractionParent[u];
        }
    }
    if (total < std::numeric_limits<Int>::min() || total > std::numeric_limits<Int>::max())
        throw std::overflow_error("最小有向全域木の費用がIntの範囲外です");
    result.cost = Int(total);
    return result;
}
}
