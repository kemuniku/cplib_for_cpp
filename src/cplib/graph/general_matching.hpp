#pragma once
#include <cplib/graph/graph.hpp>
#include <bit>

namespace cplib::detail {
// 葉追加可能な木の集合併合。小木の祖先マスクと境界集合の直接更新により、
// 通常のUnion-Findの逆Ackermann因子を加えず全操作O(N+操作数)に保つ。
struct MatchingTreeUnion {
    struct MicroTree {
        std::vector<Int> vertices;
        Int outside = 0;
        UInt live = 0;
    };

    Int limit;
    std::vector<Int> parent, depth, blockId, position, representative, size, first, last, next, top;
    std::vector<UInt> ancestors;
    std::vector<bool> inserted, deleted;
    std::vector<MicroTree> blocks;

    explicit MatchingTreeUnion(Int n)
        : limit(std::max<Int>(8, std::bit_width(UInt(std::max<Int>(n, 2))))), parent(n), depth(n),
          blockId(n), position(n), representative(n), size(n), first(n), last(n), next(n), top(n),
          ancestors(n), inserted(n), deleted(n), blocks{{{0}, 0, 1}} {
        inserted[0] = true;
        ancestors[0] = 1;
    }

    Int microRoot(Int v) const {
        Int b = blockId[v];
        UInt mask = ancestors[v] & blocks[b].live;
        return mask ? blocks[b].vertices[std::bit_width(mask) - 1] : blocks[b].outside;
    }

    void makeMacro(Int v) {
        if (size[v])
            return;
        representative[v] = v;
        size[v] = 1;
        first[v] = last[v] = top[v] = v;
        next[v] = -1;
    }

    void mergeMacro(Int u, Int v) {
        Int a = representative[u], b = representative[v];
        if (a == b)
            return;
        if (size[a] < size[b])
            std::swap(a, b);
        for (Int x = first[b]; x != -1; x = next[x])
            representative[x] = a;
        next[last[a]] = first[b];
        last[a] = last[b];
        size[a] += size[b];
        if (depth[top[b]] < depth[top[a]])
            top[a] = top[b];
    }

    Int root(Int v) {
        if (!inserted[v])
            return v;
        Int x = v, result = microRoot(x);
        if (blockId[x] == blockId[result])
            return result;
        x = top[representative[result]];
        for (;;) {
            result = microRoot(x);
            if (blockId[x] == blockId[result])
                return result;
            mergeMacro(result, x);
            x = top[representative[result]];
        }
    }

    void rebuildBlock(Int b, std::vector<Int> vertices, Int outside) {
        blocks[b] = {std::move(vertices), outside, 0};
        const auto &values = blocks[b].vertices;
        for (Int i = 0; i < Int(values.size()); ++i) {
            Int v = values[i];
            blockId[v] = b;
            position[v] = i;
        }
        for (Int i = 0; i < Int(values.size()); ++i) {
            Int v = values[i];
            UInt bit = UInt(1) << i;
            ancestors[v] = bit;
            if (v && blockId[parent[v]] == b)
                ancestors[v] |= ancestors[parent[v]];
            if (!deleted[v])
                blocks[b].live |= bit;
        }
    }

    // 満杯の小木を双方b/4以上に分割する。葉追加は償却O(1)。
    void splitBlock(Int b) {
        auto vertices = blocks[b].vertices;
        Int count = vertices.size(), outside = blocks[b].outside;
        std::vector<Int> sizes(count);
        for (Int i = count - 1; i >= 0; --i) {
            ++sizes[i];
            Int v = vertices[i];
            if (v && blockId[parent[v]] == b)
                sizes[position[parent[v]]] += sizes[i];
        }
        Int pivot = outside;
        for (Int i = 0; i < count; ++i)
            if (sizes[i] * 2 >= count)
                pivot = vertices[i];
        UInt selected = 0;
        Int total = 0;
        for (Int i = 0; i < count; ++i) {
            Int v = vertices[i];
            if (v && parent[v] == pivot && total * 4 < count) {
                selected |= UInt(1) << i;
                total += sizes[i];
            }
        }
        std::vector<Int> left, right;
        for (Int v : vertices)
            (ancestors[v] & selected ? right : left).push_back(v);
        assert(!left.empty() && !right.empty());
        Int nb = blocks.size();
        blocks.emplace_back();
        for (Int v : right)
            blockId[v] = nb;
        rebuildBlock(b, std::move(left), outside);
        rebuildBlock(nb, std::move(right), pivot);
        makeMacro(pivot);
    }

    void grow(Int p, Int v) {
        assert(!inserted[v] && inserted[p]);
        inserted[v] = true;
        parent[v] = p;
        depth[v] = depth[p] + 1;
        Int b = blockId[p], i = blocks[b].vertices.size();
        blockId[v] = b;
        position[v] = i;
        ancestors[v] = ancestors[p] | (UInt(1) << i);
        blocks[b].vertices.push_back(v);
        blocks[b].live |= UInt(1) << i;
        if (i + 1 == limit)
            splitBlock(b);
    }

    void joinParent(Int v) {
        assert(v && inserted[v]);
        deleted[v] = true;
        blocks[blockId[v]].live &= ~(UInt(1) << position[v]);
    }
};

struct MatchingSearch {
    struct Edge {
        Int u = 0, v = 0;
    };

    struct Contraction {
        Int v, base;
    };

    Int n, time = 0, deadline = 0, stamp = 0, order = 0, saved = 0, head = 0;
    std::vector<std::vector<Int>> graph, members;
    std::vector<Int> mate, potential, label, initial, queue;
    std::vector<Edge> link;
    std::vector<std::vector<Edge>> events;
    std::vector<Contraction> contractions;
    MatchingTreeUnion tree;

    explicit MatchingSearch(Int n)
        : n(n), graph(n + 1), members(n + 1), mate(n + 1), potential(n + 1), label(n + 1),
          initial(n + 1), link(n + 1), events(n / 2 + 1), tree(n + 1) {
    }

    Int base(Int v) {
        return tree.root(initial[v]);
    }

    void extend(Int x, Int y) {
        Int z = mate[y];
        tree.grow(x, y);
        tree.grow(y, z);
        label[y] = -1;
        potential[y] = time;
        link[z] = {x, y};
        label[z] = label[x];
        potential[z] = time + 1;
        queue.push_back(z);
    }

    void contract(Int x, Int y) {
        Int a = base(x), b = base(y);
        --stamp;
        label[mate[a]] = label[mate[b]] = stamp;
        for (;;) {
            if (mate[b])
                std::swap(a, b);
            a = base(link[a].u);
            if (label[mate[a]] == stamp)
                break;
            label[mate[a]] = stamp;
        }
        Int ancestor = a;
        for (Int endpoint : {x, y}) {
            Int v = base(endpoint);
            while (v != ancestor) {
                Int w = mate[v], parent = base(link[v].u);
                link[w] = {x, y};
                label[w] = label[x];
                potential[w] = 1 + 2 * time - potential[w];
                queue.push_back(w);
                tree.joinParent(v);
                tree.joinParent(w);
                contractions.push_back({v, ancestor});
                contractions.push_back({w, ancestor});
                v = parent;
            }
        }
    }

    // 時刻バケットで双対変数を更新し、最短増加路の長さをO(V+E)で求める。
    bool shortest() {
        time = 0;
        deadline = n + 1;
        stamp = -1;
        saved = 0;
        queue.clear();
        head = 0;
        contractions.clear();
        tree = MatchingTreeUnion(n + 1);
        for (Int u = 0; u <= n; ++u) {
            initial[u] = u;
            label[u] = 0;
            potential[u] = 1;
        }
        for (auto &bucket : events)
            bucket.clear();
        for (Int u = 1; u <= n; ++u)
            if (!mate[u]) {
                tree.grow(0, u);
                label[u] = u;
                queue.push_back(u);
            }
        for (;;) {
            while (head < Int(queue.size())) {
                Int x = queue[head++];
                for (Int y : graph[x]) {
                    if (label[y] > 0) {
                        Int at = (potential[x] + potential[y]) / 2;
                        if (label[x] != label[y]) {
                            if (at == time)
                                goto found;
                            deadline = std::min(deadline, at);
                        } else if (base(x) != base(y)) {
                            if (at == time)
                                contract(x, y);
                            else if (at <= n / 2)
                                events[at].push_back({x, y});
                        }
                    } else if (!label[y]) {
                        Int at = potential[x] + 1;
                        if (at == time)
                            extend(x, y);
                        else if (at <= n / 2)
                            events[at].push_back({x, y});
                    }
                }
            }
            for (;;) {
                ++time;
                saved = contractions.size();
                if (time > n / 2)
                    return false;
                if (time == deadline)
                    goto found;
                bool changed = false;
                for (auto [x, y] : events[time]) {
                    if (label[y] > 0) {
                        if (potential[x] + potential[y] != 2 * time || base(x) == base(y))
                            continue;
                        if (label[x] != label[y])
                            goto found;
                        contract(x, y);
                        changed = true;
                    } else if (!label[y]) {
                        extend(x, y);
                        changed = true;
                    }
                }
                if (changed)
                    break;
            }
        }
    found:
        for (Int u = 1; u <= n; ++u) {
            if (label[u] > 0)
                potential[u] -= time;
            else if (label[u] < 0)
                potential[u] = 1 + time - potential[u];
        }
        return true;
    }

    void rematch(Int v, Int w) {
        std::vector<Edge> stack{{v, w}};
        while (!stack.empty()) {
            auto [x, y] = stack.back();
            stack.pop_back();
            Int old = mate[x];
            mate[x] = y;
            if (mate[old] != x)
                continue;
            auto e = link[x];
            if (e.v == base(e.v)) {
                mate[old] = e.u;
                stack.push_back({e.u, old});
            } else {
                stack.push_back({e.v, e.u});
                stack.push_back({e.u, e.v});
            }
        }
    }

    bool augment(Int start, Int startBase) {
        struct Frame {
            Int x, bx, edge = 0;
            std::vector<Int> pending;
            Int group = 0, member = 0;
        };

        std::vector<Frame> stack{{start, startBase, 0, {}}};
        while (!stack.empty()) {
            std::size_t i = stack.size() - 1;
            if (stack[i].group < Int(stack[i].pending.size())) {
                Int b = stack[i].pending[stack[i].group];
                if (stack[i].member < Int(members[b].size())) {
                    Int v = members[b][stack[i].member++], bx = base(b);
                    stack.push_back({v, bx, 0, {}});
                } else {
                    ++stack[i].group;
                    stack[i].member = 0;
                }
                continue;
            }
            Int x = stack[i].x, bx = stack[i].bx;
            if (stack[i].edge == Int(graph[x].size())) {
                stack.pop_back();
                continue;
            }
            Int y = graph[x][stack[i].edge++];
            if (potential[x] + potential[y] != 0)
                continue;
            Int by = base(y);
            if (label[by] > 0) {
                if (label[bx] >= label[by])
                    continue;
                std::vector<Int> pending;
                Int v = by;
                while (v != bx) {
                    Int w = base(mate[v]), parent = base(link[v].u);
                    link[w] = {x, y};
                    pending.push_back(w);
                    tree.joinParent(v);
                    tree.joinParent(w);
                    v = parent;
                }
                stack[i].pending.assign(pending.rbegin(), pending.rend());
                stack[i].group = stack[i].member = 0;
            } else if (!label[by]) {
                label[by] = -1;
                Int z = mate[by];
                if (!z) {
                    rematch(x, y);
                    rematch(y, x);
                    return true;
                }
                Int bz = base(z);
                tree.grow(initial[x], by);
                tree.grow(by, bz);
                link[bz] = {x, y};
                label[bz] = order++;
                stack[i].pending = {bz};
                stack[i].group = stack[i].member = 0;
            }
        }
        return false;
    }

    // 正の双対値を持つ花を残し、頂点素な最短増加路の極大集合をO(V+E)で反転。
    Int maximal() {
        for (Int u = 0; u <= n; ++u)
            initial[u] = u;
        for (Int i = 0; i < saved; ++i) {
            auto e = contractions[i];
            initial[e.v] = e.base;
        }
        std::vector<Int> path;
        for (Int u = 1; u <= n; ++u) {
            Int v = u;
            while (initial[v] != v) {
                path.push_back(v);
                v = initial[v];
            }
            for (Int x : path)
                initial[x] = v;
            path.clear();
        }
        tree = MatchingTreeUnion(n + 1);
        order = 1;
        for (Int u = 0; u <= n; ++u) {
            label[u] = 0;
            members[u].clear();
        }
        for (Int u = 1; u <= n; ++u)
            members[initial[u]].push_back(u);
        Int result = 0;
        for (Int u = 1; u <= n; ++u) {
            if (mate[u])
                continue;
            Int b = initial[u];
            if (label[b])
                continue;
            tree.grow(0, b);
            label[b] = order++;
            for (Int v : members[b])
                if (augment(v, b)) {
                    ++result;
                    break;
                }
        }
        assert(result > 0);
        return result;
    }
};
}

namespace cplib {
// Gabow法。多重辺を許容し、自己ループと重みを無視する。
// 孤立点を除く処理を含めてO(V+E sqrt V)時間、O(V+E)領域。
// 入力は変更しない。静的グラフはbuild()後に呼ぶ。孤立点を除いた探索部分はO(E sqrt(V))。
template <UnDirectedGraph G> std::vector<std::pair<Int, Int>> maximum_matching(const G &graph) {
    if constexpr (G::is_static)
        graph.static_graph_initialized_check();
    std::vector<Int> id(graph.len), original{-1};
    for (Int u = 0; u < graph.len; ++u)
        for (auto [v, cost] : graph.to_and_cost(u))
            if (u < v)
                for (Int x : {u, v})
                    if (!id[x]) {
                        id[x] = original.size();
                        original.push_back(x);
                    }
    if (original.size() == 1)
        return {};
    Int n = original.size() - 1;
    detail::MatchingSearch search(n);
    for (Int u = 0; u < graph.len; ++u)
        for (auto [v, cost] : graph.to_and_cost(u))
            if (u < v) {
                search.graph[id[u]].push_back(id[v]);
                search.graph[id[v]].push_back(id[u]);
            }
    Int count = 0;
    while (count * 2 + 1 < n && search.shortest())
        count += search.maximal();
    std::vector<std::pair<Int, Int>> result;
    for (Int u = 1; u <= n; ++u)
        if (search.mate[u] > u)
            result.emplace_back(original[u], original[search.mate[u]]);
    return result;
}
}
