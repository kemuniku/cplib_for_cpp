#pragma once
#include <cplib/graph/graph.hpp>
#include <cplib/tree/private/lca_build.hpp>
#include <bit>

namespace cplib {
// Schieber–Vishkin法とEuler tour/ladders。構築O(N)、LCA/LA照会O(1)。
class LowestCommonAncestor {
    Int n = 0, root = 0, kind = 0;
    bool ordered = false, laEnabled = true;
    std::vector<std::uint32_t> data, prefixLCA, pathAttach, laIndex;
    std::vector<std::uint8_t> prefixBranch;
    std::vector<std::int32_t> parents, headParent, laJump, laLadder;
    std::vector<Int> laBase;

    static int log2(std::uint32_t x) {
        assert(x);
        return 31 - std::countl_zero(x);
    }

    // Euler tourのジャンプと、最長経路を5倍に延長したladderを構築する。時間・空間O(N)
    void buildLA(const std::vector<Int> &order) {
        auto vertex = [&](Int i) { return ordered ? i : order[i]; };
        auto dep = [&](Int v) -> Int { return data[3 * v + 2]; };
        std::vector<std::int32_t> leaf(n);
        std::iota(leaf.begin(), leaf.end(), 0);
        for (Int i = n - 1; i > 0; --i) {
            Int v = vertex(i), p = parents[v];
            if (dep(leaf[v]) > dep(leaf[p]))
                leaf[p] = leaf[v];
        }
        laBase.resize(n);
        Int total = 0;
        for (Int i = 0; i < n; ++i) {
            Int v = vertex(i);
            if (v == root || leaf[v] != leaf[parents[v]]) {
                Int bottom = leaf[v],
                    length = std::min(5 * (dep(bottom) - dep(v) + 1), dep(bottom) + 1);
                laBase[bottom] = total + dep(bottom);
                total += length;
            }
        }
        laLadder.resize(total);
        for (Int i = 0; i < n; ++i) {
            Int v = vertex(i);
            if (v == root || leaf[v] != leaf[parents[v]]) {
                Int bottom = leaf[v],
                    length = std::min(5 * (dep(bottom) - dep(v) + 1), dep(bottom) + 1),
                    start = laBase[bottom] - dep(bottom), u = bottom;
                for (Int j = 0; j < length; ++j) {
                    laLadder[start + j] = u;
                    u = parents[u];
                }
            }
        }
        for (Int v = 0; v < n; ++v)
            laBase[v] = laBase[leaf[v]];
        std::vector<std::int32_t> child(n, -1), sibling(n);
        for (Int v = 0; v < n; ++v)
            if (v != root) {
                Int p = parents[v];
                sibling[v] = child[p];
                child[p] = v;
            }
        laIndex.resize(n);
        laJump.resize(2 * n);
        std::vector<std::int32_t> stack;
        stack.reserve(dep(leaf[root]) + 1);
        stack.push_back(root);
        Int index = 1;
        laIndex[root] = index;
        laJump[index] = root;
        while (!stack.empty()) {
            Int v = stack.back(), u = child[v];
            if (u != -1) {
                child[v] = sibling[u];
                stack.push_back(u);
                laIndex[u] = ++index;
            } else {
                stack.pop_back();
                if (stack.empty())
                    break;
                ++index;
            }
            Int jump = index & -index;
            laJump[index] = stack[std::max<Int>(0, Int(stack.size()) - 1 - jump)];
        }
    }

    // 同じ接点を持つ頂点間のLCAを、縦パスのビット演算で返す。O(1)
    [[gnu::noinline]] Int lcaInsideBranch(Int u, Int v) const {
        std::uint32_t a = data[3 * u], b = data[3 * v];
        Int x = u, y = v;
        if (a != b) {
            auto common = data[3 * u + 1] & data[3 * v + 1] & (~0u << log2(a ^ b));
            auto lowA = data[3 * u + 1] ^ common, lowB = data[3 * v + 1] ^ common;
            if (lowA) {
                int k = log2(lowA);
                x = headParent[(a & (~0u << k)) | (1u << k)];
            }
            if (lowB) {
                int k = log2(lowB);
                y = headParent[(b & (~0u << k)) | (1u << k)];
            }
        }
        return ordered ? std::min(x, y) : (data[3 * x + 2] <= data[3 * y + 2] ? x : y);
    }

public:
    LowestCommonAncestor() = default;

    // 親配列から構築。parent[root]は不使用。no_laでLAの前計算を省略する。
    explicit LowestCommonAncestor(std::span<const Int> parent, Int root_, bool no_la = false)
        : n(parent.size()), root(root_), laEnabled(!no_la), parents(n) {
        assert(n <= std::numeric_limits<std::int32_t>::max() && 0 <= root && root < n);
        bool invalid = false, unordered = root != 0, nonPath = root != 0, nonStar = false;
        for (Int v = 0; v < n; ++v)
            if (v != root) {
                Int p = parent[v];
                invalid |= p < 0 || p >= n;
                unordered |= p >= v;
                nonPath |= p != v - 1;
                nonStar |= p != root;
                parents[v] = p;
            }
        assert(!invalid);
        ordered = !unordered;
        parents[root] = -1;
        if (!nonPath) {
            kind = 1;
            return;
        }
        if (!nonStar) {
            kind = 2;
            return;
        }
        data.resize(3 * n);
        prefixBranch.resize(n);
        prefixLCA.resize(64 * 64);
        std::vector<Int> order;
        if (!ordered) {
            std::vector<Int> head(n, -1), next(n);
            for (Int v = 0; v < n; ++v)
                if (v != root) {
                    next[v] = head[parent[v]];
                    head[parent[v]] = v;
                }
            order.reserve(n);
            order.push_back(root);
            for (std::size_t i = 0; i < order.size(); ++i)
                for (Int v = head[order[i]]; v != -1; v = next[v])
                    order.push_back(v);
            assert(Int(order.size()) == n);
        }
        auto vertex = [&](Int i) { return ordered ? i : order[i]; };
        std::vector<std::int32_t> size(n);
        headParent.resize(n + 1);
        Int deepest =
            ordered ? detail::lca_native::build<true>(parents.data(), nullptr, n, root, size.data(),
                                                      data.data(), prefixBranch.data(),
                                                      headParent.data())
                    : detail::lca_native::build<false>(parents.data(), order.data(), n, root,
                                                       size.data(), data.data(),
                                                       prefixBranch.data(), headParent.data());
        prefixLCA[0] = root;
        for (Int i = 1; i < std::min<Int>(n, 64); ++i) {
            Int v = vertex(i), p = prefixBranch[parents[v]];
            prefixLCA[i * 64 + i] = v;
            for (Int j = 0; j < i; ++j) {
                auto ancestor = prefixLCA[p * 64 + j];
                prefixLCA[i * 64 + j] = ancestor;
                prefixLCA[j * 64 + i] = ancestor;
            }
        }
        if (data[3 * deepest + 2] >= 64 && std::popcount(data[3 * deepest + 1]) > 2) {
            pathAttach.assign(n, ~0u);
            for (Int v = deepest; v != -1; v = parents[v])
                pathAttach[v] = v;
            for (Int i = 1; i < n; ++i) {
                Int v = vertex(i);
                if (pathAttach[v] == ~0u)
                    pathAttach[v] = pathAttach[parents[v]];
            }
            kind = 3;
        }
        if (!no_la)
            buildLA(order);
    }

    // 頂点数を返す。森の場合は追加した根を含む。O(1)
    Int numVertices() const {
        return parents.size();
    }

    // 頂点vの親を返す。根の場合は-1を返す。O(1)
    Int parentOf(Int v) const {
        return parents[v];
    }

    // 根から頂点vまでの辺数を返す。O(1)
    Int depth(Int v) const {
        assert(0 <= v && v < n);
        return kind == 1 ? v : kind == 2 ? Int(v != root) : data[3 * v + 2];
    }

    // vからk辺上の祖先を返す。LA有効での構築が必要。kが負または根を超える場合は-1。O(1)
    Int la(Int v, Int k) const {
        assert(laEnabled && 0 <= v && v < n);
        if (k < 0)
            return -1;
        Int targetDepth = depth(v) - k;
        if (targetDepth < 0)
            return -1;
        if (k == 0)
            return v;
        if (k == 1)
            return parents[v];
        if (kind == 1)
            return v - k;
        Int step = Int(1) << log2(k >> 1), index = (Int(laIndex[v]) & -step) | step,
            jump = laJump[index];
        return laLadder[laBase[jump] - targetDepth];
    }

    // 頂点uとvの最小共通祖先を返す。O(1)
    Int lca(Int u, Int v) const {
        assert(0 <= u && u < n && 0 <= v && v < n);
        if (kind != 0) {
            if (kind == 3) {
                Int a = pathAttach[u], b = pathAttach[v];
                if (a != b)
                    return ordered ? std::min(a, b) : (data[3 * a + 2] <= data[3 * b + 2] ? a : b);
            } else if (kind == 1)
                return std::min(u, v);
            else
                return u == v ? u : root;
        }
        Int a = prefixBranch[u], b = prefixBranch[v];
        return a != b ? prefixLCA[a * 64 + b] : lcaInsideBranch(u, v);
    }

    // 頂点uとvを結ぶパスの辺数を返す。O(1)
    Int dist(Int u, Int v) const {
        return depth(u) + depth(v) - 2 * depth(lca(u, v));
    }

    // startingからgoalへd辺進んだ頂点を返す。dが負またはパスの辺数を超える場合は-1。O(1)
    Int la(Int starting, Int goal, Int d) const {
        assert(laEnabled);
        Int ancestor = lca(starting, goal), up = depth(starting) - depth(ancestor),
            length = up + depth(goal) - depth(ancestor);
        if (d < 0 || d > length)
            return -1;
        return d <= up ? la(starting, d) : la(goal, length - d);
    }

    // 根をxとしたときのyとzの最小共通祖先を返す。O(1)
    Int median(Int x, Int y, Int z) const {
        return lca(x, y) ^ lca(y, z) ^ lca(x, z);
    }
};

// 根付き木の親配列から構築する。parent[root]は参照しない。N < 2^31。時間・空間O(N)
// no_la=trueでLA用の前計算を省略する。
inline auto initLCAFromParent(std::span<const Int> parent, Int root, bool no_la = false) {
    return LowestCommonAncestor(parent, root, no_la);
}

namespace detail {
inline auto lca_undirected_adj(const std::vector<std::vector<Int>> &adj) {
    std::vector<std::vector<Int>> result(adj.size());
    for (Int v = 0; v < Int(adj.size()); ++v)
        for (Int u : adj[v]) {
            assert(0 <= u && u < Int(adj.size()));
            result[v].push_back(u);
            result[u].push_back(v);
        }
    return result;
}

template <class G>
    requires(UnDirectedGraph<G> || DirectedGraph<G>)
auto lca_undirected_adj(const G &g) {
    std::vector<std::vector<Int>> result(g.len);
    for (Int v = 0; v < g.len; ++v)
        for (auto [u, c] : g.to_and_cost(v)) {
            result[v].push_back(u);
            result[u].push_back(v);
        }
    return result;
}

template <class Adj> auto lca_from_adj(const Adj &adj, Int root, bool forest, bool no_la) {
    Int n;
    if constexpr (requires { adj.len; })
        n = adj.len;
    else
        n = adj.size();
    std::vector<Int> parent(n + Int(forest), -2), stack;
    if (forest)
        parent[n] = -1;
    else {
        assert(0 <= root && root < n);
        parent[root] = -1;
        stack.push_back(root);
    }
    for (Int start = 0; start < std::max<Int>(1, n); ++start) {
        if (forest) {
            if (start == n || parent[start] != -2)
                continue;
            parent[start] = n;
            stack.push_back(start);
        } else if (start != 0)
            break;
        while (!stack.empty()) {
            Int v = stack.back();
            stack.pop_back();
            auto visit = [&](Int u) {
                if (parent[u] == -2) {
                    parent[u] = v;
                    stack.push_back(u);
                }
            };
            if constexpr (UnDirectedGraph<Adj>) {
                for (auto [u, c] : adj.to_and_cost(v))
                    visit(u);
            } else
                for (Int u : adj[v])
                    visit(u);
        }
    }
    return initLCAFromParent(parent, forest ? n : root, no_la);
}
}

// 辺の向きと重みを無視すると木になるgから構築する。時間・空間O(N + M)、木ではO(N)
// no_la=trueでLA用の前計算を省略する。
template <class G>
    requires(UnDirectedGraph<G> || DirectedGraph<G>)
auto initLCA(const G &g, Int root, bool no_la = false) {
    if constexpr (UnDirectedGraph<G>)
        return detail::lca_from_adj(g, root, false, no_la);
    else
        return detail::lca_from_adj(detail::lca_undirected_adj(g), root, false, no_la);
}

// 木の隣接リストから辺の向きを無視して構築する。時間・空間O(N + M)、木ではO(N)
// no_la=trueでLA用の前計算を省略する。
inline auto initLCA(const std::vector<std::vector<Int>> &adj, Int root, bool no_la = false) {
    return detail::lca_from_adj(detail::lca_undirected_adj(adj), root, false, no_la);
}

// 森に根Nを追加し、各成分の最小番号の頂点と結ぶ。辺の向きと重みは無視する。O(N + M)
// no_la=trueでLA用の前計算を省略する。
template <class G>
    requires(UnDirectedGraph<G> || DirectedGraph<G>)
auto initLCAFromForest(const G &g, bool no_la = false) {
    if constexpr (UnDirectedGraph<G>)
        return detail::lca_from_adj(g, g.len, true, no_la);
    else
        return detail::lca_from_adj(detail::lca_undirected_adj(g), g.len, true, no_la);
}

// 森の隣接リストの向きを無視し、根Nを追加して各成分の最小頂点と結ぶ。O(N + M)
// no_la=trueでLA用の前計算を省略する。
inline auto initLCAFromForest(const std::vector<std::vector<Int>> &adj, bool no_la = false) {
    return detail::lca_from_adj(detail::lca_undirected_adj(adj), adj.size(), true, no_la);
}

inline Int numVertices(const LowestCommonAncestor &t) {
    return t.numVertices();
}

inline Int parentOf(const LowestCommonAncestor &t, Int v) {
    return t.parentOf(v);
}

inline Int depth(const LowestCommonAncestor &t, Int v) {
    return t.depth(v);
}

inline Int la(const LowestCommonAncestor &t, Int v, Int k) {
    return t.la(v, k);
}

inline Int la(const LowestCommonAncestor &t, Int s, Int g, Int d) {
    return t.la(s, g, d);
}

inline Int lca(const LowestCommonAncestor &t, Int u, Int v) {
    return t.lca(u, v);
}

inline Int dist(const LowestCommonAncestor &t, Int u, Int v) {
    return t.dist(u, v);
}

inline Int median(const LowestCommonAncestor &t, Int x, Int y, Int z) {
    return t.median(x, y, z);
}
}
