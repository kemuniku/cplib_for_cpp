#pragma once
#include <cplib/utils/backwards_index.hpp>
#include <cplib/graph/planar_graph.hpp>
#include <cplib/tree/prufer.hpp>
#include <cplib/math/isprime.hpp>
#include <random>
#include <unordered_set>

namespace cplib {

inline std::mt19937_64 generator(std::random_device{}());

inline void seed(UInt value) {
    generator.seed(value);
}

inline Int draw(Int lo, Int hi) {
    assert(lo <= hi);
    return std::uniform_int_distribution<Int>(lo, hi)(generator);
}

template <class T> void shuffle(std::vector<T> &v) {
    std::shuffle(v.begin(), v.end(), generator);
}

struct PairHash {
    std::size_t operator()(const std::pair<Int, Int> &p) const {
        auto a = std::hash<Int>{}(p.first), b = std::hash<Int>{}(p.second);
        return a ^ (b + 0x9e3779b97f4a7c15ULL + (a << 6) + (a >> 2));
    }
};

// 一様な整数列。重複なしでは密度に応じて全列シャッフルと棄却法を使い分ける。
template <class L, class R>
std::vector<Int> randomseq(Int n, ClosedSlice<L, R> slice, bool unique = false) {
    Int l = slice.a, r = slice.b, length = r - l + 1;
    assert(n >= 0 && length >= 1);
    std::vector<Int> result;
    if (unique) {
        assert(n <= length);
        if (n >= length / 2) {
            result.resize(length);
            std::iota(result.begin(), result.end(), l);
            shuffle(result);
            result.resize(n);
        } else {
            std::unordered_set<Int> seen;
            while (Int(result.size()) != n) {
                Int x = draw(l, r);
                if (seen.insert(x).second)
                    result.push_back(x);
            }
        }
    } else {
        result.reserve(n);
        for (Int i = 0; i < n; ++i)
            result.push_back(draw(l, r));
    }
    return result;
}

// 長さn,総和がsumである各要素が非負整数である数列を一様ランダムに返す
inline std::vector<Int> randomseq_from_sum(Int n, Int sum) {
    assert(n >= 0 && sum >= 0);
    if (n == 0) {
        assert(sum == 0);
        return {};
    }
    auto cuts = randomseq(n - 1, closed_slice(Int(1), n + sum - 1), true);
    std::sort(cuts.begin(), cuts.end());
    std::vector<Int> result;
    Int now = 0;
    for (Int x : cuts) {
        result.push_back(x - now - 1);
        now = x;
    }
    result.push_back(n + sum - 1 - now);
    return result;
}

// 長さnの括弧列を返す。1は"("、-1は")"を表す。
inline std::vector<Int> random_parenthesis_sequence(Int n) {
    assert(n >= 0 && n % 2 == 0);
    Int N = n / 2, M = n / 2;
    std::vector<Int> result;
    result.reserve(n);
    for (Int i = 0; i < n; ++i) {
        Int tmp = draw(1, (N - M + 1) * (N + M));
        if (tmp <= (N - M) * (N + 1)) {
            result.push_back(-1);
            --N;
        } else {
            result.push_back(1);
            --M;
        }
    }
    return result;
}

inline std::string random_parenthesis_string(Int n) {
    std::string out;
    for (Int x : random_parenthesis_sequence(n))
        out += x == 1 ? '(' : ')';
    return out;
}

// 括弧列から二分木を復元
inline UnWeightedUnDirectedGraph make_binary_tree_from_sequence(const std::vector<Int> &ps) {
    Int n = ps.size() / 2;
    assert(n >= 1);
    std::vector<Int> stack, memo(2 * n, -1);
    for (Int i = 0; i < 2 * n; ++i) {
        if (ps[i] == 1)
            stack.push_back(i);
        else {
            assert(!stack.empty());
            memo[stack.back()] = i;
            stack.pop_back();
        }
    }
    auto g = initUnWeightedUnDirectedGraph(n);
    Int now = 0;
    auto dfs = [&](auto &&self, Int l, Int r, Int no) -> void {
        Int to = memo[l];
        if (to != l + 1) {
            ++now;
            g.add_edge(no, now);
            self(self, l + 1, to, now);
        }
        if (to + 1 != r) {
            ++now;
            g.add_edge(no, now);
            self(self, to + 1, r, now);
        }
    };
    dfs(dfs, 0, 2 * n, 0);
    assert(g.edge_count() == n - 1);
    return g;
}

// n頂点の二分木を一様ランダムに返す
inline auto random_binary_tree(Int n) {
    assert(n >= 1);
    return make_binary_tree_from_sequence(random_parenthesis_sequence(2 * n));
}

// n頂点の木を一様ランダムに返す
inline UnWeightedUnDirectedGraph random_tree(Int n) {
    assert(n >= 1);
    if (n == 1)
        return initUnWeightedUnDirectedGraph(1);
    if (n == 2) {
        auto g = initUnWeightedUnDirectedGraph(2);
        g.add_edge(0, 1);
        return g;
    }
    return prufer_decode(randomseq(n - 2, closed_slice(Int(0), n - 1)));
}

// sliceに含まれる素数を一様ランダムに返す
template <class L, class R> Int random_prime(ClosedSlice<L, R> slice) {
    Int l = slice.a, r = slice.b;
    if (r - l + 1 <= 1500) {
        bool found = false;
        for (Int x = l; x <= r; ++x)
            if (isprime(x)) {
                found = true;
                break;
            }
        assert(found);
    }
    for (;;) {
        Int x = draw(l, r);
        if (isprime(x))
            return x;
    }
}

template <class L, class R>
std::vector<Int> random_prime_sequence(Int n, ClosedSlice<L, R> slice, bool unique = false) {
    std::vector<Int> out;
    if (unique) {
        if (Int(slice.b) - Int(slice.a) + 1 <= 1000000000) {
            std::vector<Int> primes;
            for (Int x = slice.a; x <= Int(slice.b); ++x)
                if (isprime(x))
                    primes.push_back(x);
            for (Int i : randomseq(n, closed_slice(Int(0), Int(primes.size()) - 1), true))
                out.push_back(primes[i]);
        } else {
            std::unordered_set<Int> seen;
            while (Int(out.size()) != n) {
                Int x = random_prime(slice);
                if (seen.insert(x).second)
                    out.push_back(x);
            }
        }
    } else
        for (Int i = 0; i < n; ++i)
            out.push_back(random_prime(slice));
    return out;
}

// ランダムな単純グラフを作成。
inline UnWeightedUnDirectedGraph random_simple_graph(Int n, Int m) {
    assert(n >= 0 && m >= 0 && m <= n * (n - 1) / 2);
    auto g = initUnWeightedUnDirectedGraph(n);
    if (n * (n - 1) <= 10000000) {
        std::vector<std::pair<Int, Int>> candidates;
        for (Int i = 0; i < n - 1; ++i)
            for (Int j = i + 1; j < n; ++j)
                candidates.emplace_back(i, j);
        shuffle(candidates);
        for (Int i = 0; i < m; ++i)
            g.add_edge(candidates[i].first, candidates[i].second);
    } else {
        std::unordered_set<std::pair<Int, Int>, PairHash> seen;
        for (Int i = 0; i < m; ++i)
            for (;;) {
                Int u = draw(0, n - 2), v = draw(u + 1, n - 1);
                if (seen.emplace(u, v).second) {
                    g.add_edge(u, v);
                    break;
                }
            }
    }
    return g;
}

// 全頂点対を試す元の平面グラフ生成。期待 O(N³ log N)、一様分布ではない。
inline UnWeightedUnDirectedGraph random_planar_graph(Int n, Int m) {
    assert(n >= 0);
    Int maximum = n < 3 ? n * (n - 1) / 2 : 3 * n - 6;
    assert(0 <= m && m <= maximum);
    auto g = initUnWeightedUnDirectedGraph(n);
    if (!m)
        return g;
    std::vector<std::pair<Int, Int>> candidates;
    for (Int u = 0; u < n; ++u)
        for (Int v = u + 1; v < n; ++v)
            candidates.emplace_back(u, v);
    shuffle(candidates);
    for (auto [u, v] : candidates) {
        g.add_edge(u, v);
        if (is_planar_graph(g)) {
            if (g.edge_count() == m)
                return g;
        } else {
            g.edge_info.pop_back();
            g.edges[u].pop_back();
            g.edges[v].pop_back();
        }
    }
    assert(g.edge_count() == m);
    return g;
}

// ランダムな単純連結グラフを生成。ただし、一様ランダムでない。
inline UnWeightedUnDirectedGraph random_connected_graph(Int n, Int m) {
    assert(n >= 1 && m >= n - 1);
    auto g = random_tree(n);
    std::unordered_set<std::pair<Int, Int>, PairHash> seen;
    if (n * (n - 1) <= 10000000) {
        for (Int i = 0; i < n - 1; ++i)
            for (Int j = i + 1; j < n; ++j)
                seen.emplace(i, j);
        for (Int i = 0; i < n; ++i)
            for (auto e : g.adjacency(i))
                seen.erase({i, e.dst});
        assert(m - n + 1 <= Int(seen.size()));
        std::vector<std::pair<Int, Int>> candidates(seen.begin(), seen.end());
        shuffle(candidates);
        for (Int i = 0; i < m - n + 1; ++i)
            g.add_edge(candidates[i].first, candidates[i].second);
    } else {
        for (Int i = 0; i < n; ++i)
            for (auto e : g.adjacency(i))
                seen.emplace(i, e.dst);
        assert(n * (n - 1) - (n - 1) >= m);
        for (Int i = 0; i < m - n + 1; ++i)
            for (;;) {
                Int u = draw(0, n - 2), v = draw(u + 1, n - 1);
                if (seen.emplace(u, v).second) {
                    g.add_edge(u, v);
                    break;
                }
            }
    }
    return g;
}

// 1の数がoneであるような長さnの01列を一様ランダムに返す
inline std::vector<Int> random_01sequence(Int n, Int one) {
    assert(0 <= one && one <= n);
    auto indices = randomseq(one, closed_slice(Int(0), n - 1), true);
    std::vector<Int> out(n);
    for (Int i : indices)
        out[i] = 1;
    return out;
}

// sliceに含まれる文字からなる長さnの文字列を一様ランダムに返す
inline std::string random_string(Int n, ClosedSlice<char, char> slice) {
    assert(n >= 0 && slice.a <= slice.b);
    std::string out;
    out.reserve(n);
    for (Int i = 0; i < n; ++i)
        out += char(draw(static_cast<unsigned char>(slice.a), static_cast<unsigned char>(slice.b)));
    return out;
}

// alphabetの位置を一様に選んで長さnの文字列を返す。重複文字の重みも保持する。
inline std::string random_string(Int n, const std::string &alphabet) {
    assert(n >= 0 && !alphabet.empty());
    std::string out;
    out.reserve(n);
    for (Int i = 0; i < n; ++i)
        out += alphabet[draw(0, alphabet.size() - 1)];
    return out;
}
}
