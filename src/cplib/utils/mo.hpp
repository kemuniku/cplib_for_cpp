#pragma once
#include <cplib/common.hpp>
#include <cmath>
#include <cplib/graph/graph.hpp>
#include <tuple>

namespace cplib {
class Mo {
    Int N, size = 0;
    std::vector<std::vector<Int>> qli;

public:
    Int width;

    // 想定長・クエリ数から幅を選び初期化する。O(N/width)。
    Mo(Int n, Int q, Int w = -1)
        : N(n),
          width(w < 0
                    ? std::max<Int>(1, Int(double(n) / std::max(1.0, std::sqrt(double(q) * 2 / 3))))
                    : w) {
        assert(width > 0);
        qli.resize(n / width + 1);
    }

    // 座標を圧縮して登録する。償却O(1)、負のrやl>rも許す。
    // lとクエリ番号は非負20bit、rは符号付き24bitで格納する。
    void insert(Int l, Int r) {
        assert(0 <= l && l <= N && l < (1 << 20));
        assert(-(1 << 23) <= r && r < (1 << 23) && r <= N);
        assert(size < (1 << 20));
        UInt packed = (UInt(r) << 40) | (UInt(l) << 20) | UInt(size++);
        qli[l / width].push_back(std::bit_cast<Int>(packed));
    }

    // ソートO(Q log Q)、端点移動O(N²/width+Q*width)。
    template <class AL, class AR, class DL, class DR, class Remember>
    void run(AL add_left, AR add_right, DL delete_left, DR delete_right, Remember remember) {
        Int nl = 0, nr = 0;
        constexpr Int mask = (1 << 20) - 1;
        for (std::size_t i = 0; i < qli.size(); ++i) {
            auto &bucket = qli[i];
            std::sort(bucket.begin(), bucket.end());
            if (i & 1)
                std::reverse(bucket.begin(), bucket.end());
            for (Int x : bucket) {
                Int ri = x >> 40, li = (x >> 20) & mask, idx = x & mask;
                while (nl > li)
                    add_left(--nl);
                while (nr < ri)
                    add_right(nr++);
                while (nl < li)
                    delete_left(nl++);
                while (nr > ri)
                    delete_right(--nr);
                remember(idx);
            }
        }
    }
};

inline Mo initMo(Int n, Int q, Int width = -1) {
    return {n, q, width};
}

inline void insert(Mo &m, Int l, Int r) {
    m.insert(l, r);
}

template <class... Callbacks> void run(Mo &m, Callbacks &&...callbacks) {
    m.run(std::forward<Callbacks>(callbacks)...);
}
}

namespace cplib {
class TreeMo {
    std::vector<Int> tour, first;

    struct Query {
        Int left, right, idx;
    };

    std::vector<Query> queries;

public:
    Int width;

    template <class G> TreeMo(const G &g, Int Q, Int root = 0, Int w = 0) {
        Int n;
        if constexpr (GraphTypes<G>)
            n = g.len;
        else
            n = g.size();
        assert(n > 0 && 0 <= root && root < n && Q >= 0 && w >= 0);
        first.assign(n, -1);
        std::vector<std::pair<Int, Int>> stack{{root, -1}};
        while (!stack.empty()) {
            Int u = stack.back().first, p = stack.back().second;
            stack.pop_back();
            if (u == -1) {
                tour.push_back(p);
                continue;
            }
            assert(0 <= u && u < n && first[u] == -1);
            first[u] = tour.size();
            tour.push_back(u);
            auto next = [&](Int v) {
                if (v != p) {
                    stack.emplace_back(-1, u);
                    stack.emplace_back(v, u);
                }
            };
            if constexpr (GraphTypes<G>) {
                for (auto [v, c] : g.to_and_cost(u)) {
                    (void)c;
                    next(v);
                }
            } else
                for (Int v : g[u])
                    next(v);
        }
        assert(Int(tour.size()) == 2 * n - 1);
        width = w ? std::min<Int>(w, tour.size())
                  : std::max<Int>(
                        1, Int(double(tour.size()) / std::max(1.0, std::sqrt(double(Q) * 2 / 3))));
    }

    // 両端を含む無向パス(u, v)を登録し、0始まりの登録番号を返す。償却O(1)。
    Int insert(Int u, Int v) {
        assert(0 <= u && u < Int(first.size()) && 0 <= v && v < Int(first.size()));
        Int a = first[u], b = first[v], id = queries.size();
        queries.push_back({std::min(a, b), std::max(a, b), id});
        return id;
    }

    // 両端を含む一本のパスを維持し、登録番号をrememberに渡す。終了時はrootの一頂点に戻す。
    // 呼び出し前の状態はrootだけのパスとする。初期頂点を追加するコールバックは呼ばない。
    // add(u,v)は端点uに辺(u,v)と頂点vを追加し、del(u,v)は辺(u,v)と端点vを削除する。
    // 回答はパスの向きや追加順によらないこと。rememberは状態を変更せず、登録番号で結果を保存する。
    // 実行中に登録内容を変更しないこと。幅BでソートO(Q log Q)、追加・削除O(N²/B + QB + N)回。
    template <class Add, class Delete, class Remember>
    void run(Add add, Delete del, Remember remember) {
        if (queries.empty())
            return;
        std::sort(queries.begin(), queries.end(), [&](auto a, auto b) {
            Int ba = a.left / width, bb = b.left / width;
            return std::tuple(ba, (ba & 1) ? -a.right : a.right, a.idx) <
                   std::tuple(bb, (bb & 1) ? -b.right : b.right, b.idx);
        });
        std::vector<bool> contains(first.size());
        contains[tour[0]] = true;
        Int left = 0, right = 0;
        auto move = [&](Int &p, Int target) {
            while (p != target) {
                Int step = p < target ? 1 : -1, u = tour[p], v = tour[p + step];
                if (contains[v]) {
                    del(v, u);
                    contains[u] = false;
                } else {
                    add(u, v);
                    contains[v] = true;
                }
                p += step;
            }
        };
        for (auto q : queries) {
            move(left, q.left);
            move(right, q.right);
            remember(q.idx);
        }
        move(left, 0);
        move(right, 0);
    }
};

template <class G> auto initTreeMo(const G &g, Int q, Int root = 0, Int width = 0) {
    return TreeMo(g, q, root, width);
}

inline Int insert(TreeMo &m, Int u, Int v) {
    return m.insert(u, v);
}

template <class... F> void run(TreeMo &m, F &&...f) {
    m.run(std::forward<F>(f)...);
}
}
