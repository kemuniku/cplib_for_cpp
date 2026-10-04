#pragma once
#include <cplib/utils/private/auto_rollback.hpp>
#include <cmath>
#include <tuple>

namespace cplib {
class RollbackMo {
    Int n, width;

    struct Query {
        Int left, right, idx;
    };

    std::vector<Query> queries;

public:
    // 長さNの列に対するRollback Moを初期化する。O(1)。
    // wが0ならQ>0でN/sqrt(Q)、Q=0でsqrt(N)を目安に幅を決める。Qは想定クエリ数。
    explicit RollbackMo(Int N, Int Q = 0, Int w = 0) : n(N) {
        assert(N >= 0 && Q >= 0 && w >= 0);
        width = std::min(std::max<Int>(1, N),
                         w > 0 ? w
                               : std::max<Int>(1, Int(Q > 0 ? double(N) / std::sqrt(double(Q))
                                                            : std::sqrt(double(N)))));
    }

    // 半開区間[l, r)を登録し、0始まりの登録番号を返す。空区間も許可する。償却O(1)。
    Int insert(Int l, Int r) {
        assert(0 <= l && l <= r && r <= n);
        Int id = queries.size();
        queries.push_back({l, r, id});
        return id;
    }

    template <class Add, class Rollback, class Answer>
    void run(Add add, Rollback rollback, Answer answer) const {
        auto qs = queries;
        std::sort(qs.begin(), qs.end(), [&](auto a, auto b) {
            return std::tuple(a.left / width, a.right, a.idx) <
                   std::tuple(b.left / width, b.right, b.idx);
        });
        Int block = -1, boundary = 0, right = 0;
        for (auto q : qs) {
            Int nb = q.left / width;
            if (nb != block) {
                while (right > boundary) {
                    rollback();
                    --right;
                }
                block = nb;
                boundary = q.left - q.left % width;
                boundary += std::min(width, n - boundary);
                right = boundary;
            }
            if (q.right <= boundary) {
                for (Int i = q.left; i < q.right; ++i)
                    add(i);
                answer(q.idx);
                for (Int i = q.left; i < q.right; ++i)
                    rollback();
            } else {
                while (right < q.right)
                    add(right++);
                for (Int l = boundary; l > q.left;)
                    add(--l);
                answer(q.idx);
                for (Int i = q.left; i < boundary; ++i)
                    rollback();
            }
        }
        while (right > boundary) {
            rollback();
            --right;
        }
    }

    template <class Add, class Answer> void runAutoRollback(Add add, Answer answer) const {
        detail::run_auto_rollback(*this, add, answer);
    }
};

inline auto initRollbackMo(Int n, Int q = 0, Int width = 0) {
    return RollbackMo(n, q, width);
}
}
