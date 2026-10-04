#pragma once
#include <cplib/utils/private/auto_rollback.hpp>
#include <unordered_map>
#include <variant>

namespace cplib {
template <class T = void> class OfflineDynamicQueries {
    using Value = std::conditional_t<std::is_void_v<T>, std::monostate, T>;

    struct Interval {
        Int left, right, idx;
    };

    std::unordered_map<Int, Int> active;
    std::vector<Interval> intervals;
    std::vector<Int> outputs, indices;
    std::vector<Value> values;

public:
    // 操作idxを有効にする。有効な番号の重複追加は禁止、取り消し後の再追加は可能。期待O(1)。
    void add(Int idx)
        requires std::is_void_v<T>
    {
        assert(!active.contains(idx));
        active[idx] = outputs.size();
    }

    // 操作idxと情報valueを登録する。取り消し後は別の値で再追加可能。期待・償却O(1)。
    void add(Int idx, Value value)
        requires(!std::is_void_v<T>)
    {
        assert(!active.contains(idx));
        Int i = values.size();
        active[idx] = i;
        values.push_back(std::move(value));
        indices.push_back(idx);
        intervals.push_back({Int(outputs.size()), -1, i});
    }

    // 有効な操作idxを取り消す。期待・償却O(1)。
    void remove(Int idx) {
        assert(active.contains(idx));
        Int i = active.at(idx);
        if constexpr (std::is_void_v<T>) {
            if (i < Int(outputs.size()))
                intervals.push_back({i, Int(outputs.size()), idx});
        } else
            intervals[i].right = outputs.size();
        active.erase(idx);
    }

    // 現在の状態の出力を登録し、実行時にidxをanswerへ渡す。償却O(1)。
    void output(Int idx) {
        outputs.push_back(idx);
    }

    template <class Apply, class Rollback, class Answer>
    void run(Apply apply, Rollback rollback, Answer answer) const {
        Int q = outputs.size();
        if (!q)
            return;
        Int size = 1;
        while (size < q)
            size *= 2;
        std::vector<std::vector<Int>> nodes(size * 2);
        auto insert = [&](Int l, Int r, Int idx) {
            for (l += size, r += size; l < r; l /= 2, r /= 2) {
                if (l & 1)
                    nodes[l++].push_back(idx);
                if (r & 1)
                    nodes[--r].push_back(idx);
            }
        };
        for (auto i : intervals)
            insert(i.left, i.right < 0 ? q : i.right, i.idx);
        if constexpr (std::is_void_v<T>)
            for (auto [idx, l] : active)
                insert(l, q, idx);
        auto visit = [&](auto &&self, Int node, Int l, Int r) -> void {
            if (l >= q)
                return;
            for (Int idx : nodes[node]) {
                if constexpr (std::is_void_v<T>)
                    apply(idx);
                else
                    apply(indices[idx], values[idx]);
            }
            if (r - l == 1)
                answer(outputs[l]);
            else {
                Int mid = (l + r) / 2;
                self(self, node * 2, l, mid);
                self(self, node * 2 + 1, mid, r);
            }
            for (std::size_t i = 0; i < nodes[node].size(); ++i)
                rollback();
        };
        visit(visit, 1, 0, size);
    }

    template <class Apply, class Answer> void runAutoRollback(Apply apply, Answer answer) const {
        detail::run_auto_rollback(*this, apply, answer);
    }
};

template <class T = void> auto initOfflineDynamicQueries() {
    return OfflineDynamicQueries<T>();
}
}
