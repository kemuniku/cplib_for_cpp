#pragma once
#include <cplib/collections/avlset.hpp>

namespace cplib {
// 行・列の二つの AVL 多重集合で、最寄りの壁を O(log N) で検索。
class GridSearcher {
    using Position = std::pair<Int, Int>;
    using Maybe = std::optional<Position>;
    AvlSortedMultiSet<Position> row_, col_;

    Maybe vertical(Int i, Int j, bool down) const {
        auto out = down ? col_.gt({j, i}) : col_.lt({j, i});
        if (!out || out->first != j)
            return {};
        return Position{out->second, out->first};
    }

    Maybe horizontal(Int i, Int j, bool right) const {
        auto out = right ? row_.gt({i, j}) : row_.lt({i, j});
        if (!out || out->first != i)
            return {};
        return out;
    }

public:
    void incl(Int i, Int j) {
        row_.incl({i, j});
        col_.incl({j, i});
    }

    void incl(Position p) {
        incl(p.first, p.second);
    }

    void excl(Int i, Int j) {
        row_.excl({i, j});
        col_.excl({j, i});
    }

    void excl(Position p) {
        excl(p.first, p.second);
    }

    bool contains(Int i, Int j) const {
        return row_.contains({i, j});
    }

    bool contains(Position p) const {
        return row_.contains(p);
    }

    Int len() const {
        return row_.len();
    }

    // 同じjでiより小さい座標の最寄りの壁を返す。
    // 見つからなければstd::nullopt。
    Maybe up(Int i, Int j) const {
        return vertical(i, j, false);
    }

    // 同じjでiより大きい座標の最寄りの壁を返す。
    // 見つからなければstd::nullopt。
    Maybe down(Int i, Int j) const {
        return vertical(i, j, true);
    }

    // 同じiでjより小さい座標の最寄りの壁を返す。
    // 見つからなければstd::nullopt。
    Maybe left(Int i, Int j) const {
        return horizontal(i, j, false);
    }

    // 同じiでjより大きい座標の最寄りの壁を返す。
    // 見つからなければstd::nullopt。
    Maybe right(Int i, Int j) const {
        return horizontal(i, j, true);
    }

    // 同じjでiより小さい座標の最寄りの壁を探し、壁の一つ下のマスを返す。
    // 見つからなければstd::nullopt。
    Maybe up_move(Int i, Int j) const {
        auto p = up(i, j);
        if (p)
            ++p->first;
        return p;
    }

    // 同じjでiより大きい座標の最寄りの壁を探し、壁の一つ上のマスを返す。
    // 見つからなければstd::nullopt。
    Maybe down_move(Int i, Int j) const {
        auto p = down(i, j);
        if (p)
            --p->first;
        return p;
    }

    // 同じiでjより小さい座標の最寄りの壁を探し、壁の一つ右のマスを返す。
    // 見つからなければstd::nullopt。
    Maybe left_move(Int i, Int j) const {
        auto p = left(i, j);
        if (p)
            ++p->second;
        return p;
    }

    // 同じiでjより大きい座標の最寄りの壁を探し、壁の一つ左のマスを返す。
    // 見つからなければstd::nullopt。
    Maybe right_move(Int i, Int j) const {
        auto p = right(i, j);
        if (p)
            --p->second;
        return p;
    }

    auto updownleftright(Int i, Int j) const {
        return std::array<Maybe, 4>{up(i, j), down(i, j), left(i, j), right(i, j)};
    }

    auto updownleftright_move(Int i, Int j) const {
        return std::array<Maybe, 4>{up_move(i, j), down_move(i, j), left_move(i, j),
                                    right_move(i, j)};
    }

    auto updownleftright_get(Int i, Int j) const {
        std::vector<Position> out;
        for (auto p : updownleftright(i, j))
            if (p)
                out.push_back(*p);
        return out;
    }

    auto updownleftright_move_get(Int i, Int j) const {
        std::vector<Position> out;
        for (auto p : updownleftright_move(i, j))
            if (p)
                out.push_back(*p);
        return out;
    }

#define CPLIB_GRID_PAIR_OVERLOAD(name)                                                             \
    auto name(Position p) const {                                                                  \
        return name(p.first, p.second);                                                            \
    }
    CPLIB_GRID_PAIR_OVERLOAD(up)
    CPLIB_GRID_PAIR_OVERLOAD(down)
    CPLIB_GRID_PAIR_OVERLOAD(left)
    CPLIB_GRID_PAIR_OVERLOAD(right)
    CPLIB_GRID_PAIR_OVERLOAD(up_move)
    CPLIB_GRID_PAIR_OVERLOAD(down_move)
    CPLIB_GRID_PAIR_OVERLOAD(left_move)
    CPLIB_GRID_PAIR_OVERLOAD(right_move)
    CPLIB_GRID_PAIR_OVERLOAD(updownleftright)
    CPLIB_GRID_PAIR_OVERLOAD(updownleftright_move)
    CPLIB_GRID_PAIR_OVERLOAD(updownleftright_get)
    CPLIB_GRID_PAIR_OVERLOAD(updownleftright_move_get)
#undef CPLIB_GRID_PAIR_OVERLOAD
};

inline GridSearcher initGridSearcher() {
    return {};
}
}
