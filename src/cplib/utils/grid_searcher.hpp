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
    // 座標を壁として追加する。登録数をNとして計算量 O(log(N+1))。
    void incl(Int i, Int j) {
        row_.incl({i, j});
        col_.incl({j, i});
    }

    // 座標を壁として追加する。登録数をNとして計算量 O(log(N+1))。
    void incl(Position p) {
        incl(p.first, p.second);
    }

    // 座標にある壁を一つ削除する。登録数をNとして計算量 O(log(N+1))。
    void excl(Int i, Int j) {
        row_.excl({i, j});
        col_.excl({j, i});
    }

    // 座標にある壁を一つ削除する。登録数をNとして計算量 O(log(N+1))。
    void excl(Position p) {
        excl(p.first, p.second);
    }

    // 座標が壁か判定する。登録数をNとして計算量 O(log(N+1))。
    bool contains(Int i, Int j) const {
        return row_.contains({i, j});
    }

    // 座標が壁か判定する。登録数をNとして計算量 O(log(N+1))。
    bool contains(Position p) const {
        return row_.contains(p);
    }

    // 登録された壁の個数を返す。計算量 O(1)。
    Int size() const {
        return row_.size();
    }

    // 同じjでiより小さい座標の最寄りの壁を返す。
    // 見つからなければstd::nullopt。登録数をNとして計算量 O(log(N+1))。
    Maybe up(Int i, Int j) const {
        return vertical(i, j, false);
    }

    // 同じjでiより大きい座標の最寄りの壁を返す。
    // 見つからなければstd::nullopt。登録数をNとして計算量 O(log(N+1))。
    Maybe down(Int i, Int j) const {
        return vertical(i, j, true);
    }

    // 同じiでjより小さい座標の最寄りの壁を返す。
    // 見つからなければstd::nullopt。登録数をNとして計算量 O(log(N+1))。
    Maybe left(Int i, Int j) const {
        return horizontal(i, j, false);
    }

    // 同じiでjより大きい座標の最寄りの壁を返す。
    // 見つからなければstd::nullopt。登録数をNとして計算量 O(log(N+1))。
    Maybe right(Int i, Int j) const {
        return horizontal(i, j, true);
    }

    // 同じjでiより小さい座標の最寄りの壁を探し、壁の一つ下のマスを返す。
    // 見つからなければstd::nullopt。登録数をNとして計算量 O(log(N+1))。
    Maybe up_move(Int i, Int j) const {
        auto p = up(i, j);
        if (p)
            ++p->first;
        return p;
    }

    // 同じjでiより大きい座標の最寄りの壁を探し、壁の一つ上のマスを返す。
    // 見つからなければstd::nullopt。登録数をNとして計算量 O(log(N+1))。
    Maybe down_move(Int i, Int j) const {
        auto p = down(i, j);
        if (p)
            --p->first;
        return p;
    }

    // 同じiでjより小さい座標の最寄りの壁を探し、壁の一つ右のマスを返す。
    // 見つからなければstd::nullopt。登録数をNとして計算量 O(log(N+1))。
    Maybe left_move(Int i, Int j) const {
        auto p = left(i, j);
        if (p)
            ++p->second;
        return p;
    }

    // 同じiでjより大きい座標の最寄りの壁を探し、壁の一つ左のマスを返す。
    // 見つからなければstd::nullopt。登録数をNとして計算量 O(log(N+1))。
    Maybe right_move(Int i, Int j) const {
        auto p = right(i, j);
        if (p)
            --p->second;
        return p;
    }

    // 上下左右それぞれの最寄りの壁を返す。登録数をNとして計算量 O(log(N+1))。
    auto updownleftright(Int i, Int j) const {
        return std::array<Maybe, 4>{up(i, j), down(i, j), left(i, j), right(i, j)};
    }

    // 上下左右それぞれの壁の手前のマスを返す。登録数をNとして計算量 O(log(N+1))。
    auto updownleftright_move(Int i, Int j) const {
        return std::array<Maybe, 4>{up_move(i, j), down_move(i, j), left_move(i, j),
                                    right_move(i, j)};
    }

    // 上下左右の最寄りの壁のうち存在するものを列挙する。登録数をNとして計算量 O(log(N+1))。
    auto updownleftright_get(Int i, Int j) const {
        std::vector<Position> out;
        for (auto p : updownleftright(i, j))
            if (p)
                out.push_back(*p);
        return out;
    }

    // 上下左右の壁の手前のマスのうち存在するものを列挙する。登録数をNとして計算量 O(log(N+1))。
    auto updownleftright_move_get(Int i, Int j) const {
        std::vector<Position> out;
        for (auto p : updownleftright_move(i, j))
            if (p)
                out.push_back(*p);
        return out;
    }

// 座標の組を受け取る各オーバーロードの意味・計算量は対応する二引数版と同じ。
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

// 空の壁集合を初期化する。計算量 O(1)。
inline GridSearcher initGridSearcher() {
    return {};
}
}
