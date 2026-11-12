#pragma once
#include <cplib/collections/private/range_sort_impl.hpp>

namespace cplib {
template <class T> class RangeSortSegmentTree : public detail::RangeSortEngine<T, true> {
    using Base = detail::RangeSortEngine<T, true>;

public:
    // 相異なる[0,limit)のキーと同数の値から構築する。演算O(1)なら要素数Nとして期待O(N)時間・O(N)空間。
    template <class Op>
    RangeSortSegmentTree(std::span<const Int> keys, std::span<const T> values, Int limit, Op op,
                         T e)
        : Base(keys, values, limit, op, e) {
    }

    using Base::len;
    using Base::get;
    using Base::sort;

    // 現在の位置iの値を返す。キー上限UとしてO(log(N+1)+log(U+1))。
    T operator[](Int i) const {
        return this->value(i);
    }

    // 末尾から指定した位置の値を返す。キー上限UとしてO(log(N+1)+log(U+1))。
    T operator[](BackwardsIndex i) const {
        return this->value(len() - i.value);
    }

    // 現在の位置iの値を更新する。キー上限UとしてO(log(N+1)+log(U+1))。
    void set(Int i, T v) {
        this->update(i, std::move(v));
    }

    // 末尾から指定した位置の値を更新する。キー上限UとしてO(log(N+1)+log(U+1))。
    void set(BackwardsIndex i, T v) {
        set(len() - i.value, std::move(v));
    }

    // 閉区間の値を現在の並び順で集約する。キー上限UとしてO(log(N+1)+log(U+1))。
    template <class L, class R> T get(ClosedSlice<L, R> s) const {
        return get(resolve_index(len(), s.a), resolve_index(len(), s.b) + 1);
    }

    // 閉区間の値を現在の並び順で集約する。キー上限UとしてO(log(N+1)+log(U+1))。
    template <class L, class R> T operator[](ClosedSlice<L, R> s) const {
        return get(s);
    }

    // 閉区間をキー順に並べ替える。キー上限Uとして償却O(log(N+1)+log(U+1))。
    template <class L, class R> void sort(ClosedSlice<L, R> s, SortOrder order = Ascending) {
        sort(resolve_index(len(), s.a), resolve_index(len(), s.b) + 1, order);
    }
};

// 要素数Nのキーと値から区間ソートとモノイド積を扱う木を構築する。演算O(1)なら期待O(N)時間・O(N)空間。
template <class Values, class Op, class T>
auto initRangeSortSegmentTree(std::span<const Int> keys, const Values &values, Int limit, Op op,
                              T e) {
    return RangeSortSegmentTree<T>(keys, values, limit, op, e);
}
}
