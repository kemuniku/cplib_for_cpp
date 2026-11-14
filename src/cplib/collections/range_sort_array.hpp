#pragma once
#include <cplib/collections/private/range_sort_impl.hpp>

namespace cplib {
class RangeSortArray : public detail::RangeSortEngine<detail::RangeSortEmpty, false> {
    using Base = detail::RangeSortEngine<detail::RangeSortEmpty, false>;

public:
    // 相異なる[0,limit)のキー列から構築する。要素数Nとして期待O(N)時間・O(N)空間。
    RangeSortArray(std::span<const Int> keys, Int limit)
        : Base(keys, {}, limit, [](auto, auto) { return detail::RangeSortEmpty{}; }, {}) {
    }

    // 現在の位置iのキーを返す。キー上限をUとしてO(log(N+1)+log(U+1))。
    Int get(Int i) const {
        return key(i);
    }

    // 現在の位置iのキーを返す。キー上限をUとしてO(log(N+1)+log(U+1))。
    Int operator[](Int i) const {
        return get(i);
    }

    // 末尾から指定した位置のキーを返す。キー上限をUとしてO(log(N+1)+log(U+1))。
    Int operator[](BackwardsIndex i) const {
        return get(size() - i.value);
    }

    // 位置iのキーを変更する。他のキーとの重複は禁止。キー上限Uとして期待O(log(N+1)+log(U+1))。
    void update(Int i, Int key) {
        Base::update(i, key, {});
    }

    // 位置iのキーを変更する。他のキーとの重複は禁止。キー上限Uとして期待O(log(N+1)+log(U+1))。
    void set(Int i, Int key) {
        update(i, key);
    }

    // 末尾から指定した位置のキーを変更する。他のキーとの重複は禁止。キー上限Uとして期待O(log(N+1)+log(U+1))。
    void set(BackwardsIndex i, Int key) {
        update(size() - i.value, key);
    }

    using Base::sort;

    // 閉区間をキー順に並べ替える。キー上限Uとして償却O(log(N+1)+log(U+1))。
    template <class L, class R> void sort(ClosedSlice<L, R> s, SortOrder order = Ascending) {
        sort(resolve_index(size(), s.a), resolve_index(size(), s.b) + 1, order);
    }

    // 現在のキー列を配列で返す。O(N)。
    std::vector<Int> toSeq() const {
        std::vector<Int> out;
        out.reserve(size());
        for (Int k : items())
            out.push_back(k);
        return out;
    }
};

// 相異なる[0,limit)のキー列から構築する。要素数Nとして期待O(N)時間・O(N)空間。
inline auto initRangeSortArray(std::span<const Int> keys, Int limit) {
    return RangeSortArray(keys, limit);
}
}
