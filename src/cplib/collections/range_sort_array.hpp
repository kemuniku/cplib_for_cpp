#pragma once
#include <cplib/collections/private/range_sort_impl.hpp>

namespace cplib {
class RangeSortArray : public detail::RangeSortEngine<detail::RangeSortEmpty, false> {
    using Base = detail::RangeSortEngine<detail::RangeSortEmpty, false>;

public:
    RangeSortArray(std::span<const Int> keys, Int limit)
        : Base(keys, {}, limit, [](auto, auto) { return detail::RangeSortEmpty{}; }, {}) {
    }

    Int get(Int i) const {
        return key(i);
    }

    Int operator[](Int i) const {
        return get(i);
    }

    Int operator[](BackwardsIndex i) const {
        return get(len() - i.value);
    }

    void update(Int i, Int key) {
        Base::update(i, key, {});
    }

    void set(Int i, Int key) {
        update(i, key);
    }

    void set(BackwardsIndex i, Int key) {
        update(len() - i.value, key);
    }

    using Base::sort;

    template <class L, class R> void sort(ClosedSlice<L, R> s, SortOrder order = Ascending) {
        sort(resolve_index(len(), s.a), resolve_index(len(), s.b) + 1, order);
    }

    std::vector<Int> toSeq() const {
        std::vector<Int> out;
        out.reserve(len());
        for (Int k : items())
            out.push_back(k);
        return out;
    }
};

inline auto initRangeSortArray(std::span<const Int> keys, Int limit) {
    return RangeSortArray(keys, limit);
}
}
