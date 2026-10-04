#pragma once
#include <cplib/common.hpp>
#include <iterator>

namespace cplib::detail {
// 値を遅延変換する範囲。構築O(1)、各要素は関数評価1回、追加領域O(1)。
template <class Range, class Function> class TransformRange {
    Range range;
    Function function;

public:
    TransformRange(Range r, Function f) : range(std::move(r)), function(std::move(f)) {
    }

    struct iterator {
        decltype(std::begin(std::declval<const Range &>())) current;
        const Function *function;

        auto operator*() const {
            return (*function)(*current);
        }

        iterator &operator++() {
            ++current;
            return *this;
        }

        bool operator==(const iterator &other) const {
            return current == other.current;
        }
    };

    auto begin() const {
        return iterator{std::begin(range), &function};
    }

    auto end() const {
        return iterator{std::end(range), &function};
    }
};

template <class Range, class Function> auto transform_range(Range r, Function f) {
    return TransformRange<Range, Function>(std::move(r), std::move(f));
}

// 条件を満たす要素だけを遅延列挙する。走査O(元の要素数)、追加領域O(1)。
template <class Range, class Predicate> class FilterRange {
    Range range;
    Predicate predicate;

public:
    FilterRange(Range r, Predicate p) : range(std::move(r)), predicate(std::move(p)) {
    }

    struct iterator {
        decltype(std::begin(std::declval<const Range &>())) current, last;
        const Predicate *predicate;

        void skip() {
            while (current != last && !(*predicate)(*current))
                ++current;
        }

        auto operator*() const {
            return *current;
        }

        iterator &operator++() {
            ++current;
            skip();
            return *this;
        }

        bool operator==(const iterator &other) const {
            return current == other.current;
        }
    };

    auto begin() const {
        iterator it{std::begin(range), std::end(range), &predicate};
        it.skip();
        return it;
    }

    auto end() const {
        return iterator{std::end(range), std::end(range), &predicate};
    }
};

template <class Range, class Predicate> auto filter_range(Range r, Predicate p) {
    return FilterRange<Range, Predicate>(std::move(r), std::move(p));
}
}
