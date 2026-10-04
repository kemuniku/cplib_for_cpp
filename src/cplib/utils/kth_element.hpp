#pragma once
#include <cplib/common.hpp>
#include <random>

namespace cplib {
namespace detail {
inline thread_local std::mt19937_64 selection_rng{std::random_device{}()};
}

// 入力をコピーして乱択選択する。期待O(n)時間・O(n)領域。
template <class Range> auto kth_element(const Range &input, Int k) {
    using T = typename Range::value_type;
    std::vector<T> x(input.begin(), input.end());
    Int now = 0, r = x.size();
    assert(k >= 0 && k < r);
    std::vector<int> c(r);
    for (;;) {
        Int l = 0, s = 1;
        std::swap(x[std::uniform_int_distribution<Int>(0, r - 1)(detail::selection_rng)], x[r - 1]);
        c[r - 1] = 0;
        for (Int i = 0; i < r - 1; ++i) {
            int tmp = x[i] < x[r - 1] ? -1 : x[r - 1] < x[i] ? 1 : 0;
            if (tmp < 0)
                ++l;
            else if (!tmp)
                ++s;
            c[i] = tmp;
        }
        if (now + l > k) {
            Int idx = 0;
            for (Int i = 0; i < r - 1; ++i)
                if (c[i] < 0)
                    std::swap(x[i], x[idx++]);
            r = idx;
        } else if (now + l + s > k)
            return x[r - 1];
        else {
            Int idx = 0;
            for (Int i = 0; i < r - 1; ++i)
                if (c[i] > 0)
                    std::swap(x[i], x[idx++]);
            now += l + s;
            r = idx;
        }
    }
}

// 3-way partitionで入力を並べ替えて乱択選択する。期待O(n)時間・O(1)追加領域。
// 引数numsの要素の順番が変わる。
template <class T> T kth_element_break(std::vector<T> &nums, Int k) {
    assert(k >= 0 && k < Int(nums.size()));
    Int left = 0, right = Int(nums.size()) - 1;
    while (left <= right) {
        if (left == right)
            return nums[left];
        T pivot = nums[std::uniform_int_distribution<Int>(left, right)(detail::selection_rng)];
        Int lt = left, i = left, gt = right;
        while (i <= gt) {
            if (nums[i] < pivot)
                std::swap(nums[lt++], nums[i++]);
            else if (pivot < nums[i])
                std::swap(nums[i], nums[gt--]);
            else
                ++i;
        }
        if (lt <= k && k <= gt)
            return nums[lt];
        else if (k < lt)
            right = lt - 1;
        else
            left = gt + 1;
    }
    std::abort();
}
}
