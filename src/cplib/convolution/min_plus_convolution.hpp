#pragma once
#include <cplib/utils/monotone_minima.hpp>
#include <cplib/utils/smawk.hpp>
#include <type_traits>

namespace cplib {
namespace detail {
// b - a <= d - c を整数のオーバーフローを避けて判定する。
template <class T> bool differenceLeq(T a, T b, T c, T d) {
    if constexpr (std::is_integral_v<T>) {
        if (b < a) {
            if (d >= c)
                return true;
            return UInt(a) - UInt(b) >= UInt(c) - UInt(d);
        }
        if (d < c)
            return false;
        return UInt(b) - UInt(a) <= UInt(d) - UInt(c);
    } else
        return b - a <= d - c;
}

template <bool Convex, class T> void checkConvolutionShape(const std::vector<T> &a) {
#ifndef NDEBUG
    for (std::size_t i = 2; i < a.size(); ++i) {
        if constexpr (Convex)
            assert(differenceLeq(a[i - 2], a[i - 1], a[i - 1], a[i]));
        else
            assert(differenceLeq(a[i - 1], a[i], a[i - 2], a[i - 1]));
    }
#else
    (void)a;
#endif
}

template <bool UseSmawk, class T>
std::vector<T> convexArbitrary(const std::vector<T> &a, const std::vector<T> &b) {
    checkConvolutionShape<true>(a);
    if (a.empty() || b.empty())
        return {};
    Int n = a.size(), m = b.size(), h = n + m - 1;
    std::vector<T> out(h);
    if constexpr (UseSmawk)
        if (n <= 256 / m) {
            for (Int k = 0; k < h; ++k) {
                Int lo = std::max(Int(0), k - m + 1), hi = std::min(k, n - 1);
                T value = a[lo] + b[k - lo];
                for (Int i = lo + 1; i <= hi; ++i)
                    value = std::min(value, a[i] + b[k - i]);
                out[k] = value;
            }
            return out;
        }
    auto better = [&](Int row, Int oldCol, Int newCol) {
        if (newCol > row)
            return false;
        if (oldCol > row)
            return true;
        if (oldCol < row - n + 1)
            return true;
        if (newCol < row - n + 1)
            return false;
        return a[row - newCol] + b[newCol] < a[row - oldCol] + b[oldCol];
    };
    std::vector<Int> indices;
    if constexpr (UseSmawk)
        indices = smawk(h, m, better);
    else
        indices = monotoneMinima(h, m, better);
    for (Int k = 0; k < h; ++k)
        out[k] = a[k - indices[k]] + b[indices[k]];
    return out;
}

template <class T>
void concaveSmawk(const std::vector<T> &a, const std::vector<T> &b, Int first, Int step, Int count,
                  Int offset, Int width, std::vector<Int> &columns, std::vector<Int> &indices) {
    Int reduced = offset + width, size = 0;
    for (Int p = offset; p < reduced; ++p) {
        Int col = columns[p];
        while (size > 0) {
            Int row = first + (size - 1) * step, old = columns[reduced + size - 1];
            if (!(a[row - col] + b[col] < a[row - old] + b[old]))
                break;
            --size;
        }
        if (size < count)
            columns[reduced + size++] = col;
    }
    if (count > 1)
        concaveSmawk(a, b, first + step, step * 2, count / 2, reduced, size, columns, indices);
    Int left = 0;
    for (Int i = 0; i < count; i += 2) {
        Int row = first + i * step, right = size - 1;
        if (i + 1 < count) {
            right = left;
            while (columns[reduced + right] != indices[row + step])
                ++right;
        }
        Int best = columns[reduced + left];
        T value = a[row - best] + b[best];
        for (Int p = left + 1; p <= right; ++p) {
            Int col = columns[reduced + p];
            T candidate = a[row - col] + b[col];
            if (candidate < value) {
                best = col;
                value = candidate;
            }
        }
        indices[row] = best;
        left = right;
    }
}

template <class T>
void concaveDivide(const std::vector<T> &a, const std::vector<T> &b, Int top, Int bottom, Int left,
                   Int right, std::vector<T> &answer, std::vector<Int> &columns,
                   std::vector<Int> &indices) {
    Int n = a.size(), t = std::max(top, left), d = std::min(bottom, right + n - 1),
        l = std::max(left, t - n + 1), r = std::min(right, d);
    if (t >= d || l >= r)
        return;
    if (d - t <= 1024 / (r - l)) {
        for (Int row = t; row < d; ++row) {
            T value = answer[row];
            for (Int col = std::max(l, row - n + 1); col < std::min(r, row + 1); ++col)
                value = std::min(value, a[row - col] + b[col]);
            answer[row] = value;
        }
    } else if (r - 1 <= t && d - 1 < l + n) {
        for (Int p = 0; p < r - l; ++p)
            columns[p] = r - 1 - p;
        concaveSmawk(a, b, t, 1, d - t, 0, r - l, columns, indices);
        for (Int row = t; row < d; ++row) {
            Int col = indices[row];
            answer[row] = std::min(answer[row], a[row - col] + b[col]);
        }
    } else if (d - t >= r - l) {
        Int mid = (t + d) / 2;
        concaveDivide(a, b, t, mid, l, r, answer, columns, indices);
        concaveDivide(a, b, mid, d, l, r, answer, columns, indices);
    } else {
        Int mid = (l + r) / 2;
        concaveDivide(a, b, t, d, l, mid, answer, columns, indices);
        concaveDivide(a, b, t, d, mid, r, answer, columns, indices);
    }
}
}

// 凸同士は差分のマージでO(N+M)。
template <class T>
std::vector<T> minPlusConvolutionConvexConvex(const std::vector<T> &a, const std::vector<T> &b) {
    detail::checkConvolutionShape<true>(a);
    detail::checkConvolutionShape<true>(b);
    if (a.empty() || b.empty())
        return {};
    std::vector<T> out(a.size() + b.size() - 1);
    std::size_t i = 0, j = 0;
    out[0] = a[0] + b[0];
    for (std::size_t k = 1; k < out.size(); ++k) {
        if (j + 1 == b.size() || (i + 1 < a.size() && a[i + 1] + b[j] < a[i] + b[j + 1]))
            ++i;
        else
            ++j;
        out[k] = a[i] + b[j];
    }
    return out;
}

// 2種類の探索を独立に提供。MonotoneMinima版O((N+M)log(N+M))、SMAWK版O(N+M)。
// 凸な a と任意の b の min-plus 畳み込み。時間 O((N + M) log(N + M))。
template <class T>
auto minPlusConvolutionConvexArbitraryMonotoneMinima(const std::vector<T> &a,
                                                     const std::vector<T> &b) {
    return detail::convexArbitrary<false>(a, b);
}

// 凸な a と任意の b の min-plus 畳み込み。時間 O(N + M)。
template <class T>
auto minPlusConvolutionConvexArbitrarySmawk(const std::vector<T> &a, const std::vector<T> &b) {
    return detail::convexArbitrary<true>(a, b);
}

// 凹同士は候補区間の両端を比較する。O(N+M)。
template <class T>
std::vector<T> minPlusConvolutionConcaveConcave(const std::vector<T> &a, const std::vector<T> &b) {
    detail::checkConvolutionShape<false>(a);
    detail::checkConvolutionShape<false>(b);
    if (a.empty() || b.empty())
        return {};
    Int n = a.size(), m = b.size();
    std::vector<T> out(n + m - 1);
    for (Int k = 0; k < Int(out.size()); ++k) {
        Int lo = std::max(Int(0), k - m + 1), hi = std::min(k, n - 1);
        out[k] = std::min(a[lo] + b[k - lo], a[hi] + b[k - hi]);
    }
    return out;
}

// 有効領域を長方形に分割し、反転した列のSMAWKを適用。O((N+M)log(N+M))。
template <class T>
std::vector<T> minPlusConvolutionConcaveArbitrary(const std::vector<T> &a,
                                                  const std::vector<T> &b) {
    detail::checkConvolutionShape<false>(a);
    if (a.empty() || b.empty())
        return {};
    Int n = a.size(), m = b.size(), h = n + m - 1;
    std::vector<T> out(h);
    for (Int k = 0; k < h; ++k) {
        Int j = std::min(k, m - 1);
        out[k] = a[k - j] + b[j];
    }
    if (n == 1 || m == 1)
        return out;
    std::vector<Int> columns(m + 2 * h), indices(h);
    detail::concaveDivide(a, b, 0, h, 0, m, out, columns, indices);
    return out;
}
}
