#pragma once
#include <cplib/graph/graph.hpp>
#include <cplib/tree/prufer.hpp>
#include <coroutine>
#include <exception>
#include <optional>
#include <functional>
#include <array>
#include <iterator>

namespace cplib {
// C++20の逐次生成範囲。yieldした値は次の反復まで有効、途中終了時も探索状態を解放する。
template <class T> class Generator {
public:
    struct promise_type {
        const T *current = nullptr;
        std::exception_ptr error;

        Generator get_return_object() {
            return Generator(std::coroutine_handle<promise_type>::from_promise(*this));
        }

        std::suspend_always initial_suspend() noexcept {
            return {};
        }

        std::suspend_always final_suspend() noexcept {
            return {};
        }

        std::suspend_always yield_value(const T &value) noexcept {
            current = &value;
            return {};
        }

        void return_void() noexcept {
        }

        void unhandled_exception() {
            error = std::current_exception();
        }
    };

    using Handle = std::coroutine_handle<promise_type>;

    struct Iterator {
        using value_type = T;
        using difference_type = std::ptrdiff_t;
        using iterator_category = std::input_iterator_tag;
        Handle handle;

        const T &operator*() const {
            return *handle.promise().current;
        }

        const T *operator->() const {
            return handle.promise().current;
        }

        Iterator &operator++() {
            handle.resume();
            if (handle.done() && handle.promise().error)
                std::rethrow_exception(handle.promise().error);
            return *this;
        }

        void operator++(int) {
            ++*this;
        }

        bool operator==(std::default_sentinel_t) const {
            return !handle || handle.done();
        }
    };

    Generator() = default;

    explicit Generator(Handle h) : handle_(h) {
    }

    Generator(const Generator &) = delete;
    Generator &operator=(const Generator &) = delete;

    Generator(Generator &&b) noexcept : handle_(std::exchange(b.handle_, {})) {
    }

    Generator &operator=(Generator &&b) noexcept {
        if (this != &b) {
            if (handle_)
                handle_.destroy();
            handle_ = std::exchange(b.handle_, {});
        }
        return *this;
    }

    ~Generator() {
        if (handle_)
            handle_.destroy();
    }

    Iterator begin() {
        Iterator it{handle_};
        if (handle_)
            ++it;
        return it;
    }

    std::default_sentinel_t end() const {
        return {};
    }

private:
    Handle handle_{};
};

namespace detail::itertools {
// 左辺値は参照し、一時入力は生成器が所有する。入力の寿命と不要な全件コピーを両立する。
template <class R> struct Input {
    using Base = std::remove_reference_t<R>;
    using value_type = typename std::remove_cv_t<Base>::value_type;
    using Store =
        std::conditional_t<std::is_lvalue_reference_v<R>, std::reference_wrapper<Base>, Base>;
    Store store;

    explicit Input(R &&r) : store(std::forward<R>(r)) {
    }

    const Base &get() const {
        if constexpr (std::is_lvalue_reference_v<R>)
            return store.get();
        else
            return store;
    }

    Int size() const {
        return get().size();
    }

    decltype(auto) operator[](Int i) const {
        return get()[i];
    }
};

template <class R> auto hold(R &&r) {
    return Input<R>(std::forward<R>(r));
}
template <class I> using Value = typename I::value_type;

template <class I> Generator<std::vector<Value<I>>> permutations(I v) {
    std::vector<Int> idx(v.size());
    std::iota(idx.begin(), idx.end(), 0);
    std::vector<Value<I>> out(v.size());
    do {
        for (Int i = 0; i < v.size(); ++i)
            out[i] = v[idx[i]];
        co_yield out;
    } while (std::next_permutation(idx.begin(), idx.end()));
}

template <class I> Generator<std::vector<Value<I>>> distinct_permutations(I v) {
    std::vector<Value<I>> out;
    for (Int i = 0; i < v.size(); ++i)
        out.push_back(v[i]);
    std::sort(out.begin(), out.end());
    do {
        co_yield out;
    } while (std::next_permutation(out.begin(), out.end()));
}

template <class I> Generator<std::vector<Value<I>>> combinations(I v, Int r) {
    Int n = v.size();
    if (r == 0) {
        co_yield std::vector<Value<I>>{};
        co_return;
    }
    if (r < 0 || r > n)
        co_return;
    std::vector<Int> idx(r);
    std::iota(idx.begin(), idx.end(), 0);
    std::vector<Value<I>> out(r);
    while (true) {
        for (Int i = 0; i < r; ++i)
            out[i] = v[idx[i]];
        co_yield out;
        Int i = r - 1;
        while (i >= 0 && idx[i] == i + n - r)
            --i;
        if (i < 0)
            break;
        ++idx[i];
        for (Int j = i + 1; j < r; ++j)
            idx[j] = idx[j - 1] + 1;
    }
}

template <Int R, class I> Generator<std::array<Value<I>, R>> combinations(I v) {
    static_assert(R >= 0);
    Int n = v.size();
    std::array<Value<I>, R> out{};
    if constexpr (R == 0) {
        co_yield out;
    } else if (R <= n) {
        std::array<Int, R> idx{};
        for (Int i = 0; i < R; ++i) {
            idx[i] = i;
            out[i] = v[i];
        }
        co_yield out;
        while (true) {
            Int i = R - 1;
            while (i >= 0 && idx[i] == i + n - R)
                --i;
            if (i < 0)
                break;
            ++idx[i];
            out[i] = v[idx[i]];
            for (Int j = i + 1; j < R; ++j) {
                idx[j] = idx[j - 1] + 1;
                out[j] = v[idx[j]];
            }
            co_yield out;
        }
    }
}

template <class I, class F> Generator<Value<I>> combinations_withf(I v, Int r, F f) {
    Int n = v.size();
    assert(r >= 0);
    if (r == 0 || r > n)
        co_return;
    std::vector<Int> idx(r);
    std::iota(idx.begin(), idx.end(), 0);
    std::vector<Value<I>> pref(r);
    pref[0] = v[0];
    for (Int i = 1; i < r; ++i)
        pref[i] = f(pref[i - 1], v[i]);
    co_yield pref.back();
    while (true) {
        Int i = r - 1;
        while (i >= 0 && idx[i] == i + n - r)
            --i;
        if (i < 0)
            break;
        ++idx[i];
        pref[i] = i == 0 ? v[idx[i]] : f(pref[i - 1], v[idx[i]]);
        for (Int j = i + 1; j < r; ++j) {
            idx[j] = idx[j - 1] + 1;
            pref[j] = f(pref[j - 1], v[idx[j]]);
        }
        co_yield pref.back();
    }
}

template <class I> Generator<std::vector<Value<I>>> product(I v, Int repeat) {
    assert(repeat >= 0);
    if (!repeat) {
        co_yield std::vector<Value<I>>{};
        co_return;
    }
    if (!v.size())
        co_return;
    std::vector<Int> idx(repeat);
    std::vector<Value<I>> out(repeat);
    bool more = true;
    while (more) {
        for (Int i = 0; i < repeat; ++i)
            out[i] = v[idx[i]];
        co_yield out;
        for (Int i = 0; i < repeat; ++i) {
            if (++idx[i] == v.size()) {
                idx[i] = 0;
                if (i == repeat - 1)
                    more = false;
            } else
                break;
        }
    }
}

template <class I> Generator<std::vector<Int>> bounded_sequences(I limits) {
    for (Int i = 0; i < limits.size(); ++i)
        assert(limits[i] >= 0);
    std::vector<Int> out(limits.size());
    while (true) {
        co_yield out;
        Int i = limits.size() - 1;
        while (i >= 0 && out[i] == limits[i])
            out[i--] = 0;
        if (i < 0)
            break;
        ++out[i];
    }
}

template <class I> Generator<std::vector<Int>> bounded_sum_sequences(Int sum, I bounds) {
    Int n = bounds.size();
    std::vector<Int> lo(n + 1), hi(n + 1);
    for (Int i = n - 1; i >= 0; --i) {
        if (bounds[i].first >= bounds[i].second)
            co_return;
        lo[i] = lo[i + 1] + bounds[i].first;
        hi[i] = hi[i + 1] + bounds[i].second - 1;
    }
    if (sum < lo[0] || sum > hi[0])
        co_return;
    if (!n) {
        co_yield std::vector<Int>{};
        co_return;
    }
    std::vector<Int> out(n), remaining(n + 1), upper(n);
    remaining[0] = sum;
    Int depth = 0;
    out[0] = std::max(Int(bounds[0].first), sum - hi[1]);
    upper[0] = std::min(Int(bounds[0].second) - 1, sum - lo[1]);
    while (depth >= 0) {
        if (out[depth] > upper[depth]) {
            if (--depth >= 0)
                ++out[depth];
        } else if (depth == n - 1) {
            co_yield out;
            ++out[depth];
        } else {
            remaining[depth + 1] = remaining[depth] - out[depth];
            ++depth;
            out[depth] = std::max(Int(bounds[depth].first), remaining[depth] - hi[depth + 1]);
            upper[depth] =
                std::min(Int(bounds[depth].second) - 1, remaining[depth] - lo[depth + 1]);
        }
    }
}

inline Generator<std::vector<Int>> monotone(Int n, Int l, Int r, Int step,
                                            std::optional<Int> target) {
    assert(n >= 0);
    if (n == 0) {
        if (!target || *target == 0)
            co_yield std::vector<Int>{};
        co_return;
    }
    if (l >= r || (step && n > r - l))
        co_return;
    std::vector<Int> out(n), prefix(n + 1);
    Int depth = 0;
    out[0] = l;
    while (depth >= 0) {
        Int count = n - depth, upper = r - 1 - step * (count - 1);
        bool exhausted = out[depth] > upper;
        if (!exhausted && target) {
            Int minimum = prefix[depth] + count * out[depth] + step * (count * (count - 1) / 2);
            Int maximum = prefix[depth] + out[depth] + (count - 1) * (r - 1) -
                          step * ((count - 1) * (count - 2) / 2);
            if (minimum > *target)
                exhausted = true;
            else if (maximum < *target) {
                out[depth] += *target - maximum;
                continue;
            }
        }
        if (exhausted) {
            if (--depth >= 0)
                ++out[depth];
        } else if (depth == n - 1) {
            co_yield out;
            ++out[depth];
        } else {
            prefix[depth + 1] = prefix[depth] + out[depth];
            out[depth + 1] = out[depth] + step;
            ++depth;
        }
    }
}

template <class I>
Generator<std::vector<typename Value<I>::value_type>> cartesian_product(I choices) {
    using T = typename Value<I>::value_type;
    Int n = choices.size();
    for (Int i = 0; i < n; ++i)
        if (choices[i].empty())
            co_return;
    std::vector<Int> idx(n);
    std::vector<T> out(n);
    while (true) {
        for (Int i = 0; i < n; ++i)
            out[i] = choices[i][idx[i]];
        co_yield out;
        Int i = n - 1;
        while (i >= 0 && idx[i] == Int(choices[i].size()) - 1)
            idx[i--] = 0;
        if (i < 0)
            break;
        ++idx[i];
    }
}
}

// 順列・組合せ・直積を遅延列挙する。固定長組合せと差分畳み込みも別経路で維持する。
// 入力位置の順列をPython itertools.permutationsと同じ順で列挙する。入力に同値要素があれば値の重複を保持する。
template <class R> auto permutations(R &&v) {
    return detail::itertools::permutations(detail::itertools::hold(std::forward<R>(v)));
}

template <class R> auto distinct_permutations(R &&v) {
    return detail::itertools::distinct_permutations(detail::itertools::hold(std::forward<R>(v)));
}

template <class R> auto combinations(R &&v, Int r) {
    return detail::itertools::combinations(detail::itertools::hold(std::forward<R>(v)), r);
}

template <Int N, class R> auto combinations(R &&v) {
    return detail::itertools::combinations<N>(detail::itertools::hold(std::forward<R>(v)));
}

// combinations(v, r)の各組合せをfで左畳み込みした値を列挙する。fが乗算なら各組合せの総積。
// 接頭辞畳み込みを保持し、変化のあった添字以降のみ再計算する差分更新版。
template <class R, class F> auto combinations_withf(R &&v, Int r, F f) {
    return detail::itertools::combinations_withf(detail::itertools::hold(std::forward<R>(v)), r,
                                                 std::move(f));
}

template <class R> auto product(R &&v, Int repeat) {
    return detail::itertools::product(detail::itertools::hold(std::forward<R>(v)), repeat);
}

// 各位置の候補から1つずつ選ぶ。右端から候補の添字を進める。
// 位置が0個なら空列を1件、空の候補があれば0件。候補内の重複は保持。
template <class R> auto cartesian_product(R &&v) {
    return detail::itertools::cartesian_product(detail::itertools::hold(std::forward<R>(v)));
}

// 左・右累積演算。初期値ありは長さn+1、なしは長さn。O(n)。
template <class T, class F> std::vector<T> accumulated(const std::vector<T> &v, F f, T first) {
    std::vector<T> out(v.size() + 1);
    out[0] = first;
    for (std::size_t i = 0; i < v.size(); ++i)
        out[i + 1] = f(out[i], v[i]);
    return out;
}

template <class T, class F> void accumulate(std::vector<T> &v, F f) {
    for (std::size_t i = 1; i < v.size(); ++i)
        v[i] = f(v[i - 1], v[i]);
}

template <class T, class F> std::vector<T> accumulated(std::vector<T> v, F f) {
    cplib::accumulate(v, f);
    return v;
}

template <class T, class F> std::vector<T> accumulatedr(const std::vector<T> &v, F f, T first) {
    std::vector<T> out(v.size() + 1);
    out.back() = first;
    for (Int i = Int(v.size()) - 1; i >= 0; --i)
        out[i] = f(v[i], out[i + 1]);
    return out;
}

template <class T, class F> void accumulater(std::vector<T> &v, F f) {
    for (Int i = Int(v.size()) - 2; i >= 0; --i)
        v[i] = f(v[i], v[i + 1]);
}

template <class T, class F> std::vector<T> accumulatedr(std::vector<T> v, F f) {
    cplib::accumulater(v, f);
    return v;
}

// 非減少の整数分割を元と同じ順序で列挙する。追加領域O(n)。
inline Generator<std::vector<Int>> partitions(Int n) {
    assert(n >= 0);
    if (!n) {
        co_yield std::vector<Int>{};
        co_return;
    }
    std::vector<Int> a(n + 1);
    Int k = 1;
    a[1] = n;
    while (k) {
        Int x = a[k - 1] + 1, y = a[k] - 1;
        --k;
        while (x <= y) {
            a[k++] = x;
            y -= x;
        }
        a[k] = x + y;
        co_yield std::vector<Int>(a.begin(), a.begin() + k + 1);
    }
}

// 非負整数列upperに対し、0 <= b[i] <= upper[i]を満たす列を辞書順に列挙する。
// 空入力では空列を1件返す。1件あたりO(upper.size())、追加領域O(upper.size())。
template <class R> auto bounded_sequences(R &&upper) {
    return detail::itertools::bounded_sequences(detail::itertools::hold(std::forward<R>(upper)));
}

// 総和sum、bounds[i].first <= a[i] < bounds[i].secondの数列を辞書順に列挙する。
// 空のboundsはsum == 0のときだけ空列を返す。空の範囲は解なし。
// 各範囲の端点・累積和・差の計算はIntに収まること。
template <class R> auto bounded_sum_sequences(Int sum, R &&bounds) {
    return detail::itertools::bounded_sum_sequences(
        sum, detail::itertools::hold(std::forward<R>(bounds)));
}

// 長さn、総和sum、各要素が[l,r)の数列を辞書順に列挙する。
inline auto bounded_sum_sequences(Int n, Int sum, Int l, Int r) {
    assert(n >= 0);
    return bounded_sum_sequences(sum, std::vector<std::pair<Int, Int>>(n, {l, r}));
}

// 長さ n、各要素が [l, r) の広義単調増加列を辞書順に列挙。
inline auto nondecreasing_sequences(Int n, Int l, Int r) {
    return detail::itertools::monotone(n, l, r, 0, std::nullopt);
}

// 総和sを指定する版。端点・総和の中間計算はIntに収まること。
inline auto nondecreasing_sequences(Int n, Int s, Int l, Int r) {
    return detail::itertools::monotone(n, l, r, 0, s);
}

// 長さ n、各要素が [l, r) の狭義単調増加列を辞書順に列挙。
inline auto strictly_increasing_sequences(Int n, Int l, Int r) {
    return detail::itertools::monotone(n, l, r, 1, std::nullopt);
}

// 総和sを指定する狭義単調増加列の版。端点・総和の中間計算はIntに収まること。
inline auto strictly_increasing_sequences(Int n, Int s, Int l, Int r) {
    return detail::itertools::monotone(n, l, r, 1, s);
}

// 初出順のグループ番号で集合分割を列挙する。k=-1なら個数制限なし。
inline Generator<std::vector<Int>> set_partitions_id(Int n, Int k = -1) {
    assert(n >= 0 && k >= -1);
    if (n == 0) {
        if (k == -1 || k == 0)
            co_yield std::vector<Int>{};
        co_return;
    }
    if (k == 0 || k > n)
        co_return;
    std::vector<Int> a(n), groups(n + 1);
    Int depth = 0;
    while (depth >= 0) {
        Int upper = groups[depth];
        if (k >= 0)
            upper = std::min(upper, k - 1);
        if (a[depth] > upper) {
            if (--depth >= 0)
                ++a[depth];
        } else {
            groups[depth + 1] = std::max(groups[depth], a[depth] + 1);
            if (k >= 0 && groups[depth + 1] + n - depth - 1 < k)
                ++a[depth];
            else if (depth == n - 1) {
                co_yield a;
                ++a[depth];
            } else
                a[++depth] = 0;
        }
    }
}

// グループそのものを差分更新し、探索で空になった配列の容量も再利用する。
// 0..<n の集合分割を、各グループの要素の列で返す。k == -1 は個数指定なし。
// 各グループ内は昇順、グループ間は最小要素の昇順。列挙順は set_partitions_id と同じ。
// n == 0 は k == -1 または k == 0 のときだけ空列を1件返す。
inline Generator<std::vector<std::vector<Int>>> set_partitions(Int n, Int k = -1) {
    assert(n >= 0 && k >= -1);
    if (n == 0) {
        if (k == -1 || k == 0)
            co_yield std::vector<std::vector<Int>>{};
        co_return;
    }
    if (k == 0 || k > n)
        co_return;
    std::vector<std::vector<Int>> groups, spare(n);
    groups.reserve(n);
    std::vector<Int> choice(n);
    Int depth = 0;
    while (depth >= 0) {
        Int upper = groups.size();
        if (k >= 0)
            upper = std::min(upper, k - 1);
        if (choice[depth] > upper) {
            if (depth == 0)
                break;
            --depth;
        } else {
            Int count = std::max(Int(groups.size()), choice[depth] + 1);
            if (k >= 0 && count + n - depth - 1 < k) {
                ++choice[depth];
                continue;
            }
            if (choice[depth] == Int(groups.size())) {
                groups.emplace_back();
                groups.back().swap(spare[depth]);
            }
            groups[choice[depth]].push_back(depth);
            if (depth == n - 1) {
                co_yield groups;
            } else {
                choice[++depth] = 0;
                continue;
            }
        }
        Int id = choice[depth];
        groups[id].pop_back();
        if (groups[id].empty()) {
            spare[depth] = std::move(groups[id]);
            groups.pop_back();
        }
        ++choice[depth];
    }
}

// 最小の未使用要素から相手を選び、順序による重複のないペア分割を列挙する。
// n == 0は空列を1件、奇数なら0件。
inline Generator<std::vector<std::pair<Int, Int>>> pairings(Int n) {
    assert(n >= 0);
    if (!n) {
        co_yield std::vector<std::pair<Int, Int>>{};
        co_return;
    }
    if (n % 2)
        co_return;
    std::vector<bool> used(n);
    std::vector<std::pair<Int, Int>> pairs(n / 2);
    Int depth = 0;
    pairs[0] = {0, 0};
    used[0] = true;
    while (depth >= 0) {
        auto &[u, v] = pairs[depth];
        ++v;
        while (v < n && used[v])
            ++v;
        if (v == n) {
            used[u] = false;
            if (--depth >= 0)
                used[pairs[depth].second] = false;
        } else if (depth == Int(pairs.size()) - 1) {
            co_yield pairs;
        } else {
            used[v] = true;
            Int next = 0;
            while (used[next])
                ++next;
            pairs[++depth] = {next, next};
            used[next] = true;
        }
    }
}

// 非空連続区間の境界を列挙する。空列の場合の唯一の境界は{0}。
inline Generator<std::vector<Int>> contiguous_partitions(Int n, Int k = -1) {
    assert(n >= 0 && k >= -1);
    if (n == 0) {
        if (k == -1 || k == 0)
            co_yield std::vector<Int>(1, 0);
        co_return;
    }
    if (k == 0 || k > n)
        co_return;
    Int first = k == -1 ? 1 : k, last = k == -1 ? n : k;
    for (Int count = first; count <= last; ++count)
        for (const auto &cuts : strictly_increasing_sequences(count - 1, 1, n)) {
            std::vector<Int> out;
            out.reserve(count + 1);
            out.push_back(0);
            out.insert(out.end(), cuts.begin(), cuts.end());
            out.push_back(n);
            co_yield out;
        }
}

// 正しい括弧列を辞書順に列挙する。追加領域O(n)。
// n 組（長さ 2*n）の正しい括弧列を辞書順に列挙。n == 0 は空文字列。
inline Generator<std::string> parenthesis_sequences(Int n) {
    assert(n >= 0);
    if (n == 0) {
        co_yield std::string{};
        co_return;
    }
    std::string a(2 * n, '\0');
    std::vector<Int> balance(2 * n + 1), choice(2 * n);
    Int depth = 0;
    while (depth >= 0) {
        if (choice[depth] >= 2) {
            --depth;
        } else {
            bool opening = choice[depth]++ == 0;
            if (opening && (depth + balance[depth]) / 2 == n)
                continue;
            if (!opening && balance[depth] == 0)
                continue;
            a[depth] = opening ? '(' : ')';
            balance[depth + 1] = balance[depth] + (opening ? 1 : -1);
            if (depth == Int(a.size()) - 1) {
                co_yield a;
            } else
                choice[++depth] = 0;
        }
    }
}

// Prüfer列によるラベル付き木の列挙。元と同じ復号法を使用する。
// 頂点番号 0..<n の木を重複なく列挙。n >= 1。
// Prüfer 列を使用し、n >= 2 では n^(n-2) 件。
inline Generator<UnWeightedUnDirectedGraph> labeled_trees(Int n) {
    assert(n >= 1);
    if (n == 1) {
        co_yield initUnWeightedUnDirectedGraph(1);
    } else {
        std::vector<Int> vertices(n);
        std::iota(vertices.begin(), vertices.end(), 0);
        for (const auto &code : product(vertices, n - 2))
            co_yield prufer_decode(code);
    }
}

// 辺集合の組合せから単純無向グラフを列挙する。
// 頂点番号 0..<n、辺数 m の単純無向グラフ。m == -1 は辺数指定なし。
// 同型でも頂点番号が異なるグラフは区別する。
inline Generator<UnWeightedUnDirectedGraph> simple_graphs(Int n, Int m = -1) {
    assert(n >= 0 && m >= -1);
    std::vector<std::pair<Int, Int>> edges;
    for (Int u = 0; u < n; ++u)
        for (Int v = u + 1; v < n; ++v)
            edges.emplace_back(u, v);
    Int first = m == -1 ? 0 : m,
        last = m == -1 ? Int(edges.size()) : std::min(m, Int(edges.size()));
    for (Int count = first; count <= last; ++count)
        for (const auto &selected : combinations(edges, count)) {
            auto g = initUnWeightedUnDirectedGraph(n);
            for (auto [u, v] : selected)
                g.add_edge(u, v);
            co_yield g;
        }
}

namespace detail::itertools {
template <class I> Generator<std::vector<Int>> topological_orders(I adj) {
    Int n = adj.size();
    std::vector<Int> indegree(n);
    for (Int i = 0; i < n; ++i)
        for (Int v : adj[i]) {
            assert(v >= 0 && v < n);
            ++indegree[v];
        }
    if (!n) {
        co_yield std::vector<Int>{};
        co_return;
    }
    std::vector<bool> used(n);
    std::vector<Int> order(n), next(n);
    Int depth = 0;
    while (depth >= 0) {
        Int u = next[depth];
        while (u < n && (used[u] || indegree[u] != 0))
            ++u;
        if (u == n) {
            if (--depth >= 0) {
                Int previous = order[depth];
                used[previous] = false;
                for (Int v : adj[previous])
                    ++indegree[v];
            }
        } else {
            next[depth] = u + 1;
            order[depth] = u;
            if (depth == n - 1) {
                co_yield order;
            } else {
                used[u] = true;
                for (Int v : adj[u])
                    --indegree[v];
                next[++depth] = 0;
            }
        }
    }
}
}

// 入次数を巻き戻しながらトポロジカル順序を辞書順で列挙する。
// 隣接リストのトポロジカル順序を辞書順に列挙する。有向閉路がある場合は0件。
template <class R> auto topological_orders(R &&adj) {
    return detail::itertools::topological_orders(detail::itertools::hold(std::forward<R>(adj)));
}

// cplib の有向グラフ版。重みは無視。静的グラフは build 済みであること。
template <class T, bool S> auto topological_orders(const BasicGraph<T, S, true> &g) {
    std::vector<std::vector<Int>> adj(g.len);
    for (Int u = 0; u < g.len; ++u)
        for (auto [v, c] : g.to_and_cost(u)) {
            (void)c;
            adj[u].push_back(v);
        }
    return topological_orders(std::move(adj));
}

template <class T, bool S> auto topological_orders(BasicGraph<T, S, true> &g) {
    return topological_orders(std::as_const(g));
}

template <class T, bool S> auto topological_orders(BasicGraph<T, S, true> &&g) {
    return topological_orders(std::as_const(g));
}

// L1ノルムの上限から次の要素の範囲を絞り、辞書順に列挙する。
// 長さn、sum(abs(a[i])) <= sの整数列を辞書順に列挙する。
// s < 0は0件。sはIntの最大値未満であること。
inline Generator<std::vector<Int>> integer_vectors_l1(Int n, Int s) {
    assert(n >= 0 && s < std::numeric_limits<Int>::max());
    if (s < 0)
        co_return;
    if (!n) {
        co_yield std::vector<Int>{};
        co_return;
    }
    std::vector<Int> a(n), remaining(n + 1);
    remaining[0] = s;
    a[0] = -s;
    Int depth = 0;
    while (depth >= 0) {
        if (a[depth] > remaining[depth]) {
            if (--depth >= 0)
                ++a[depth];
        } else if (depth == n - 1) {
            co_yield a;
            ++a[depth];
        } else {
            remaining[depth + 1] = remaining[depth] - std::abs(a[depth]);
            ++depth;
            a[depth] = -remaining[depth];
        }
    }
}

template <class T, class U> struct Counterexample {
    T input;
    U actual, expected;
};

// 範囲または逐次生成器を先頭から比較し、最初の不一致で中断する。
template <class R, class F, class N>
    requires requires(R &r) {
        std::begin(r);
        std::end(r);
    }
// 最初にsolveとnaiveの戻り値が!=となる入力と両出力を返す。
// 全件一致ならstd::nullopt。コールバックは入力を変更しないこと。
auto find_counterexample(R &&cases, F solve, N naive) {
    using T = std::remove_cvref_t<decltype(*std::begin(cases))>;
    using U = std::remove_cvref_t<std::invoke_result_t<F, const T &>>;
    for (const auto &input : cases) {
        U actual = solve(input), expected = naive(input);
        if (actual != expected)
            return std::optional<Counterexample<T, U>>(
                {input, std::move(actual), std::move(expected)});
    }
    return std::optional<Counterexample<T, U>>{};
}

template <class G, class F, class N>
    requires std::invocable<G>
auto find_counterexample(G cases, F solve, N naive) {
    auto range = cases();
    return find_counterexample(range, std::move(solve), std::move(naive));
}

// 要素削除と0方向への二分縮小を、変化がなくなるまで貪欲に繰り返す。
// failsが真の整数列を、要素削除・0への置換・0方向への移動で縮小する。
// 変化するたびに先頭から再試行する貪欲法で、最小の反例とは限らない。
// failsは決定的で入力を変更しないこと。初期入力でも真であること。
template <class F> std::vector<Int> shrink_counterexample(std::vector<Int> input, F fails) {
    assert(fails(input));
    while (true) {
        bool changed = false;
        for (Int i = 0; i < Int(input.size()); ++i) {
            auto candidate = input;
            candidate.erase(candidate.begin() + i);
            if (fails(candidate)) {
                input = std::move(candidate);
                changed = true;
                break;
            }
        }
        if (changed)
            continue;
        for (Int i = 0; i < Int(input.size()); ++i) {
            Int value = input[i];
            if (!value)
                continue;
            Int delta = value;
            while (delta) {
                auto candidate = input;
                candidate[i] = value - delta;
                if (fails(candidate)) {
                    input = std::move(candidate);
                    changed = true;
                    break;
                }
                delta /= 2;
            }
            if (changed)
                break;
        }
        if (!changed)
            break;
    }
    return input;
}
}

namespace cplib {
namespace detail::itertools {
inline Generator<std::vector<Int>> rooted_tree_levels(Int n) {
    assert(n >= 1);
    std::vector<Int> levels(n);
    std::iota(levels.begin(), levels.end(), 0);
    while (true) {
        co_yield levels;
        Int p = n - 1;
        while (p > 0 && levels[p] == 1)
            --p;
        if (!p)
            break;
        Int q = p - 1;
        while (levels[q] != levels[p] - 1)
            --q;
        Int period = p - q;
        for (Int i = p; i < n; ++i)
            levels[i] = levels[i - period];
    }
}

inline UnWeightedUnDirectedGraph tree_from_levels(const std::vector<Int> &levels) {
    UnWeightedUnDirectedGraph g(levels.size());
    std::vector<Int> path(levels.size());
    for (Int v = 1; v < Int(levels.size()); ++v) {
        g.add_edge(path[levels[v] - 1], v);
        path[levels[v]] = v;
    }
    return g;
}

inline bool unrooted_representative(const std::vector<Int> &levels) {
    Int n = levels.size(), halfStart = -1;
    for (Int first = 1; first < n;) {
        Int last = first + 1;
        while (last < n && levels[last] > 1)
            ++last;
        Int size = last - first;
        if (size > n / 2)
            return false;
        if (n % 2 == 0 && size == n / 2)
            halfStart = first;
        first = last;
    }
    if (halfStart >= 0) {
        Int half = n / 2;
        for (Int i = 0; i < half; ++i) {
            Int index = i < halfStart ? i : i + half, left = levels[index],
                right = levels[halfStart + i] - 1;
            if (left != right)
                return left > right;
        }
    }
    return true;
}
}

// n 頂点の根付き木を、根を保つ同型を除いて列挙。n >= 1、根は頂点 0。
// 子の順序は区別せず、頂点番号は深さ優先順。1件あたり O(n)、追加領域 O(n)。
inline Generator<UnWeightedUnDirectedGraph> rooted_trees(Int n) {
    for (auto &levels : detail::itertools::rooted_tree_levels(n))
        co_yield detail::itertools::tree_from_levels(levels);
}

// n 頂点の根なし木を同型を除いて列挙。n >= 1、頂点番号は 0..<n。
// rooted_trees と同じ根付き木列挙から重心で代表を選ぶ。ハッシュ衝突なし。
// 根付き木の同型類数を R(n) として全体 O(n * R(n))、追加領域 O(n)。
inline Generator<UnWeightedUnDirectedGraph> unlabeled_trees(Int n) {
    for (auto &levels : detail::itertools::rooted_tree_levels(n))
        if (detail::itertools::unrooted_representative(levels))
            co_yield detail::itertools::tree_from_levels(levels);
}
}
