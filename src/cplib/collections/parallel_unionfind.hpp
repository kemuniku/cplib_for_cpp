#pragma once
#include <cplib/common.hpp>

namespace cplib {
// 4冪長の区間ごとのUnion-Find。構築O(N log N)。
class ParallelUnionFind {
    std::vector<std::int32_t> data;
    std::vector<Int> offsets;
    Int n = 0, components = 0;

    Int find(Int x) {
        Int result = x;
        while (data[result] >= 0)
            result = data[result];
        while (x != result) {
            Int next = data[x];
            data[x] = result;
            x = next;
        }
        return result;
    }

    template <class Callback> void uniteBlock(Int p, Int dis, Int layer, const Callback &onMerge) {
        Int base = offsets[layer], x = find(base + p), y = find(base + p + dis);
        if (x == y)
            return;
        if (data[x] > data[y])
            std::swap(x, y);
        if (layer == 0) {
            if constexpr (!std::is_same_v<Callback, std::nullptr_t>) {
                if constexpr (requires { bool(onMerge); }) {
                    if (onMerge)
                        onMerge(x, y);
                } else
                    onMerge(x, y);
            }
        }
        data[x] += data[y];
        data[y] = x;
        if (layer == 0) {
            --components;
            return;
        }
        Int step = Int(1) << (2 * layer - 2);
        for (Int i = 0; i < 4; ++i)
            uniteBlock(p + i * step, dis, layer - 1, onMerge);
    }

public:
    explicit ParallelUnionFind(Int N) : n(N), components(N) {
        assert(N >= 0);
        Int total = 0, width = 1;
        while (width <= N) {
            offsets.push_back(total);
            Int size = N - width + 1;
            assert(size <= std::numeric_limits<std::int32_t>::max() - total);
            total += size;
            if (width > N / 4)
                break;
            width *= 4;
        }
        data.assign(total, -1);
    }

    Int root(Int x) {
        assert(0 <= x && x < n);
        return find(x);
    }

    bool issame(Int x, Int y) {
        return root(x) == root(y);
    }

    Int siz(Int x) {
        return -Int(data[root(x)]);
    }

    Int count() const {
        return components;
    }

    std::vector<Int> roots() const {
        std::vector<Int> result;
        result.reserve(components);
        for (Int x = 0; x < n; ++x)
            if (data[x] < 0)
                result.push_back(x);
        return result;
    }

    // Q回の区間結合の合計O(N log N α(N)+Q log N)。結合直前にcallback(残る根,消える根)。
    // 各i∈[0,length)についてa+iとb+iを結合し、実際の結合回数を返します。
    // onMerge(x,y)は結合直前に呼び、結合後はxが代表になります。
    // コールバックからこの構造への結合操作は行わないでください。
    template <class Callback = std::nullptr_t>
    Int unite(Int a, Int b, Int length, Callback onMerge = nullptr) {
        assert(0 <= a && a <= n && 0 <= b && b <= n && 0 <= length && length <= n - a &&
               length <= n - b);
        if (a == b || length == 0)
            return 0;
        Int la = std::min(a, b), dis = std::max(a, b) - la, layer = 0, width = 1;
        while (width <= (length - 1) / 4) {
            ++layer;
            width *= 4;
        }
        Int before = components, remaining = length;
        while (remaining > 0) {
            remaining = std::max<Int>(0, remaining - width);
            uniteBlock(la + remaining, dis, layer, onMerge);
        }
        return before - components;
    }

    template <class Callback = std::nullptr_t>
        requires(!std::integral<Callback>)
    bool unite(Int a, Int b, Callback onMerge = nullptr) {
        return unite(a, b, 1, onMerge) != 0;
    }

    ParallelUnionFind copy() const {
        return *this;
    }
};

inline auto initParallelUnionFind(Int n) {
    return ParallelUnionFind(n);
}

inline Int root(ParallelUnionFind &p, Int x) {
    return p.root(x);
}

inline bool issame(ParallelUnionFind &p, Int x, Int y) {
    return p.issame(x, y);
}

inline Int siz(ParallelUnionFind &p, Int x) {
    return p.siz(x);
}

inline Int count(const ParallelUnionFind &p) {
    return p.count();
}

inline auto roots(const ParallelUnionFind &p) {
    return p.roots();
}

inline auto copy(const ParallelUnionFind &p) {
    return p.copy();
}

template <class... Args> auto unite(ParallelUnionFind &p, Args &&...args) {
    return p.unite(std::forward<Args>(args)...);
}
}
