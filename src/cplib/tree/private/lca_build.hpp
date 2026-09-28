#pragma once

#include <cstdint>
namespace cplib::detail::lca_native {
template<bool Ordered>
static std::int64_t build(const std::int32_t* __restrict parent,
        const std::int64_t* __restrict order, std::int64_t n, std::int64_t root,
        std::int32_t* __restrict size, std::uint32_t* __restrict data,
        std::uint8_t* __restrict branches, std::int32_t* __restrict head) {
    // 葉の区間と縦パスの情報を、配列同士が重ならない条件で構築する。
    for (std::int64_t i = n - 1; i > 0; --i) {
        if (i > 16) {
            const std::int64_t future = Ordered ? i - 16 : order[i - 16];
            __builtin_prefetch(size + parent[future], 1, 3);
        }
        const std::int64_t v = Ordered ? i : order[i];
        const std::int32_t count = size[v] == 0 ? 1 : size[v];
        size[v] = count;
        size[parent[v]] += count;
    }
    if (size[root] == 0) size[root] = 1;
    const std::uint32_t root_label = 1U << (31 - __builtin_clz(static_cast<unsigned>(size[root])));
    data[3 * root] = data[3 * root + 1] = root_label;
    data[3 * root + 2] = 0;
    branches[root] = 0;
    head[root_label] = -1;
    std::int64_t deepest = root;
    std::uint32_t max_depth = 0;
    for (std::int64_t i = 1; i < n; ++i) {
        // 後で参照する親の情報を先にキャッシュへ読み込む。
        if (i + 16 < n) {
            const std::int64_t future = Ordered ? i + 16 : order[i + 16];
            const std::int64_t p = parent[future];
            __builtin_prefetch(data + 3 * p, 0, 3);
            __builtin_prefetch(branches + p, 0, 3);
            __builtin_prefetch(size + p, 1, 3);
        }
        const std::int64_t v = Ordered ? i : order[i];
        const std::int64_t p = parent[v];
        const std::uint32_t end = size[p];
        const std::uint32_t begin = end - size[v];
        size[p] = begin;
        size[v] = end;
        const unsigned shift = 31 - __builtin_clz(begin ^ end);
        const std::uint32_t label = end & (~0U << shift);
        const std::uint32_t bit = label & -label;
        branches[v] = i < 64 ? i : branches[p];
        data[3 * v] = label;
        data[3 * v + 1] = data[3 * p + 1] | bit;
        const std::uint32_t depth = data[3 * p + 2] + 1;
        data[3 * v + 2] = depth;
        if (depth > max_depth) { max_depth = depth; deepest = v; }
        if (label != data[3 * p]) head[label] = p;
    }
    return deepest;
}
}
