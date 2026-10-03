#pragma once
#include <cplib/common.hpp>
#include <memory>
#include <bit>

namespace cplib {
// 16分木の追記専用プールを版間で共有する。二点更新の共通経路は一度だけコピー。
class PersistentUnionFind {
    using Node = std::array<std::int32_t, 16>;

    struct Pool {
        std::vector<Node> nodes{Node{}};

        std::int32_t copyNode(std::int32_t node, Int shift) {
            Node data = nodes[node];
            if (node == 0 && shift == 0)
                data.fill(-1);
            assert(nodes.size() < std::size_t(std::numeric_limits<std::int32_t>::max()));
            std::int32_t result = nodes.size();
            nodes.push_back(data);
            return result;
        }

        std::int32_t setOne(std::int32_t node, Int shift, Int index, std::int32_t value) {
            auto result = copyNode(node, shift);
            Int slot = (index >> shift) & 15;
            if (shift == 0)
                nodes[result][slot] = value;
            else {
                auto child = setOne(nodes[result][slot], shift - 4, index, value);
                nodes[result][slot] = child;
            }
            return result;
        }

        std::int32_t setTwo(std::int32_t node, Int shift, Int x, Int y, std::int32_t xv,
                            std::int32_t yv) {
            auto result = copyNode(node, shift);
            Int xs = (x >> shift) & 15, ys = (y >> shift) & 15;
            if (shift == 0) {
                nodes[result][xs] = xv;
                nodes[result][ys] = yv;
            } else if (xs == ys) {
                auto child = setTwo(nodes[result][xs], shift - 4, x, y, xv, yv);
                nodes[result][xs] = child;
            } else {
                auto xc = setOne(nodes[result][xs], shift - 4, x, xv);
                nodes[result][xs] = xc;
                auto yc = setOne(nodes[result][ys], shift - 4, y, yv);
                nodes[result][ys] = yc;
            }
            return result;
        }
    };

    Int size = 0, height = 0;
    std::int32_t node = 0;
    std::shared_ptr<Pool> pool;

    std::int32_t get(Int index) const {
        auto v = node;
        Int shift = height;
        while (v != 0 && shift > 0) {
            v = pool->nodes[v][(index >> shift) & 15];
            shift -= 4;
        }
        return v == 0 ? -1 : pool->nodes[v][index & 15];
    }

    std::pair<Int, std::int32_t> rootAndSize(Int x) const {
        assert(0 <= x && x < size);
        auto value = get(x);
        while (value >= 0) {
            x = value;
            value = get(x);
        }
        return {x, value};
    }

public:
    Int count = 0;

    explicit PersistentUnionFind(Int n) : size(n), pool(std::make_shared<Pool>()), count(n) {
        assert(n >= 0 && n <= std::numeric_limits<std::int32_t>::max());
        if (n > 1)
            height = ((std::bit_width(UInt(n - 1)) - 1) / 4) * 4;
    }

    Int root(Int x) const {
        return rootAndSize(x).first;
    }

    bool issame(Int x, Int y) const {
        return root(x) == root(y);
    }

    Int siz(Int x) const {
        return -Int(rootAndSize(x).second);
    }

    // Union by size。根探索O(log² N)、追加ノードO(log N)。以前の版は変更しない。
    PersistentUnionFind unite(Int x, Int y) const {
        auto a = rootAndSize(x), b = rootAndSize(y);
        auto result = *this;
        if (a.first == b.first)
            return result;
        if (a.second > b.second)
            std::swap(a, b);
        result.node = pool->setTwo(node, height, a.first, b.first, a.second + b.second,
                                   std::int32_t(a.first));
        --result.count;
        return result;
    }
};

inline auto initPersistentUnionFind(Int n) {
    return PersistentUnionFind(n);
}

inline Int root(const PersistentUnionFind &p, Int x) {
    return p.root(x);
}

inline bool issame(const PersistentUnionFind &p, Int x, Int y) {
    return p.issame(x, y);
}

inline Int siz(const PersistentUnionFind &p, Int x) {
    return p.siz(x);
}

inline auto unite(const PersistentUnionFind &p, Int x, Int y) {
    return p.unite(x, y);
}
}
