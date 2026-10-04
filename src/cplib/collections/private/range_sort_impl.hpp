#pragma once
#include <cplib/utils/backwards_index.hpp>
#include <cplib/utils/itertools.hpp>
#include <unordered_set>
#include <functional>
#include <cplib/collections/private/sort_order.hpp>

namespace cplib {

namespace detail {
struct RangeSortEmpty {};

// 相異なる[0,keyLimit)内の整数キーによる区間ソートを扱う。ソート時にキーと値は一緒に移動する。
// Products=trueでは現在の並び順のモノイド積を保持する。mergeは左、右の順に結合し、eを単位元とする。
// Nは要素数、UはkeyLimit。mergeがO(1)なら構築とQ操作の合計はO((N+Q)(log(N+1)+log(U+1)))、空間O(N)。
// ソートの計算量は償却で、不要なノードは再利用する。空区間の積はe、ソートは何もしない。
template <class T, bool Products> class RangeSortEngine {
    struct Product {
        T ascending{}, descending{};
    };

    struct Node : std::conditional_t<Products, Product, RangeSortEmpty> {
        Int left = 0, right = 0, count = 0, key = 0, bit = 0;
    };

    std::vector<Node> nodes;
    std::vector<Int> freeNodes, roots, next, starts;
    std::vector<bool> reversed;
    std::vector<T> products;
    Int keyLimit, treeBase_ = 1;
    std::function<T(T, T)> merge;
    T e;
#ifndef NDEBUG
    std::unordered_set<Int> activeKeys;
#endif
    // 空きノードを再利用し、なければ確保します。償却O(1)。
    Int newNode() {
        if (!freeNodes.empty()) {
            Int n = freeNodes.back();
            freeNodes.pop_back();
            return n;
        }
        nodes.emplace_back();
        return nodes.size() - 1;
    }

    // ノードの値を解放して再利用候補にします。O(1)。
    void releaseNode(Int n) {
        nodes[n] = {};
        freeNodes.push_back(n);
    }

    // 子から要素数と最小キーを更新する。Products=trueでは両方向の積も更新する。mergeがO(1)ならO(1)。
    void pull(Int n) {
        Int a = nodes[n].left, b = nodes[n].right;
        nodes[n].key = nodes[a].key;
        nodes[n].count = nodes[a].count + nodes[b].count;
        if constexpr (Products) {
            nodes[n].ascending = merge(nodes[a].ascending, nodes[b].ascending);
            nodes[n].descending = merge(nodes[b].descending, nodes[a].descending);
        }
    }

    // 指定した分岐ビットと空でない子から内部ノードを生成します。償却O(1)。
    Int branch(Int a, Int b, Int bit) {
        Int n = newNode();
        nodes[n].left = a;
        nodes[n].right = b;
        nodes[n].bit = bit;
        pull(n);
        return n;
    }

    // キー1個の葉を生成します。償却O(1)。
    Int singleton(Int key, T value) {
        Int n = newNode();
        nodes[n].count = 1;
        nodes[n].key = key;
        nodes[n].bit = -1;
        if constexpr (Products)
            nodes[n].ascending = nodes[n].descending = value;
        return n;
    }

    // 小さいキーからk個と残りに破壊的分割します。O(log(U+1))。
    std::pair<Int, Int> split(Int node, Int k) {
        if (!k)
            return {0, node};
        if (k == nodes[node].count)
            return {node, 0};
        Int l = nodes[node].left, r = nodes[node].right, c = nodes[l].count;
        if (k == c) {
            releaseNode(node);
            return {l, r};
        }
        if (k < c) {
            auto [a, b] = split(l, k);
            nodes[node].left = b;
            pull(node);
            return {a, node};
        }
        auto [a, b] = split(r, k - c);
        nodes[node].right = a;
        pull(node);
        return {node, b};
    }

    // 分岐のない経路を省略したキーの木を破壊的に融合します。
    Int meld(Int a, Int b) {
        Int ab = nodes[a].bit, bb = nodes[b].bit;
        UInt diff = UInt(nodes[a].key) ^ UInt(nodes[b].key);
        Int bit = diff ? Int(std::bit_width(diff)) - 1 : -1;
        if (bit > std::max(ab, bb))
            return nodes[a].key < nodes[b].key ? branch(a, b, bit) : branch(b, a, bit);
        if (ab < bb)
            return meld(b, a);
        if (ab > bb) {
            if ((UInt(nodes[b].key) >> ab) & 1) {
                Int child = meld(nodes[a].right, b);
                nodes[a].right = child;
            } else {
                Int child = meld(nodes[a].left, b);
                nodes[a].left = child;
            }
        } else {
            assert(ab >= 0);
            Int l = meld(nodes[a].left, nodes[b].left), r = meld(nodes[a].right, nodes[b].right);
            nodes[a].left = l;
            nodes[a].right = r;
            releaseNode(b);
        }
        pull(a);
        return a;
    }

    // indexを含むブロックの始点を求めます。O(log(N+1))。
    Int blockStart(Int index) const {
        if (roots[index])
            return index;
        for (Int n = index + treeBase_; n > 1; n >>= 1)
            if ((n & 1) && starts[n - 1] != -1)
                return starts[n - 1];
        return -1;
    }

    // ブロック境界を更新し、値の変わらない祖先で打ち切ります。O(log(N+1))。
    void setStart(Int i, Int v) {
        Int n = i + treeBase_;
        starts[n] = v;
        while ((n >>= 1) > 0) {
            Int value = std::max(starts[n * 2], starts[n * 2 + 1]);
            if (starts[n] == value)
                break;
            starts[n] = value;
        }
    }

    void setProduct(Int i, T v) {
        if constexpr (Products) {
            Int n = i + treeBase_;
            products[n] = v;
            while ((n >>= 1) > 0)
                products[n] = merge(products[n * 2], products[n * 2 + 1]);
        }
    }

    // ブロック全体の積を外側のセグ木へ反映します。O(log(N+1))。
    void refresh(Int i) {
        if constexpr (Products) {
            auto &n = nodes[roots[i]];
            setProduct(i, reversed[i] ? n.descending : n.ascending);
        }
    }

    T product(Int l, Int r) const {
        T a = e, b = e;
        for (l += treeBase_, r += treeBase_; l < r; l >>= 1, r >>= 1) {
            if (l & 1)
                a = merge(a, products[l++]);
            if (r & 1)
                b = merge(products[--r], b);
        }
        return merge(a, b);
    }

    // 配列上のiの直前にブロック境界を作る。O(log(N+1)+log(U+1))。
    void cut(Int i) {
        if (i == len() || roots[i])
            return;
        Int start = blockStart(i), root = roots[start], k = i - start;
        auto [a, b] = split(root, reversed[start] ? nodes[root].count - k : k);
        if (reversed[start])
            std::swap(a, b);
        roots[start] = a;
        roots[i] = b;
        reversed[i] = reversed[start];
        next[i] = next[start];
        next[start] = i;
        setStart(i, i);
        refresh(start);
        refresh(i);
    }

    // 現在の位置iの葉の添字を返す。O(log(N+1)+log(U+1))。
    Int locate(Int i) const {
        assert(0 <= i && i < len());
        Int start = blockStart(i), node = roots[start], k = i - start;
        if (reversed[start])
            k = nodes[node].count - 1 - k;
        while (nodes[node].count > 1) {
            Int l = nodes[node].left;
            if (k < nodes[l].count)
                node = l;
            else {
                k -= nodes[l].count;
                node = nodes[node].right;
            }
        }
        return node;
    }

    // 指定キーの値と祖先の積を更新します。O(log(U+1))。
    void updateValue(Int node, Int key, T value) {
        if (nodes[node].count == 1) {
            nodes[node].ascending = nodes[node].descending = value;
            return;
        }
        updateValue((UInt(key) >> nodes[node].bit) & 1 ? nodes[node].right : nodes[node].left, key,
                    value);
        pull(node);
    }

    // キー昇順の順位[l,r)の積を、指定された方向で返します。O(log(U+1))。
    T nodeProduct(Int node, Int l, Int r, bool rev) const {
        if (!l && r == nodes[node].count)
            return rev ? nodes[node].descending : nodes[node].ascending;
        Int a = nodes[node].left, b = nodes[node].right, c = nodes[a].count;
        if (r <= c)
            return nodeProduct(a, l, r, rev);
        if (l >= c)
            return nodeProduct(b, l - c, r - c, rev);
        T x = nodeProduct(a, l, c, rev), y = nodeProduct(b, 0, r - c, rev);
        return rev ? merge(y, x) : merge(x, y);
    }

    // 1ブロック内の位置[l,r)の積を返します。O(log(U+1))。
    T blockProduct(Int start, Int l, Int r) const {
        return reversed[start] ? nodeProduct(roots[start], next[start] - r, next[start] - l, true)
                               : nodeProduct(roots[start], l - start, r - start, false);
    }

public:
    // 与えた並び順で構築する。Nは要素数、Uはlimit。O(N)時間・空間。
    // キーは[0,limit)で相異なること。Products=trueの場合はkeysとvaluesを同じ長さとする。
    // キーの重複はassert有効時に検査する（ハッシュ集合操作は期待O(1)）。
    template <class Op>
    RangeSortEngine(std::span<const Int> keys, std::span<const T> values, Int limit, Op op,
                    T identity)
        : roots(keys.size()), next(keys.size()), reversed(keys.size()), keyLimit(limit), merge(op),
          e(identity) {
        assert(limit >= 0);
        if constexpr (Products)
            assert(keys.size() == values.size());
        while (treeBase_ < len())
            treeBase_ *= 2;
        starts.assign(treeBase_ * 2, -1);
        if constexpr (Products)
            products.assign(treeBase_ * 2, e);
        nodes.reserve(std::max<Int>(1, len() * 2));
        nodes.emplace_back();
        for (Int i = 0; i < len(); ++i) {
            Int key = keys[i];
            assert(0 <= key && key < limit);
#ifndef NDEBUG
            assert(activeKeys.insert(key).second);
#endif
            roots[i] = singleton(key, Products ? values[i] : e);
            next[i] = i + 1;
            starts[treeBase_ + i] = i;
            if constexpr (Products)
                products[treeBase_ + i] = values[i];
        }
        for (Int i = treeBase_ - 1; i > 0; --i) {
            starts[i] = std::max(starts[i * 2], starts[i * 2 + 1]);
            if constexpr (Products)
                products[i] = merge(products[i * 2], products[i * 2 + 1]);
        }
    }

    // 要素数をO(1)で返す。
    Int len() const {
        return roots.size();
    }

    Int size() const {
        return len();
    }

    // 現在の位置iのキーを返す。O(log(N+1)+log(U+1))。
    Int key(Int i) const {
        return nodes[locate(i)].key;
    }

    // 現在の位置iの値を取得する。O(log(N+1)+log(U+1))。
    T value(Int i) const
        requires Products
    {
        return nodes[locate(i)].ascending;
    }

    // キーを保って値を変更する。O(log(N+1)+log(U+1))。
    void update(Int i, T value)
        requires Products
    {
        Int k = key(i), start = blockStart(i);
        updateValue(roots[start], k, value);
        refresh(start);
    }

    // 位置iのキーと値を変更する。O(log(N+1)+log(U+1))。他の要素とのキー重複は禁止。
    void update(Int i, Int k, T value) {
        assert(0 <= i && i < len() && 0 <= k && k < keyLimit);
#ifndef NDEBUG
        Int old = key(i);
        assert(k == old || !activeKeys.contains(k));
        activeKeys.erase(old);
        activeKeys.insert(k);
#endif
        cut(i);
        cut(i + 1);
        releaseNode(roots[i]);
        roots[i] = singleton(k, value);
        reversed[i] = false;
        refresh(i);
    }

    // 現在の並び順で[l,r)の積を返す。空区間は単位元。O(log(N+1)+log(U+1))。
    T get(Int l, Int r) const
        requires Products
    {
        assert(0 <= l && l <= r && r <= len());
        if (l == r)
            return e;
        Int a = blockStart(l), b = blockStart(r - 1);
        if (a == b)
            return blockProduct(a, l, r);
        T out = blockProduct(a, l, next[a]);
        if (next[a] < b)
            out = merge(out, product(next[a], b));
        return merge(out, blockProduct(b, b, r));
    }

    // 全体の積をO(1)で返す。
    T get_all() const
        requires Products
    {
        return products[1];
    }

    // [l,r)をキー順にソートする。償却O(log(N+1)+log(U+1))。空区間は変更しない。
    void sort(Int l, Int r, SortOrder order = Ascending) {
        assert(0 <= l && l <= r && r <= len());
        if (r - l <= 1)
            return;
        cut(l);
        cut(r);
        Int root = roots[l];
        for (Int start = next[l]; start < r; start = next[start]) {
            root = meld(root, roots[start]);
            roots[start] = 0;
            setProduct(start, e);
            setStart(start, -1);
        }
        roots[l] = root;
        next[l] = r;
        reversed[l] = order == Descending;
        refresh(l);
    }

    // 現在の配列順にキーを列挙する。全体O(N)。列挙中の変更は禁止。
    Generator<Int> items() const {
        std::vector<Int> stack;
        for (Int start = 0; start < len(); start = next[start]) {
            stack.push_back(roots[start]);
            while (!stack.empty()) {
                Int node = stack.back();
                stack.pop_back();
                if (nodes[node].count == 1)
                    co_yield nodes[node].key;
                else if (reversed[start]) {
                    stack.push_back(nodes[node].left);
                    stack.push_back(nodes[node].right);
                } else {
                    stack.push_back(nodes[node].right);
                    stack.push_back(nodes[node].left);
                }
            }
        }
    }
};
}
}
