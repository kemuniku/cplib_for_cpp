#pragma once
#include <cplib/geometry/base.hpp>

namespace cplib {
namespace detail::half_plane_count {
using Wide = __int128;
using UWide = unsigned __int128;

inline UWide magnitude(Wide value) {
    return value < 0 ? UWide(0) - UWide(value) : UWide(value);
}

// 各因子の絶対値は 2^64 未満。積同士を減算せずに符号付きで比較する。
inline int compare_products(Wide a, Wide b, Wide c, Wide d) {
    const UWide lhs = magnitude(a) * magnitude(b), rhs = magnitude(c) * magnitude(d);
    const bool lhs_negative = lhs != 0 && ((a < 0) != (b < 0));
    const bool rhs_negative = rhs != 0 && ((c < 0) != (d < 0));
    if (lhs_negative != rhs_negative)
        return lhs_negative ? -1 : 1;
    const int comparison = (lhs > rhs) - (lhs < rhs);
    return lhs_negative ? -comparison : comparison;
}

template <class T> bool xy_less(Point<T> lhs, Point<T> rhs) {
    return lhs.x < rhs.x || (lhs.x == rhs.x && lhs.y < rhs.y);
}

template <class T> struct Box {
    T xmin, xmax, ymin, ymax;

    explicit Box(std::span<const Point<T>> points)
        : xmin(points.front().x), xmax(xmin), ymin(points.front().y), ymax(ymin) {
        for (const auto &p : points) {
            xmin = std::min(xmin, p.x);
            xmax = std::max(xmax, p.x);
            ymin = std::min(ymin, p.y);
            ymax = std::max(ymax, p.y);
        }
    }

    bool narrow_projection(Wide a, Wide b) const {
        constexpr UWide limit = std::numeric_limits<Int>::max();
        const UWide aa = magnitude(a), bb = magnitude(b);
        if (aa > limit || bb > limit)
            return false;
        const UWide x = aa * std::max(magnitude(xmin), magnitude(xmax));
        const UWide y = bb * std::max(magnitude(ymin), magnitude(ymax));
        return x <= limit && y <= limit - x;
    }
};

enum class Relation { outside, inside, mixed };

template <class Number> struct Query {
    Number a, b, c;

    template <class T> bool contains(T x, T y) const {
        if constexpr (sizeof(Number) <= sizeof(Int)) {
            return a * x + b * y <= c;
        } else {
            const Wide lhs = a * x, rhs = b * y;
            Wide projection;
            if (__builtin_add_overflow(lhs, rhs, &projection))
                return lhs < 0;
            return projection <= c;
        }
    }

    template <class T> Relation classify(const Box<T> &box) const {
        if (contains(a < 0 ? box.xmin : box.xmax, b < 0 ? box.ymin : box.ymax))
            return Relation::inside;
        if (!contains(a < 0 ? box.xmax : box.xmin, b < 0 ? box.ymax : box.ymin))
            return Relation::outside;
        return Relation::mixed;
    }

    template <class T> std::size_t scan(std::span<const Point<T>> points) const {
        std::size_t result = 0;
        for (const auto &p : points)
            result += contains(p.x, p.y);
        return result;
    }
};

template <class T>
concept Coordinate = std::is_integral_v<T> && std::is_signed_v<T> && sizeof(T) <= sizeof(Int);
}

// 2直線の四分割を再帰する静的な partition tree。頂点の重複と境界を含めて数える。
// 部分木の凸包の上下鎖を保持する。空間 O(N log N)、1クエリ O(N^log_4(3))。
// 構築は凸包ソートと有理数探索で O(N log^2 N + W N log N)（選択の平均計算量、Wは座標ビット幅）。
template <detail::half_plane_count::Coordinate T = Int> class PartitionTree {
    using Box = detail::half_plane_count::Box<T>;
    using Relation = detail::half_plane_count::Relation;
    using Wide = detail::half_plane_count::Wide;
    static constexpr std::size_t leaf_size = 32;

    struct Node {
        Box box;
        std::size_t first, last;
        std::size_t lower = 0, upper = 0, end = 0;
        std::array<std::size_t, 4> children{};
    };
    struct Slope {
        UInt numerator, denominator;
    };
    using Iterator = typename std::vector<Point<T>>::iterator;

    std::vector<Point<T>> points, chains;
    std::vector<Node> nodes;

    static Slope advance_slope(Slope base, Slope step, UInt times) {
        return {base.numerator + times * step.numerator, base.denominator + times * step.denominator};
    }

    // 整数座標の二点が作る既約な傾きは、分子・分母とも UInt に収まる。
    // Stern–Brocot 木の同方向への連続移動を指数探索でまとめる。
    template <class Compare> static void find_slope(Compare compare) {
        constexpr UInt limit = std::numeric_limits<UInt>::max();
        std::array<Slope, 2> bounds{{{0, 1}, {1, 0}}};
        while (true) {
            const int side = compare(advance_slope(bounds[0], bounds[1], 1));
            if (side == 0)
                return;
            const auto index = side < 0 ? 0 : 1;
            const auto base = bounds[index], step = bounds[1 - index];
            UInt max_steps = limit;
            if (step.numerator != 0)
                max_steps = (limit - base.numerator) / step.numerator;
            if (step.denominator != 0)
                max_steps = std::min(max_steps, (limit - base.denominator) / step.denominator);
            UInt same = 1, other = 1;
            while (true) {
                other = same + std::min(same, max_steps - same);
                const int comparison = compare(advance_slope(base, step, other));
                if (comparison == 0)
                    return;
                if (comparison != side)
                    break;
                // 隣接する分数の間には、分子・分母が mediant 未満の分数は存在しない。
                assert(other < max_steps);
                same = other;
            }
            while (other - same > 1) {
                const UInt middle = same + (other - same) / 2;
                const int comparison = compare(advance_slope(base, step, middle));
                if (comparison == 0)
                    return;
                if (comparison == side)
                    same = middle;
                else
                    other = middle;
            }
            bounds[index] = advance_slope(base, step, same);
            bounds[1 - index] = advance_slope(base, step, other);
        }
    }

    static void split(Iterator first, Iterator lower, Iterator middle, Iterator upper, Iterator last, const Box &box) {
        std::nth_element(first, middle, last, detail::half_plane_count::xy_less<T>);
        auto select = [&](auto less) {
            std::nth_element(first, lower, middle, less);
            std::nth_element(middle, upper, last, less);
            return int(less(*upper, *lower)) - int(less(*lower, *upper));
        };
        auto compare = [&](Wide a, UInt b) {
            if (box.narrow_projection(a, b)) {
                const Int x = static_cast<Int>(a), y = static_cast<Int>(b);
                return select([&](Point<T> p, Point<T> q) { return x * p.x + y * p.y < x * q.x + y * q.y; });
            }
            return select([&](Point<T> p, Point<T> q) {
                return detail::half_plane_count::compare_products(a, Wide(p.x) - q.x, -Wide(b), Wide(p.y) - q.y) < 0;
            });
        };
        // 左右の射影中央値の差は、xの係数を増やすと単調非増加。
        const int direction = compare(0, 1);
        if (direction == 0)
            return;
        // 中央値を与える二点を結ぶ方向を数回試すと、通常は有理数探索に入らずに済む。
        for (int attempt = 0; attempt < 6; ++attempt) {
            const Wide dx = Wide(upper->x) - lower->x, dy = Wide(lower->y) - upper->y;
            if (dx == 0 || (dy > 0) - (dy < 0) != direction)
                break;
            if (compare(dy, static_cast<UInt>(dx)) == 0)
                return;
        }
        select([&](Point<T> p, Point<T> q) {
            return p.x == q.x ? p.y < q.y : (direction > 0 ? p.x < q.x : p.x > q.x);
        });
        // 端で中央値が同じxなら二つの子はその縦線上だけにある。
        // それらのy範囲は重ならない（端点を除く）ので、どの直線も両方の子を分断しない。
        if (lower->x == upper->x)
            return;
        find_slope([&](Slope slope) {
            return -direction * compare(direction * Wide(slope.numerator), slope.denominator);
        });
    }

    void append_chain(std::span<const Point<T>> sorted, bool lower) {
        const auto start = chains.size();
        for (const auto &p : sorted) {
            while (chains.size() - start >= 2) {
                const auto a = chains[chains.size() - 2], b = chains.back();
                const int turn = detail::half_plane_count::compare_products(
                    Wide(b.x) - a.x, Wide(p.y) - a.y, Wide(b.y) - a.y, Wide(p.x) - a.x);
                if (lower ? turn > 0 : turn < 0)
                    break;
                chains.pop_back();
            }
            chains.push_back(p);
        }
    }

    std::size_t build(std::size_t first, std::size_t last) {
        const auto group = std::span<const Point<T>>(points).subspan(first, last - first);
        const auto index = nodes.size();
        nodes.push_back({Box(group), first, last});
        if (last - first <= leaf_size)
            return index;

        {
            std::vector<Point<T>> sorted(group.begin(), group.end());
            std::sort(sorted.begin(), sorted.end(), detail::half_plane_count::xy_less<T>);
            sorted.erase(std::unique(sorted.begin(), sorted.end(), [](Point<T> p, Point<T> q) {
                return p.x == q.x && p.y == q.y;
            }), sorted.end());
            nodes[index].lower = chains.size();
            append_chain(sorted, true);
            nodes[index].upper = chains.size();
            append_chain(sorted, false);
            nodes[index].end = chains.size();
        }

        const auto middle = first + (last - first) / 2;
        const auto lower = first + (middle - first) / 2, upper = middle + (last - middle) / 2;
        const auto begin = points.begin();
        split(begin + first, begin + lower, begin + middle, begin + upper, begin + last, nodes[index].box);
        const std::array cuts{first, lower, middle, upper, last};
        for (std::size_t i = 0; i < 4; ++i)
            nodes[index].children[i] = build(cuts[i], cuts[i + 1]);
        return index;
    }

    template <class Number>
    int edge_projection(const Point<T> &p, const Point<T> &q, const detail::half_plane_count::Query<Number> &query) const {
        if constexpr (sizeof(Number) <= sizeof(Int)) {
            const auto lhs = query.a * p.x + query.b * p.y, rhs = query.a * q.x + query.b * q.y;
            return (lhs > rhs) - (lhs < rhs);
        }
        return detail::half_plane_count::compare_products(
            query.a, Wide(p.x) - q.x, -Wide(query.b), Wide(p.y) - q.y);
    }

    template <class Number>
    Relation classify_hull(const Node &node, const detail::half_plane_count::Query<Number> &query) const {
        const auto first = chains[node.lower];
        const bool inside = query.contains(first.x, first.y);
        if (node.end - node.lower <= 12) {
            for (auto i = node.lower + 1; i < node.end; ++i)
                if (query.contains(chains[i].x, chains[i].y) != inside)
                    return Relation::mixed;
        } else {
            // 全点が同じ側かどうかは、先頭点と反対側の極値を一つ調べれば分かる。
            const bool use_upper = inside == (query.b > 0);
            auto left = use_upper ? node.upper : node.lower;
            auto right = (use_upper ? node.end : node.upper) - 1;
            while (left < right) {
                const auto middle = left + (right - left) / 2;
                const int order = edge_projection(chains[middle], chains[middle + 1], query);
                if (inside ? order < 0 : order > 0)
                    left = middle + 1;
                else
                    right = middle;
            }
            if (query.contains(chains[left].x, chains[left].y) != inside)
                return Relation::mixed;
        }
        return inside ? Relation::inside : Relation::outside;
    }

    template <class Number>
    std::size_t count(std::size_t index, const detail::half_plane_count::Query<Number> &query) const {
        const auto &node = nodes[index];
        auto relation = query.classify(node.box);
        if (relation == Relation::mixed && node.end != 0 && query.a != 0 && query.b != 0)
            relation = classify_hull(node, query);
        if (relation == Relation::inside)
            return node.last - node.first;
        if (relation == Relation::outside)
            return 0;
        if (node.end == 0)
            return query.scan(std::span<const Point<T>>(points).subspan(node.first, node.last - node.first));
        std::size_t result = 0;
        for (auto child : node.children)
            result += count(child, query);
        return result;
    }

public:
    explicit PartitionTree(std::span<const Point<T>> input) : points(input.begin(), input.end()) {
        if (!points.empty())
            build(0, points.size());
    }

    std::size_t count(Int a, Int b, __int128 c) const {
        if (nodes.empty())
            return 0;
        if (nodes.front().box.narrow_projection(a, b)) {
            if (c < std::numeric_limits<Int>::min())
                return 0;
            if (c > std::numeric_limits<Int>::max())
                return points.size();
            return count(0, detail::half_plane_count::Query<Int>{a, b, static_cast<Int>(c)});
        }
        return count(0, detail::half_plane_count::Query<Wide>{a, b, c});
    }
};
}
