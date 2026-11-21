#pragma once
#include <cplib/geometry/half_plane_intersection.hpp>
#include <iterator>

namespace cplib {

// 半平面の共通部分として表す有界凸多角形。空・点・線分も扱う。
// T が符号付き整数なら制約は ax+by+c<=0 の係数、そうでなければ有向直線の左閉半平面。
// 演算は T のまま行う。係数生成、半平面交差、各問い合わせの中間値が T に収まること。
template <class T> class HalfPlaneConvexPolygon {
    static constexpr bool integer = detail::hpi::IntegerCoefficient<T>;
    static_assert(integer || cplib::isFloatingPoint<T> || detail::hpi::IsFraction<T>::value);

public:
    using constraint_type = std::conditional_t<integer, std::array<T, 3>, Line<T>>;
    using coordinate_type = std::conditional_t<integer, Fraction<T>, T>;
    using value_type = Point<coordinate_type>;

private:
    // 面では反時計回りの辺。点・線分では退化後も交差を続けられるよう定義制約を保持する。
    std::vector<constraint_type> constraints_;
    Int dimension_ = -1;

    // 点・線分の端点 (x/w,y/w)。面の頂点は隣接辺から求める。
    using HomogeneousPoint = std::array<T, 3>;
    std::array<HomogeneousPoint, 2> endpoints_{};

    // 有界な点・線分では下側・上側の包絡線が両方存在する。
    template <class Envelope, class Bound>
    static HomogeneousPoint endpoint(const Envelope &envelope, const Bound &x, bool upper) {
        using namespace detail::hpi;
        if constexpr (integer) {
            const auto &line = functionAt(envelope, x);
            const T y = line.p * x.num + line.q * x.den;
            return {T(x.num * line.d), upper ? T(-y) : y, T(x.den * line.d)};
        } else {
            const T y = valueAt(envelope, x.value);
            return {x.value, upper ? -y : y, T(1)};
        }
    }

    static Point<T> edgeDirection(const constraint_type &h) {
        if constexpr (integer)
            return {T(-h[1]), h[0]};
        else
            return h.t - h.s;
    }

    static HomogeneousPoint intersectLines(const constraint_type &a, const constraint_type &b) {
        if constexpr (integer) {
            T x = a[1] * b[2] - a[2] * b[1];
            T y = a[2] * b[0] - a[0] * b[2];
            T w = a[0] * b[1] - a[1] * b[0];
            return w < T(0) ? HomogeneousPoint{T(-x), T(-y), T(-w)} : HomogeneousPoint{x, y, w};
        } else {
            using detail::hpi::checked;
            const auto u = edgeDirection(a), v = edgeDirection(b);
            // 垂直線は x を直接使う。傾きで解き、長さの積によるアンダーフローを避ける。
            if (u.x == T(0))
                return {a.s.x, checked(b.s.y + checked(v.y / v.x) * (a.s.x - b.s.x)), T(1)};
            const T slope = checked(u.y / u.x);
            T x = b.s.x;
            if (v.x != T(0)) {
                const T otherSlope = checked(v.y / v.x);
                const T yAtOtherStart = checked(a.s.y + slope * (b.s.x - a.s.x));
                x = checked(b.s.x + (b.s.y - yAtOtherStart) / checked(slope - otherSlope));
            }
            return {x, checked(a.s.y + slope * (x - a.s.x)), T(1)};
        }
    }

    template <class U> static constexpr bool pointCoordinate =
        std::is_same_v<U, T> || (integer && std::is_same_v<U, coordinate_type>);

    template <class U> Int location(const Point<U> &point) const {
        constexpr bool rationalQuery = integer && !std::is_same_v<U, T>;
        if constexpr (!integer) {
            detail::hpi::checked(point.x);
            detail::hpi::checked(point.y);
        } else if constexpr (rationalQuery) {
            if (point.x.den <= T(0) || point.y.den <= T(0))
                throw std::invalid_argument("包含判定には正の分母の有限分数が必要です");
        }
        if (empty())
            return -1;
        bool boundary = dimension_ < 2;
        for (const auto &h : constraints_) {
            const T side = [&]() -> T {
                if constexpr (rationalQuery)
                    return T(detail::hpi::compareRatios<T>(h[0] * point.x.num, point.x.den,
                                                         -(h[1] * point.y.num + h[2] * point.y.den), point.y.den));
                else if constexpr (integer)
                    return h[0] * point.x + h[1] * point.y + h[2];
                else
                    return detail::hpi::checked(-cross(edgeDirection(h), point - h.s));
            }();
            if (side > T(0))
                return -1;
            boundary |= side == T(0);
        }
        return boundary ? 0 : 1;
    }

public:
    HalfPlaneConvexPolygon() = default;

    // 入力は任意順。冗長な制約を除き、正の長さの境界辺を周回順に保持する。
    // 計算量・数値型の前提は half_plane_intersection と同じ。非有界なら domain_error。
    explicit HalfPlaneConvexPolygon(std::span<const constraint_type> constraints) {
        const auto boundary = [&] {
            if constexpr (integer)
                return detail::hpi::integerIntersectionBoundary<T>(constraints);
            else
                return detail::hpi::intersectionBoundary<T>(constraints);
        }();
        if (boundary.kind == hpiUnbounded)
            throw std::domain_error("非有界な半平面交差は凸多角形にできません");
        dimension_ = boundary.dimension;
        if (empty())
            return;
        if (dimension_ == 2) {
            boundary.forEdges([&](std::size_t i) { constraints_.push_back(constraints[i]); });
        } else {
            constraints_.assign(constraints.begin(), constraints.end());
            endpoints_[0] = endpoint(boundary.lower, boundary.left, false);
            if (dimension_ == 1) {
                const bool vertical = compare(boundary.left, boundary.right) == 0;
                endpoints_[1] = vertical ? endpoint(boundary.upper, boundary.right, true)
                                         : endpoint(boundary.lower, boundary.right, false);
            }
        }
    }

    HalfPlaneConvexPolygon(std::initializer_list<constraint_type> constraints)
        : HalfPlaneConvexPolygon(std::span<const constraint_type>(constraints.begin(), constraints.size())) {}

    bool empty() const { return dimension_ < 0; }
    Int dimension() const { return dimension_; }
    std::size_t size() const { return dimension_ < 2 ? dimension_ + 1 : constraints_.size(); }

    // i 番目の頂点から次の頂点への辺に対応する。空・点・線分では空。
    std::span<const constraint_type> half_planes() const {
        return dimension_ == 2 ? std::span<const constraint_type>(constraints_) : std::span<const constraint_type>{};
    }

    // 整数版は参照した頂点だけを Fraction に変換する。
    value_type operator[](std::size_t i) const {
        assert(i < size());
        const auto [x, y, w] =
            dimension_ < 2 ? endpoints_[i] : intersectLines(constraints_[i == 0 ? size() - 1 : i - 1], constraints_[i]);
        if constexpr (integer)
            return {coordinate_type(x, w), coordinate_type(y, w)};
        else
            return {x, y};
    }

    class const_iterator {
        const HalfPlaneConvexPolygon *polygon_ = nullptr;
        std::size_t index_ = 0;

        friend class HalfPlaneConvexPolygon;
        const_iterator(const HalfPlaneConvexPolygon *polygon, std::size_t index) : polygon_(polygon), index_(index) {}

    public:
        using iterator_concept = std::forward_iterator_tag;
        using iterator_category = std::input_iterator_tag;
        using value_type = HalfPlaneConvexPolygon::value_type;
        using difference_type = std::ptrdiff_t;

        const_iterator() = default;
        value_type operator*() const { return (*polygon_)[index_]; }
        const_iterator &operator++() {
            ++index_;
            return *this;
        }
        const_iterator operator++(int) {
            auto previous = *this;
            ++*this;
            return previous;
        }
        bool operator==(const const_iterator &) const = default;
    };

    const_iterator begin() const { return {this, 0}; }
    const_iterator end() const { return {this, size()}; }

    ConvexPolygon<coordinate_type> toConvexPolygon() const {
        return ConvexPolygon<coordinate_type>(begin(), end());
    }

    Polygon<coordinate_type> toPolygon() const { return Polygon<coordinate_type>(begin(), end()); }

    // 制約を直接評価する。整数点は O(N) 演算。整数版への有理点は正の分母が必要で、
    // 分数比較の除算を含む。点・線分の狭義内部は空。
    template <class U = T> requires pointCoordinate<U>
    bool contains(const Point<U> &point, bool strict = false) const {
        return location(point) >= (strict ? 1 : 0);
    }

    template <class U = T> requires pointCoordinate<U>
    bool on_edge(const Point<U> &point) const { return location(point) == 0; }

    // 空でないこと。面では O(log N) 回の座標演算・比較。
    // 整数版の線分では、端点での内積の比較に分数比較の除算を使う。
    Int support_index(Point<T> direction) const {
        assert(!empty());
        if (dimension_ == 0 || (direction.x == T(0) && direction.y == T(0)))
            return 0;
        if (dimension_ == 1) {
            const auto &[x0, y0, w0] = endpoints_[0];
            const auto &[x1, y1, w1] = endpoints_[1];
            const T first = direction.x * x0 + direction.y * y0;
            const T last = direction.x * x1 + direction.y * y1;
            if constexpr (integer)
                return detail::hpi::compareRatios(last, w1, first, w0) > 0 ? 1 : 0;
            else
                return last > first ? 1 : 0;
        }
        return detail::convexEdgeSupportIndex(Int(size()), [&](Int i) { return edgeDirection(constraints_[i]); },
                                              direction);
    }

    value_type support_point(Point<T> direction) const { return (*this)[support_index(direction)]; }

    HalfPlaneConvexPolygon cut(const constraint_type &halfPlane) const {
        if (empty())
            return {};
        auto constraints = constraints_;
        constraints.push_back(halfPlane);
        return HalfPlaneConvexPolygon(constraints);
    }

    HalfPlaneConvexPolygon intersection(const HalfPlaneConvexPolygon &other) const {
        if (empty() || other.empty())
            return {};
        auto constraints = constraints_;
        constraints.insert(constraints.end(), other.constraints_.begin(), other.constraints_.end());
        return HalfPlaneConvexPolygon(constraints);
    }
};

// 正規化済みの凸多角形から同じ領域を表す制約を作る。座標型の昇格は行わない。
// 計算量は辺数に対する半平面交差と同じ。
template <class T, class Id>
HalfPlaneConvexPolygon<T> half_plane_convex_polygon(const ConvexPolygon<T, Id> &polygon) {
    using Result = HalfPlaneConvexPolygon<T>;
    using Constraint = typename Result::constraint_type;
    if (polygon.empty())
        return {};
    std::vector<Constraint> constraints;
    constraints.reserve(std::max<std::size_t>(4, polygon.size()));
    auto append = [&](Point<T> s, Point<T> t) {
        Line<T> line{s, t};
        if constexpr (detail::hpi::IntegerCoefficient<T>)
            constraints.push_back(line.toCoefficients());
        else
            constraints.push_back(line);
    };
    auto boundX = [&](const T &left, const T &right) {
        append({left, T(1)}, {left, T(0)});
        append({right, T(0)}, {right, T(1)});
    };
    auto boundY = [&](const T &bottom, const T &top) {
        append({T(0), bottom}, {T(1), bottom});
        append({T(1), top}, {T(0), top});
    };

    if (polygon.size() == 1) {
        boundX(polygon[0].x, polygon[0].x);
        boundY(polygon[0].y, polygon[0].y);
    } else if (polygon.size() == 2) {
        const Point<T> s = polygon[0], t = polygon[1];
        append(s, t);
        append(t, s);
        if (s.x != t.x)
            boundX(s.x, t.x);
        else
            boundY(s.y, t.y);
    } else {
        for (std::size_t i = 0; i < polygon.size(); ++i)
            append(polygon[i], polygon[(i + 1) % polygon.size()]);
    }
    return Result(constraints);
}

template <class T, class U = T>
    requires requires(const HalfPlaneConvexPolygon<T> &polygon, const Point<U> &point) { polygon.contains(point); }
bool contains(const HalfPlaneConvexPolygon<T> &polygon, const Point<U> &point, bool strict = false) {
    return polygon.contains(point, strict);
}

template <class T, class U = T>
    requires requires(const HalfPlaneConvexPolygon<T> &polygon, const Point<U> &point) { polygon.on_edge(point); }
bool on_edge(const HalfPlaneConvexPolygon<T> &polygon, const Point<U> &point) {
    return polygon.on_edge(point);
}

template <class T> Int support_index(const HalfPlaneConvexPolygon<T> &polygon, Point<T> direction) {
    return polygon.support_index(direction);
}

template <class T> auto support_point(const HalfPlaneConvexPolygon<T> &polygon, Point<T> direction) {
    return polygon.support_point(direction);
}

template <class T> auto toConvexPolygon(const HalfPlaneConvexPolygon<T> &polygon) {
    return polygon.toConvexPolygon();
}

template <class T>
auto convex_cut(const HalfPlaneConvexPolygon<T> &polygon, const typename HalfPlaneConvexPolygon<T>::constraint_type &h) {
    return polygon.cut(h);
}

template <detail::hpi::IntegerCoefficient T>
auto convex_cut(const HalfPlaneConvexPolygon<T> &polygon, const Line<T> &line) {
    if (exact_equal(line.s, line.t))
        throw std::invalid_argument("半平面の端点は異なる必要があります");
    return polygon.cut(line.toCoefficients());
}

template <class T> auto convex_intersection(const HalfPlaneConvexPolygon<T> &a, const HalfPlaneConvexPolygon<T> &b) {
    return a.intersection(b);
}

} // namespace cplib
