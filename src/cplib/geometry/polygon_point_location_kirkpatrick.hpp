#pragma once
#include <cplib/geometry/polygon.hpp>
#include <boost/multiprecision/cpp_int.hpp>
#include <algorithm>
#include <array>
#include <cassert>
#include <cstdint>
#include <deque>
#include <map>
#include <numeric>
#include <set>
#include <type_traits>
#include <utility>
#include <vector>

namespace cplib {

// 三角形分割と Kirkpatrick 階層によるオンライン点位置探索。
// 構築 O(N log N)、領域 O(N)、クエリ O(log N)。初期分割は単調分割で行う。
// 整数座標専用。中間計算のオーバーフローを避けるため |x|, |y| <= 10^9 を想定する。
template <class T> class KirkpatrickPolygonPointLocator {
    static_assert(std::is_integral_v<T> && sizeof(T) <= 8);
    using Wide = __int128_t;
    using P = Point<Wide>;
    struct Triangle {
        std::array<int, 3> v;
        std::vector<int> children;
        bool inside = false;
    };
    std::vector<P> points_;
    std::vector<Triangle> nodes_;
    std::vector<int> top_;
    std::vector<std::pair<int, int>> boundary_;
    std::vector<P> original_vertices_;
    std::array<int, 4> frame_{};
    P bottom_left_{}, top_right_{};

    static Wide cross(P a, P b, P c) {
        return (b.x - a.x) * (c.y - a.y) - (b.y - a.y) * (c.x - a.x);
    }
    static Wide signed_area(const std::vector<int> &ring, const std::vector<P> &points) {
        Wide area = 0;
        for (int i = 0; i < int(ring.size()); ++i) {
            P a = points[ring[i]], b = points[ring[(i + 1) % ring.size()]];
            area += a.x * b.y - a.y * b.x;
        }
        return area;
    }
    static bool on_segment(P a, P b, P p) {
        return cross(a, b, p) == 0 && std::min(a.x, b.x) <= p.x && p.x <= std::max(a.x, b.x) &&
               std::min(a.y, b.y) <= p.y && p.y <= std::max(a.y, b.y);
    }
    bool in_triangle(const std::array<int, 3> &tri, P p) const {
        return cross(points_[tri[0]], points_[tri[1]], p) >= 0 &&
               cross(points_[tri[1]], points_[tri[2]], p) >= 0 &&
               cross(points_[tri[2]], points_[tri[0]], p) >= 0;
    }
    int add_triangle(std::array<int, 3> v, std::vector<int> children = {}, bool inside = false) {
        Wide c = cross(points_[v[0]], points_[v[1]], points_[v[2]]);
        assert(c != 0);
        if (c < 0)
            std::swap(v[1], v[2]);
        nodes_.push_back({v, std::move(children), inside});
        return int(nodes_.size()) - 1;
    }
    // 階層更新時に生じる次数 12 以下の局所領域だけを分割する。
    std::vector<std::array<int, 3>> triangulate_quadratic(std::vector<int> ring) const {
        if (signed_area(ring, points_) < 0)
            std::reverse(ring.begin(), ring.end());
        const int m = int(ring.size());
        std::vector<int> previous(m), following(m), bad(m), generation(m);
        std::vector<char> alive(m, true);
        for (int i = 0; i < m; ++i) {
            previous[i] = (i + m - 1) % m;
            following[i] = (i + 1) % m;
        }
        auto convex = [&](int i) {
            int a = ring[previous[i]], b = ring[i], c = ring[following[i]];
            return a != b && b != c && c != a && cross(points_[a], points_[b], points_[c]) > 0;
        };
        auto blocks = [&](int i, int j) {
            if (j == previous[i] || j == i || j == following[i])
                return false;
            P a = points_[ring[previous[i]]], b = points_[ring[i]], c = points_[ring[following[i]]];
            P p = points_[ring[j]];
            if (p == a || p == b || p == c)
                return false;
            return on_segment(a, c, p) ||
                   (cross(a, b, p) > 0 && cross(b, c, p) > 0 && cross(c, a, p) > 0);
        };
        std::deque<std::pair<int, int>> ears;
        auto recompute = [&](int i) {
            ++generation[i];
            bad[i] = 0;
            if (!convex(i))
                return;
            for (int j = 0; j < m; ++j)
                if (alive[j] && blocks(i, j))
                    ++bad[i];
            if (bad[i] == 0)
                ears.push_back({i, generation[i]});
        };
        for (int i = 0; i < m; ++i)
            recompute(i);
        std::vector<std::array<int, 3>> result;
        result.reserve(m - 2);
        int left = m, cursor = 0;
        while (left > 3) {
            while (!ears.empty() && (!alive[ears.front().first] ||
                   ears.front().second != generation[ears.front().first]))
                ears.pop_front();
            assert(!ears.empty());
            cursor = ears.front().first;
            ears.pop_front();
            int a = previous[cursor], c = following[cursor];
            result.push_back({ring[a], ring[cursor], ring[c]});
            alive[cursor] = false;
            following[a] = c;
            previous[c] = a;
            --left;
            for (int i = 0; i < m; ++i)
                if (alive[i] && i != a && i != c && convex(i) && bad[i] > 0 &&
                    blocks(i, cursor) && --bad[i] == 0)
                    ears.push_back({i, generation[i]});
            recompute(a);
            recompute(c);
        }
        while (!alive[cursor])
            cursor = following[cursor];
        result.push_back({ring[cursor], ring[following[cursor]], ring[previous[cursor]]});
        return result;
    }
    static std::pair<int, int> edge_key(int a, int b) {
        return {std::min(a, b), std::max(a, b)};
    }
    void triangulate_monotone_face(const std::vector<int> &face,
                                   const std::vector<Wide> &y,
                                   const std::vector<int> &ring,
                                   std::vector<std::array<int, 3>> &result) const {
        const int m = int(face.size());
        assert(m >= 3);
        auto emit = [&](int a, int b, int c) {
            assert(cross(points_[ring[a]], points_[ring[b]], points_[ring[c]]) != 0);
            result.push_back({ring[a], ring[b], ring[c]});
        };
        if (m == 3) {
            emit(face[0], face[1], face[2]);
            return;
        }
        int top = 0, bottom = 0;
        for (int i = 1; i < m; ++i) {
            if (y[face[i]] > y[face[top]]) top = i;
            if (y[face[i]] < y[face[bottom]]) bottom = i;
        }
        int at = top;
        while (at != bottom) {
            int next = (at + 1) % m;
            assert(y[face[next]] < y[face[at]]);
            at = next;
        }
        while (at != top) {
            int next = (at + 1) % m;
            assert(y[face[next]] > y[face[at]]);
            at = next;
        }
        std::vector<int> priority;
        priority.reserve(m);
        std::vector<int> side(m);
        priority.push_back(top);
        int left = (top + 1) % m, right = (top + m - 1) % m;
        for (int i = 1; i < m - 1; ++i) {
            if (left == bottom ||
                (right != bottom && y[face[right]] > y[face[left]])) {
                priority.push_back(right);
                side[right] = -1;
                right = (right + m - 1) % m;
            } else {
                priority.push_back(left);
                side[left] = 1;
                left = (left + 1) % m;
            }
        }
        priority.push_back(bottom);
        std::vector<int> stack{priority[0], priority[1]};
        for (int i = 2; i < m - 1; ++i) {
            int v = priority[i];
            if (side[v] != side[stack.back()]) {
                for (int j = 0; j + 1 < int(stack.size()); ++j) {
                    if (side[v] == 1) emit(face[stack[j + 1]], face[stack[j]], face[v]);
                    else emit(face[stack[j]], face[stack[j + 1]], face[v]);
                }
                stack = {priority[i - 1], v};
            } else {
                int last = stack.back();
                stack.pop_back();
                while (!stack.empty()) {
                    int a = v;
                    int b = side[v] == 1 ? stack.back() : last;
                    int c = side[v] == 1 ? last : stack.back();
                    if (cross(points_[ring[face[a]]], points_[ring[face[b]]],
                              points_[ring[face[c]]]) <= 0)
                        break;
                    emit(face[a], face[b], face[c]);
                    last = stack.back();
                    stack.pop_back();
                }
                stack.push_back(last);
                stack.push_back(v);
            }
        }
        int v = priority.back();
        for (int j = 0; j + 1 < int(stack.size()); ++j) {
            if (side[stack[j + 1]] == 1) emit(face[stack[j]], face[stack[j + 1]], face[v]);
            else emit(face[stack[j + 1]], face[stack[j]], face[v]);
        }
    }
    // 単純多角形を単調分割し、各面を線形時間で三角形分割する。
    // shear > 全頂点の x 座標差とし、y' = shear*y + x で頂点の高さを一意にする。
    std::vector<std::array<int, 3>> triangulate_monotone(std::vector<int> ring,
                                                           Wide shear) const {
        using Big = boost::multiprecision::int256_t;
        if (signed_area(ring, points_) < 0)
            std::reverse(ring.begin(), ring.end());
        const int m = int(ring.size());
        assert(m >= 3);
        std::vector<Wide> x(m), y(m);
        for (int i = 0; i < m; ++i) {
            P p = points_[ring[i]];
            x[i] = p.x;
            y[i] = shear * p.y + p.x;
        }
        struct SweepEdge {
            Wide xhigh, yhigh, xlow, dy, dxdown;
        };
        std::vector<SweepEdge> sweep_edges(m);
        for (int i = 0; i < m; ++i) {
            int j = (i + 1) % m;
            int high = y[i] > y[j] ? i : j, low = high == i ? j : i;
            assert(y[high] != y[low]);
            sweep_edges[i] = {x[high], y[high], x[low], y[high] - y[low], x[low] - x[high]};
        }
        Wide sweep_y = 0;
        struct SweepCompare {
            const std::vector<SweepEdge> *edges;
            const std::vector<Wide> *x;
            const Wide *sweep_y;
            int m;
            bool operator()(int a, int b) const {
                if (a == b) return false;
                auto value = [&](int id) {
                    if (id >= m)
                        return std::pair<Big, Big>{Big((*x)[id - m]), Big(1)};
                    const SweepEdge &e = (*edges)[id];
                    Wide num = e.xhigh * e.dy + e.dxdown * (e.yhigh - *sweep_y);
                    return std::pair<Big, Big>{Big(num), Big(e.dy)};
                };
                auto [an, ad] = value(a);
                auto [bn, bd] = value(b);
                Big lhs = an * bd, rhs = bn * ad;
                if (lhs != rhs) return lhs < rhs;
                Wide as = a >= m ? 0 : (*edges)[a].dxdown;
                Wide bs = b >= m ? 0 : (*edges)[b].dxdown;
                Wide da = a >= m ? 1 : (*edges)[a].dy;
                Wide db = b >= m ? 1 : (*edges)[b].dy;
                Big slope_diff = Big(as) * db - Big(bs) * da;
                return slope_diff != 0 ? slope_diff < 0 : a < b;
            }
        };
        std::set<int, SweepCompare> active(SweepCompare{&sweep_edges, &x, &sweep_y, m});
        std::vector<typename std::set<int, SweepCompare>::iterator> active_it(m, active.end());
        std::vector<int> helper(m, -1), type(m);
        std::vector<std::pair<int, int>> diagonals;
        std::vector<int> order(m);
        std::iota(order.begin(), order.end(), 0);
        std::sort(order.begin(), order.end(), [&](int a, int b) { return y[a] > y[b]; });
        for (int i = 0; i < m; ++i) {
            int prev = (i + m - 1) % m, next = (i + 1) % m;
            bool prev_below = y[prev] < y[i], next_below = y[next] < y[i];
            bool convex = cross(points_[ring[prev]], points_[ring[i]], points_[ring[next]]) > 0;
            if (prev_below && next_below) type[i] = convex ? 0 : 1; // start / split
            else if (!prev_below && !next_below) type[i] = convex ? 2 : 3; // end / merge
            else type[i] = 4; // regular
        }
        auto add_diagonal = [&](int a, int b) {
            assert(a != b);
            diagonals.push_back(edge_key(a, b));
        };
        auto insert_edge = [&](int edge, int vertex) {
            assert(y[vertex] > y[(edge + 1) % m]);
            auto [it, inserted] = active.insert(edge);
            assert(inserted);
            active_it[edge] = it;
            helper[edge] = vertex;
        };
        auto erase_edge = [&](int edge) {
            assert(active_it[edge] != active.end());
            active.erase(active_it[edge]);
            active_it[edge] = active.end();
        };
        auto left_edge = [&](int vertex) {
            auto it = active.lower_bound(m + vertex);
            assert(it != active.begin());
            return *std::prev(it);
        };
        for (int v : order) {
            sweep_y = y[v];
            int before = (v + m - 1) % m;
            switch (type[v]) {
            case 0: // start
                insert_edge(v, v);
                break;
            case 1: { // split
                int e = left_edge(v);
                add_diagonal(v, helper[e]);
                helper[e] = v;
                insert_edge(v, v);
                break;
            }
            case 2: // end
                if (type[helper[before]] == 3)
                    add_diagonal(v, helper[before]);
                erase_edge(before);
                break;
            case 3: { // merge
                if (type[helper[before]] == 3)
                    add_diagonal(v, helper[before]);
                erase_edge(before);
                int e = left_edge(v);
                if (type[helper[e]] == 3)
                    add_diagonal(v, helper[e]);
                helper[e] = v;
                break;
            }
            case 4: // regular
                if (y[(v + 1) % m] < y[v]) {
                    if (type[helper[before]] == 3)
                        add_diagonal(v, helper[before]);
                    erase_edge(before);
                    insert_edge(v, v);
                } else {
                    int e = left_edge(v);
                    if (type[helper[e]] == 3)
                        add_diagonal(v, helper[e]);
                    helper[e] = v;
                }
                break;
            }
        }
        assert(active.empty());
        std::vector<std::pair<int, int>> graph_edges;
        graph_edges.reserve(m + diagonals.size());
        for (int i = 0; i < m; ++i)
            graph_edges.push_back(edge_key(i, (i + 1) % m));
        graph_edges.insert(graph_edges.end(), diagonals.begin(), diagonals.end());
        std::sort(graph_edges.begin(), graph_edges.end());
        graph_edges.erase(std::unique(graph_edges.begin(), graph_edges.end()), graph_edges.end());
        struct Halfedge { int from, to, next = -1; };
        std::vector<Halfedge> halfedges;
        halfedges.reserve(2 * graph_edges.size());
        std::vector<std::vector<int>> outgoing(m);
        for (auto [a, b] : graph_edges) {
            int id = int(halfedges.size());
            halfedges.push_back({a, b});
            halfedges.push_back({b, a});
            outgoing[a].push_back(id);
            outgoing[b].push_back(id + 1);
        }
        auto upper = [](Wide dx, Wide dy) { return dy > 0 || (dy == 0 && dx > 0); };
        std::vector<int> position(halfedges.size());
        for (int v = 0; v < m; ++v) {
            auto &out = outgoing[v];
            std::sort(out.begin(), out.end(), [&](int a, int b) {
                P origin = points_[ring[v]];
                P pa = points_[ring[halfedges[a].to]], pb = points_[ring[halfedges[b].to]];
                Wide ax = pa.x - origin.x, ay = pa.y - origin.y;
                Wide bx = pb.x - origin.x, by = pb.y - origin.y;
                bool ha = upper(ax, ay), hb = upper(bx, by);
                if (ha != hb) return ha > hb;
                Wide turn = ax * by - ay * bx;
                if (turn != 0) return turn > 0;
                return ax * ax + ay * ay < bx * bx + by * by;
            });
            for (int i = 0; i < int(out.size()); ++i)
                position[out[i]] = i;
        }
        for (int id = 0; id < int(halfedges.size()); ++id) {
            const auto &out = outgoing[halfedges[id].to];
            int at = position[id ^ 1];
            halfedges[id].next = out[(at + int(out.size()) - 1) % out.size()];
        }
        std::vector<char> used(halfedges.size());
        std::vector<std::array<int, 3>> result;
        result.reserve(m - 2);
        for (int start = 0; start < int(halfedges.size()); ++start) {
            if (used[start]) continue;
            std::vector<int> face;
            int at = start;
            do {
                assert(!used[at]);
                used[at] = true;
                face.push_back(halfedges[at].from);
                at = halfedges[at].next;
            } while (at != start);
            Wide area = 0;
            for (int i = 0; i < int(face.size()); ++i) {
                P a = points_[ring[face[i]]], b = points_[ring[face[(i + 1) % face.size()]]];
                area += a.x * b.y - a.y * b.x;
            }
            if (area <= 0) continue; // 外側の面
            triangulate_monotone_face(face, y, ring, result);
        }
        assert(result.size() == std::size_t(m - 2));
        return result;
    }
    void build_hierarchy(std::vector<int> current) {
        const int vertex_count = int(points_.size());
        std::vector<char> fixed(vertex_count);
        for (int v : frame_)
            fixed[v] = true;
        while (true) {
            std::vector<std::vector<int>> adjacent(vertex_count), incident(vertex_count);
            for (int id : current) {
                const auto &v = nodes_[id].v;
                for (int j = 0; j < 3; ++j) {
                    incident[v[j]].push_back(id);
                    adjacent[v[j]].push_back(v[(j + 1) % 3]);
                    adjacent[v[j]].push_back(v[(j + 2) % 3]);
                }
            }
            std::vector<int> candidates;
            for (int v = 0; v < vertex_count; ++v) {
                auto &adj = adjacent[v];
                std::sort(adj.begin(), adj.end());
                adj.erase(std::unique(adj.begin(), adj.end()), adj.end());
                if (!fixed[v] && !incident[v].empty() && adj.size() <= 12)
                    candidates.push_back(v);
            }
            std::sort(candidates.begin(), candidates.end(), [&](int a, int b) {
                return adjacent[a].size() < adjacent[b].size();
            });
            std::vector<char> blocked(vertex_count), selected(vertex_count);
            std::vector<int> removed;
            for (int v : candidates) {
                if (blocked[v])
                    continue;
                selected[v] = true;
                removed.push_back(v);
                blocked[v] = true;
                for (int u : adjacent[v])
                    blocked[u] = true;
            }
            if (removed.empty()) {
                top_ = std::move(current);
                break;
            }
            std::vector<int> next_level;
            next_level.reserve(current.size());
            for (int id : current) {
                const auto &v = nodes_[id].v;
                if (!selected[v[0]] && !selected[v[1]] && !selected[v[2]])
                    next_level.push_back(add_triangle(v, {id}));
            }
            for (int center : removed) {
                const auto &star = incident[center];
                const auto &neighbors = adjacent[center];
                std::vector<std::pair<int, int>> opposite;
                opposite.reserve(star.size());
                for (int id : star) {
                    const auto &v = nodes_[id].v;
                    int a = -1, b = -1;
                    for (int u : v)
                        if (u != center) {
                            if (a < 0) a = u;
                            else b = u;
                        }
                    opposite.push_back({a, b});
                }
                std::vector<int> ring;
                ring.reserve(neighbors.size());
                int first = neighbors[0], previous = -1, at = first;
                do {
                    ring.push_back(at);
                    int following = -1;
                    for (auto [a, b] : opposite) {
                        if (a == at && b != previous) following = b;
                        if (b == at && a != previous) following = a;
                    }
                    assert(following >= 0);
                    previous = at;
                    at = following;
                    assert(ring.size() <= neighbors.size());
                } while (at != first);
                assert(ring.size() == neighbors.size());
                for (auto tri : triangulate_quadratic(ring))
                    next_level.push_back(add_triangle(tri, star));
            }
            current = std::move(next_level);
        }
    }

public:
    explicit KirkpatrickPolygonPointLocator(const Polygon<T> &polygon) {
        assert(polygon.size() >= 3);
        for (auto p : polygon)
            points_.push_back({Wide(p.x), Wide(p.y)});
        // 180 度の頂点は境界を変えずに削除できる。
        int original_count = int(points_.size()), remaining = original_count;
        std::vector<int> previous(original_count), following(original_count);
        std::vector<char> alive(original_count, true);
        std::deque<int> candidates;
        for (int i = 0; i < original_count; ++i) {
            previous[i] = (i + original_count - 1) % original_count;
            following[i] = (i + 1) % original_count;
            candidates.push_back(i);
        }
        while (!candidates.empty() && remaining > 3) {
            int at = candidates.front();
            candidates.pop_front();
            if (!alive[at] || cross(points_[previous[at]], points_[at], points_[following[at]]) != 0)
                continue;
            int a = previous[at], b = following[at];
            alive[at] = false;
            following[a] = b;
            previous[b] = a;
            --remaining;
            candidates.push_back(a);
            candidates.push_back(b);
        }
        if (remaining != original_count) {
            std::vector<P> reduced;
            reduced.reserve(remaining);
            int first = 0;
            while (!alive[first]) ++first;
            int at = first;
            do {
                reduced.push_back(points_[at]);
                at = following[at];
            } while (at != first);
            points_ = std::move(reduced);
        }
        const int n = int(points_.size());
        std::vector<int> inner(n);
        for (int i = 0; i < n; ++i)
            inner[i] = i;
        if (signed_area(inner, points_) < 0)
            std::reverse(inner.begin(), inner.end());
        assert(signed_area(inner, points_) > 0);
        original_vertices_ = points_;
        std::sort(original_vertices_.begin(), original_vertices_.end());
        for (int i = 0; i < n; ++i)
            boundary_.push_back(edge_key(inner[i], inner[(i + 1) % n]));
        std::sort(boundary_.begin(), boundary_.end());
        Wide minx = points_[0].x, maxx = minx, miny = points_[0].y, maxy = miny;
        for (P p : points_) {
            minx = std::min(minx, p.x); maxx = std::max(maxx, p.x);
            miny = std::min(miny, p.y); maxy = std::max(maxy, p.y);
        }
        Wide margin = std::max(maxx - minx, maxy - miny) + 1;
        bottom_left_ = {minx - margin, miny - margin};
        top_right_ = {maxx + margin, maxy + margin};
        // 右下から反時計回り。橋の相手は元の多角形の最右頂点。
        frame_ = {n, n + 1, n + 2, n + 3};
        points_.push_back({top_right_.x, bottom_left_.y});
        points_.push_back(top_right_);
        points_.push_back({bottom_left_.x, top_right_.y});
        points_.push_back(bottom_left_);
        Wide shear = top_right_.x - bottom_left_.x + 1;
        std::vector<int> current;
        for (auto tri : triangulate_monotone(inner, shear))
            current.push_back(add_triangle(tri, {}, true));
        std::map<std::pair<Wide, Wide>, int> point_ids;
        for (int i = 0; i < n; ++i)
            point_ids[{points_[i].x, points_[i].y}] = i;
        std::vector<P> polygon_points(points_.begin(), points_.begin() + n);
        auto hull = convex_hull(polygon_points, false);
        std::vector<int> hull_ids;
        std::vector<char> seen(n);
        for (P p : hull) {
            int id = point_ids.at({p.x, p.y});
            if (!seen[id]) {
                seen[id] = true;
                hull_ids.push_back(id);
            }
        }
        std::vector<int> position(n);
        for (int i = 0; i < n; ++i)
            position[inner[i]] = i;
        std::sort(hull_ids.begin(), hull_ids.end(), [&](int a, int b) {
            return position[a] < position[b];
        });
        const int h = int(hull_ids.size());
        assert(h >= 3);
        auto append_exterior = [&](const std::vector<int> &domain) {
            for (auto tri : triangulate_monotone(domain, shear))
                current.push_back(add_triangle(tri));
        };
        // 凸包の各辺と元の多角形との間にあるポケットを分割する。
        for (int i = 0; i < h; ++i) {
            std::vector<int> chain;
            int at = position[hull_ids[i]], end = hull_ids[(i + 1) % h];
            while (true) {
                chain.push_back(inner[at]);
                if (inner[at] == end)
                    break;
                at = (at + 1) % n;
            }
            if (chain.size() < 3 || signed_area(chain, points_) == 0)
                continue;
            if (signed_area(chain, points_) < 0)
                std::reverse(chain.begin(), chain.end());
            append_exterior(chain);
        }
        // 凸包と外枠の間を、四隅から凸包の極値点への橋で4つに分ける。
        std::array<int, 4> extreme{0, 0, 0, 0};
        for (int i = 1; i < h; ++i) {
            P p = points_[hull_ids[i]];
            P r = points_[hull_ids[extreme[0]]];
            P t = points_[hull_ids[extreme[1]]];
            P l = points_[hull_ids[extreme[2]]];
            P b = points_[hull_ids[extreme[3]]];
            if (p.x > r.x || (p.x == r.x && p.y < r.y)) extreme[0] = i;
            if (p.y > t.y || (p.y == t.y && p.x > t.x)) extreme[1] = i;
            if (p.x < l.x || (p.x == l.x && p.y > l.y)) extreme[2] = i;
            if (p.y < b.y || (p.y == b.y && p.x < b.x)) extreme[3] = i;
        }
        int total_arcs = 0;
        for (int i = 0; i < 4; ++i)
            total_arcs += (extreme[(i + 1) % 4] - extreme[i] + h) % h;
        assert(total_arcs == h);
        for (int i = 0; i < 4; ++i) {
            std::vector<int> strip{frame_[i], frame_[(i + 1) % 4]};
            int at = extreme[(i + 1) % 4];
            while (true) {
                strip.push_back(hull_ids[at]);
                if (at == extreme[i])
                    break;
                at = (at + h - 1) % h;
            }
            assert(signed_area(strip, points_) > 0);
            append_exterior(strip);
        }
        assert(current.size() == std::size_t(2 * n + 2));
        build_hierarchy(std::move(current));
    }

    // 0: 外部、1: 境界、2: 内部。
    Int contains(Point<T> query) const {
        P p{Wide(query.x), Wide(query.y)};
        if (std::binary_search(original_vertices_.begin(), original_vertices_.end(), p))
            return 1;
        if (p.x < bottom_left_.x || p.x > top_right_.x ||
            p.y < bottom_left_.y || p.y > top_right_.y)
            return 0;
        int at = -1;
        for (int id : top_)
            if (in_triangle(nodes_[id].v, p)) {
                at = id;
                break;
            }
        if (at < 0)
            return 0;
        while (!nodes_[at].children.empty()) {
            int next = -1;
            for (int id : nodes_[at].children)
                if (in_triangle(nodes_[id].v, p)) {
                    next = id;
                    break;
                }
            assert(next >= 0);
            at = next;
        }
        const Triangle &leaf = nodes_[at];
        for (int i = 0; i < 3; ++i) {
            int a = leaf.v[i], b = leaf.v[(i + 1) % 3];
            if (std::binary_search(boundary_.begin(), boundary_.end(), edge_key(a, b)) &&
                on_segment(points_[a], points_[b], p))
                return 1;
        }
        return leaf.inside ? 2 : 0;
    }
};

template <class T> auto initKirkpatrickPolygonPointLocator(const Polygon<T> &polygon) {
    return KirkpatrickPolygonPointLocator<T>(polygon);
}

template <class T> Int contains(const KirkpatrickPolygonPointLocator<T> &locator, Point<T> p) {
    return locator.contains(p);
}

} // namespace cplib
