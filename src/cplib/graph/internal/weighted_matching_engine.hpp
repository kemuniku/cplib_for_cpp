#pragma once
#include <cplib/graph/graph.hpp>
#include <cplib/utils/itertools.hpp>
#include <queue>
#include <tuple>

namespace cplib {
struct WeightedMatchingResult {
    Int weight = 0;
    std::vector<std::pair<Int, Int>> matching;
};

// 最大重みマッチングの主双対探索。
// 花を固定順の子と巡回辺で表し、基点の変更は頂点列の並べ替えを伴わない。
// 疎グラフでは頂点列をAVL木で連結・分割し、双対値と最小余裕を遅延更新する。
// 密グラフでは所属配列と成分間の最小辺を使う。
namespace detail::weighted_matching {
inline constexpr Int infinity = std::numeric_limits<Int>::max();

struct Direction {
    Int src = 0, dst = 0;

    // 辺の向きを反転する。
    Direction reverse() const {
        return {dst, src};
    }
};

struct Edge {
    Int src, dst, weight;
};

enum Color { Hidden, Idle, Even, Odd };

// 外側・内側・未到達の頂点の双対値の変化率を返す。
inline Int slope(Color c) {
    return c == Even ? -1 : c == Odd ? 1 : 0;
}

struct Component {
    bool live = false;
    Int parent = 0, count = 0, first = 0, base = 0, pivot = 0, root = 0, version = 0;
    Color color = Hidden;
    Int dual = 0, changed = 0;
    Direction previous;
    std::vector<Int> children;
    std::vector<Direction> cycle;
    std::vector<Int> ends, members, direct;
    Int bestVertex = 0, crossEdge = -1, best = infinity, cross = infinity;
};

struct Leaf {
    Int left = 0, right = 0, parent = 0, height = 0, count = 0, owner = 0, potential = 0, lazy = 0,
        offer = 0, minimum = infinity, edge = -1, arg = 0;
};

struct Event {
    Int time = infinity, kind = -1, component = 0, version = 0, edge = -1;

    bool operator>(const Event &b) const {
        return std::tie(time, kind, component, version, edge) >
               std::tie(b.time, b.kind, b.component, b.version, b.edge);
    }
};

template <bool Fast> struct Machine {
    Int n, serial = 0, queueHead = 0, time = 0, currentVisit = 0;
    std::vector<Edge> arcs;
    std::vector<std::vector<Int>> adjacency;
    std::vector<Component> components;
    std::vector<Leaf> leaves;
    std::vector<Int> free, mate, queue, marked, owner;
    std::priority_queue<Event, std::vector<Event>, std::greater<Event>> heap;
    std::vector<Event> pending;
    bool heapReady = false;
    std::vector<std::vector<Int>> between;

    explicit Machine(Int count)
        : n(count), adjacency(n + 1), components(2 * n + 1), leaves(n + 1), mate(n + 1),
          marked(2 * n + 1) {
        if constexpr (!Fast) {
            owner.resize(n + 1);
            between.assign(2 * n + 1, std::vector<Int>(2 * n + 1, -1));
        }
    }

    // AVL部分木の双対値と最小余裕を一括更新する。O(1)。
    void apply(Int root, Int delta) {
        if (!root || !delta)
            return;
        leaves[root].potential += delta;
        leaves[root].lazy += delta;
        if (leaves[root].minimum != infinity)
            leaves[root].minimum += delta;
    }

    // 保留している双対値の更新を子へ伝える。O(1)。
    void push(Int root) {
        Int d = leaves[root].lazy;
        if (!d)
            return;
        apply(leaves[root].left, d);
        apply(leaves[root].right, d);
        leaves[root].lazy = 0;
    }

    // AVL節点の高さ・要素数・最小余裕を再計算する。O(1)。
    void pull(Int root) {
        Int l = leaves[root].left, r = leaves[root].right;
        auto &a = leaves[root];
        a.height = 1 + std::max(leaves[l].height, leaves[r].height);
        a.count = 1 + leaves[l].count + leaves[r].count;
        Int best = infinity, arg = 0;
        if (a.edge >= 0) {
            best = a.potential + a.offer;
            arg = root;
        }
        for (Int c : {l, r})
            if (c && leaves[c].minimum < best) {
                best = leaves[c].minimum;
                arg = leaves[c].arg;
            }
        a.minimum = best;
        a.arg = arg;
    }

    // 子へのリンクと逆向きの親リンクを同時に更新する。
    void attach(Int root, Int child, bool right) {
        if (right)
            leaves[root].right = child;
        else
            leaves[root].left = child;
        if (child)
            leaves[child].parent = root;
    }

    // 一回の回転で、頂点列の順序を保ちながらAVL木を組み替える。
    Int rotate(Int root, bool left) {
        push(root);
        Int out = left ? leaves[root].right : leaves[root].left;
        push(out);
        if (left) {
            attach(root, leaves[out].left, true);
            attach(out, root, false);
        } else {
            attach(root, leaves[out].right, false);
            attach(out, root, true);
        }
        leaves[out].parent = 0;
        pull(root);
        pull(out);
        return out;
    }

    // 左右の高さの差を解消する。O(1)。
    Int balance(Int root) {
        pull(root);
        leaves[root].parent = 0;
        Int l = leaves[root].left, r = leaves[root].right;
        if (leaves[l].height > leaves[r].height + 1) {
            if (leaves[leaves[l].right].height > leaves[leaves[l].left].height)
                attach(root, rotate(l, true), false);
            return rotate(root, false);
        }
        if (leaves[r].height > leaves[l].height + 1) {
            if (leaves[leaves[r].left].height > leaves[leaves[r].right].height)
                attach(root, rotate(r, false), true);
            return rotate(root, true);
        }
        return root;
    }

    // 一節点を間に挟んで二つのAVL木を連結する。O(高さの差+1)。
    Int join(Int l, Int middle, Int r) {
        if (leaves[l].height > leaves[r].height + 1) {
            push(l);
            attach(l, join(leaves[l].right, middle, r), true);
            return balance(l);
        }
        if (leaves[r].height > leaves[l].height + 1) {
            push(r);
            attach(r, join(l, middle, leaves[r].left), false);
            return balance(r);
        }
        attach(middle, l, false);
        attach(middle, r, true);
        leaves[middle].parent = 0;
        pull(middle);
        return middle;
    }

    // 頂点列を先頭count個と残りに分割する。O(log V)。
    std::pair<Int, Int> split(Int root, Int count) {
        if (!root)
            return {0, 0};
        push(root);
        Int l = leaves[root].left, r = leaves[root].right, size = leaves[l].count;
        leaves[root].left = leaves[root].right = leaves[root].parent = 0;
        std::pair<Int, Int> out;
        if (count <= size) {
            auto [a, b] = split(l, count);
            out = {a, join(b, root, r)};
        } else {
            auto [a, b] = split(r, count - size - 1);
            out = {join(l, root, a), b};
        }
        if (out.first)
            leaves[out.first].parent = 0;
        if (out.second)
            leaves[out.second].parent = 0;
        return out;
    }

    // 二つの頂点列を連結する。O(log V)。
    Int concatenate(Int l, Int r) {
        if (!l) {
            if (r)
                leaves[r].parent = 0;
            return r;
        }
        if (!r) {
            leaves[l].parent = 0;
            return l;
        }
        auto [a, b] = split(l, leaves[l].count - 1);
        return join(a, b, r);
    }

    // 頂点を含む最上位の花を返す。密版O(1)、疎版O(log V)。
    Int top(Int v) const {
        if constexpr (Fast) {
            while (leaves[v].parent)
                v = leaves[v].parent;
            return leaves[v].owner;
        } else
            return owner[v];
    }

    // 時刻の項を除いた頂点の双対値を返す。密版O(1)、疎版O(log V)。
    Int potential(Int v) const {
        Int out = leaves[v].potential;
        if constexpr (Fast) {
            v = leaves[v].parent;
            while (v) {
                out += leaves[v].lazy;
                v = leaves[v].parent;
            }
        }
        return out;
    }

    // 頂点の、所属するAVL木の中での位置を返す。O(log V)。
    Int rank(Int v) const {
        Int out = leaves[leaves[v].left].count;
        while (leaves[v].parent) {
            Int p = leaves[v].parent;
            if (leaves[p].right == v)
                out += leaves[leaves[p].left].count + 1;
            v = p;
        }
        return out;
    }

    // 指定頂点を含む直接の子の添字を返す。密版O(1)、疎版O(log V)。
    Int childIndex(Int b, Int v) const {
        if constexpr (Fast) {
            Int index = rank(v) - rank(components[b].first);
            return std::upper_bound(components[b].ends.begin(), components[b].ends.end(), index) -
                   components[b].ends.begin();
        } else
            return components[b].direct[v];
    }

    // 最上位の花の頂点を列挙する。O(頂点数)。
    Generator<Int> vertices(Int b) const {
        if constexpr (Fast) {
            std::vector<Int> stack;
            Int current = components[b].root;
            while (current || !stack.empty()) {
                if (current) {
                    stack.push_back(current);
                    current = leaves[current].left;
                } else {
                    current = stack.back();
                    stack.pop_back();
                    co_yield current;
                    current = leaves[current].right;
                }
            }
        } else
            for (Int v : components[b].members)
                co_yield v;
    }

    // 一つの花の頂点双対値を一括更新する。密版O(頂点数)、疎版O(1)。
    void shift(Int b, Int delta) {
        if constexpr (Fast)
            apply(components[b].root, delta);
        else {
            for (Int v : components[b].members)
                leaves[v].potential += delta;
            if (components[b].best != infinity)
                components[b].best += delta;
        }
    }

    // 花の双対値を現在の時刻にそろえる。
    void touch(Int b) {
        if (b > n)
            components[b].dual -= 2 * slope(components[b].color) * (time - components[b].changed);
        components[b].changed = time;
    }

    // 古いイベントを無効化するため、花の世代番号を更新する。
    void invalidate(Int b) {
        components[b].version = ++serial;
    }

    // 花へ入る辺のうち最小の余裕を与える候補を返す。O(1)。
    std::pair<Int, Int> bestOffer(Int b) const {
        if constexpr (Fast) {
            Int root = components[b].root;
            return {leaves[root].minimum, leaves[root].arg};
        } else
            return {components[b].best, components[b].bestVertex};
    }

    // 初回走査のイベントはまとめて蓄え、以後はヒープに追加する。
    void schedule(Event e) {
        if (heapReady)
            heap.push(e);
        else
            pending.push_back(e);
    }

    // 未到達の花に入る候補をイベントとして登録する。疎版O(log V)。
    void scheduleOffer(Int b) {
        if constexpr (Fast) {
            if (components[b].color != Idle)
                return;
            auto [value, v] = bestOffer(b);
            invalidate(b);
            if (v)
                schedule({value, 0, b, components[b].version, leaves[v].edge});
        }
    }

    // 外側頂点からの候補を更新する。密版O(1)、疎版O(log V)。
    void offer(Int v, Int edge, Int value) {
        if (leaves[v].edge >= 0 && value >= leaves[v].offer)
            return;
        Int b = top(v);
        if constexpr (Fast) {
            std::vector<Int> path{v};
            while (leaves[path.back()].parent)
                path.push_back(leaves[path.back()].parent);
            for (auto i = path.rbegin(); i != path.rend(); ++i)
                push(*i);
            leaves[v].offer = value;
            leaves[v].edge = edge;
            for (Int x : path)
                pull(x);
        } else {
            leaves[v].offer = value;
            leaves[v].edge = edge;
            Int key = leaves[v].potential + value;
            if (key < components[b].best) {
                components[b].best = key;
                components[b].bestVertex = v;
            }
        }
        scheduleOffer(b);
    }

    // 二つの外側頂点を結ぶ辺がタイトになる時刻を返す。
    Int edgeTime(Int id) const {
        auto e = arcs[id];
        return (potential(e.src) + potential(e.dst) - 2 * e.weight) / 2;
    }

    // 密版で、外側の花同士の最小辺を更新する。O(V)。
    void activateCross(Int b) {
        if constexpr (!Fast) {
            components[b].cross = infinity;
            components[b].crossEdge = -1;
            for (Int c = 1; c < Int(components.size()); ++c) {
                if (c == b || components[c].color != Even)
                    continue;
                Int edge = between[b][c];
                if (edge < 0)
                    continue;
                Int t = edgeTime(edge);
                if (t < components[b].cross) {
                    components[b].cross = t;
                    components[b].crossEdge = edge;
                }
                if (t < components[c].cross) {
                    components[c].cross = t;
                    components[c].crossEdge = edge;
                }
            }
        }
    }

    // 花のラベルと双対値の変化率を変更し、必要な探索を登録する。
    void label(Int b, Color color) {
        Color old = components[b].color;
        touch(b);
        shift(b, (slope(old) - slope(color)) * time);
        components[b].color = color;
        invalidate(b);
        if (color == Even && old != Even) {
            // 外側同士の候補は、登録した頂点の走査中にまとめて更新する。
            if constexpr (!Fast) {
                components[b].cross = infinity;
                components[b].crossEdge = -1;
            }
            for (Int v : vertices(b))
                queue.push_back(v);
        } else if (color == Idle)
            scheduleOffer(b);
        else if (color == Odd && b > n) {
            if constexpr (Fast)
                schedule({time + components[b].dual / 2, 2, b, components[b].version, -1});
        }
    }

    // 二本の交互路を交互にたどり、共通祖先を返す。異なる木なら0。
    Int ancestor(Int a, Int b) {
        ++currentVisit;
        while (a || b) {
            if (a) {
                if (marked[a] == currentVisit)
                    return a;
                marked[a] = currentVisit;
                auto e = components[a].previous;
                if (!e.src)
                    a = 0;
                else
                    a = top(components[top(e.src)].previous.src);
            }
            std::swap(a, b);
        }
        return 0;
    }

    struct Branch {
        std::vector<Int> nodes;
        std::vector<Direction> edges;
    };

    // 外側の花から祖先までの交互路を、上向きに列挙する。
    Branch branch(Int current, Int finish) const {
        Branch out;
        while (current != finish) {
            out.nodes.push_back(current);
            auto matched = components[current].previous;
            out.edges.push_back(matched.reverse());
            current = top(matched.src);
            out.nodes.push_back(current);
            auto unmatched = components[current].previous;
            out.edges.push_back(unmatched.reverse());
            current = top(unmatched.src);
        }
        return out;
    }

    // 交互路二本と一辺を奇閉路にし、花として縮約する。
    void contract(Int x, Int y, Int common) {
        auto a = branch(top(x), common), b = branch(top(y), common);
        Int id = free.back();
        free.pop_back();
        Component flower;
        flower.live = true;
        flower.base = components[common].base;
        flower.first = components[common].first;
        flower.previous = components[common].previous;
        flower.changed = time;
        flower.children.push_back(common);
        for (Int i = Int(a.nodes.size()) - 1; i >= 0; --i) {
            flower.cycle.push_back(a.edges[i].reverse());
            flower.children.push_back(a.nodes[i]);
        }
        flower.cycle.push_back({x, y});
        for (Int i = 0; i < Int(b.nodes.size()); ++i) {
            flower.children.push_back(b.nodes[i]);
            flower.cycle.push_back(b.edges[i]);
        }
        if constexpr (!Fast)
            flower.direct.resize(n + 1);
        for (Int i = 0; i < Int(flower.children.size()); ++i) {
            Int child = flower.children[i];
            label(child, Even);
            touch(child);
            components[child].color = Hidden;
            components[child].parent = id;
            invalidate(child);
            flower.count += components[child].count;
            flower.ends.push_back(flower.count);
            if constexpr (Fast)
                flower.root = concatenate(flower.root, components[child].root);
            else
                for (Int v : components[child].members) {
                    flower.members.push_back(v);
                    flower.direct[v] = i;
                    owner[v] = id;
                }
        }
        if constexpr (Fast)
            leaves[flower.root].owner = id;
        else
            for (Int other = 1; other < Int(components.size()); ++other) {
                if (!components[other].live)
                    continue;
                Int chosen = -1, value = infinity;
                for (Int child : flower.children) {
                    Int edge = between[child][other];
                    if (edge < 0)
                        continue;
                    auto e = arcs[edge];
                    Int key = potential(e.src) + potential(e.dst) - 2 * e.weight;
                    if (key < value) {
                        value = key;
                        chosen = edge;
                    }
                }
                between[id][other] = between[other][id] = chosen;
            }
        flower.color = Even;
        components[id] = std::move(flower);
        invalidate(id);
        activateCross(id);
    }

    // 巡回辺を指定された向きで返す。
    Direction cycleEdge(Int b, Int index, Int direction) const {
        if (direction == 1)
            return components[b].cycle[index];
        Int prev = (index + Int(components[b].children.size()) - 1) % components[b].children.size();
        return components[b].cycle[prev].reverse();
    }

    // 指定頂点を花の基点にし、内部の交互路を反転する。再帰は使わない。
    void expose(Int component, Int vertex) {
        struct Task {
            Int kind, a, b;
        };

        std::vector<Task> tasks{{0, component, vertex}};
        while (!tasks.empty()) {
            auto t = tasks.back();
            tasks.pop_back();
            if (t.kind == 1) {
                mate[t.a] = t.b;
                mate[t.b] = t.a;
                continue;
            }
            Int b = t.a, v = t.b;
            if (b <= n) {
                mate[v] = 0;
                continue;
            }
            Int target = childIndex(b, v), count = components[b].children.size(),
                at = components[b].pivot, distance = (target - at + count) % count,
                direction = distance % 2 == 0 ? 1 : -1;
            while (at != target) {
                Int next = (at + direction + count) % count;
                auto e = cycleEdge(b, at, direction);
                tasks.push_back({1, e.src, e.dst});
                tasks.push_back({0, components[b].children[next], e.dst});
                tasks.push_back({0, components[b].children[at], e.src});
                at = (next + direction + count) % count;
            }
            tasks.push_back({0, components[b].children[target], v});
            components[b].pivot = target;
            components[b].base = v;
        }
    }

    // 異なる根を結ぶ二本の交互路を記録し、そのマッチ辺を反転する。
    void augment(Int x, Int y) {
        std::vector<std::pair<Int, Int>> bases;
        std::vector<Direction> selected{{x, y}};
        for (Int endpoint : {x, y}) {
            Int v = endpoint;
            while (true) {
                Int b = top(v);
                bases.emplace_back(b, v);
                auto incoming = components[b].previous;
                if (!incoming.src)
                    break;
                Int inner = top(incoming.src);
                auto edge = components[inner].previous;
                bases.emplace_back(inner, edge.dst);
                selected.push_back(edge);
                v = edge.src;
            }
        }
        for (auto [b, v] : bases)
            expose(b, v);
        for (auto e : selected) {
            mate[e.src] = e.dst;
            mate[e.dst] = e.src;
        }
    }

    // 未到達の花と、そのマッチ先を交互木に追加する。
    void grow(Int x, Int y) {
        Int b = top(y), z = mate[components[b].base];
        assert(z);
        components[b].previous = {x, y};
        label(b, Odd);
        Int next = top(z);
        components[next].previous = {components[b].base, z};
        label(next, Even);
    }

    // 内側の花を分割し、入口から基点までの偶数長の交互路を復元する。
    void expand(Int b) {
        auto entry = components[b].previous;
        Int first = childIndex(b, entry.dst), finish = components[b].pivot,
            count = components[b].children.size(), distance = (finish - first + count) % count,
            direction = distance % 2 == 0 ? 1 : -1;
        std::vector<Color> colors(count, Idle);
        std::vector<Direction> incoming(count);
        Int at = first;
        colors[at] = Odd;
        incoming[at] = entry;
        while (at != finish) {
            Int next = (at + direction + count) % count;
            colors[next] = colors[at] == Odd ? Even : Odd;
            incoming[next] = cycleEdge(b, at, direction);
            at = next;
        }
        if constexpr (Fast) {
            Int remaining = components[b].root;
            for (Int child : components[b].children) {
                auto [l, r] = split(remaining, components[child].count);
                components[child].root = l;
                leaves[l].owner = child;
                remaining = r;
            }
        } else
            for (Int child : components[b].children) {
                components[child].best = infinity;
                components[child].bestVertex = 0;
                for (Int v : components[child].members) {
                    owner[v] = child;
                    if (leaves[v].edge >= 0) {
                        Int key = leaves[v].potential + leaves[v].offer;
                        if (key < components[child].best) {
                            components[child].best = key;
                            components[child].bestVertex = v;
                        }
                    }
                }
            }
        components[b].color = Hidden;
        invalidate(b);
        for (Int i = 0; i < count; ++i) {
            Int child = components[b].children[i];
            components[child].parent = 0;
            components[child].color = Odd;
            components[child].changed = time;
            components[child].previous = incoming[i];
            label(child, colors[i]);
        }
        components[b] = Component{};
        free.push_back(b);
    }

    // 次の有効なイベントを返す。密版O(V)、疎版は一取り出しO(log V)。
    Event nextEvent() {
        Event out;
        if constexpr (Fast) {
            if (!heapReady) {
                heap = decltype(heap)(std::greater<Event>{}, std::move(pending));
                pending.clear();
                heapReady = true;
            }
            while (!heap.empty()) {
                auto e = heap.top();
                heap.pop();
                if (e.kind == 1) {
                    auto edge = arcs[e.edge];
                    if (top(edge.src) != top(edge.dst))
                        return e;
                } else {
                    Int b = e.component;
                    if (!components[b].live || components[b].parent)
                        continue;
                    if (components[b].version != e.version)
                        continue;
                    if (components[b].color == (e.kind == 0 ? Idle : Odd))
                        return e;
                }
            }
        } else
            for (Int b = 1; b < Int(components.size()); ++b) {
                auto &c = components[b];
                if (c.color == Idle) {
                    auto [value, v] = bestOffer(b);
                    if (value < out.time)
                        out = {value, 0, b, 0, leaves[v].edge};
                } else if (c.color == Even) {
                    if (c.cross < out.time)
                        out = {c.cross, 1, b, 0, c.crossEdge};
                } else if (c.color == Odd && b > n) {
                    Int t = c.changed + c.dual / 2;
                    if (t < out.time)
                        out = {t, 2, b, 0, -1};
                }
            }
        return out;
    }

    // 新しく外側になった頂点を走査し、時刻を進めずに使える辺を直ちに処理する。
    bool scan() {
        while (queueHead < Int(queue.size())) {
            Int u = queue[queueHead++], source = top(u), p = potential(u);
            // 頂点から出る辺番号を列挙する。密版では連続配置を利用する。
            for (Int id : adjacency[u]) {
                auto edge = arcs[id];
                Int dest = top(edge.dst);
                if (source == dest)
                    continue;
                if (components[dest].color == Even) {
                    Int t = (p + potential(edge.dst) - 2 * edge.weight) / 2;
                    if (t == time) {
                        Int common = ancestor(source, dest);
                        if (!common) {
                            augment(u, edge.dst);
                            return true;
                        }
                        contract(u, edge.dst, common);
                        source = top(u);
                    } else {
                        if constexpr (Fast)
                            schedule({t, 1, 0, 0, id});
                        else {
                            if (t < components[source].cross) {
                                components[source].cross = t;
                                components[source].crossEdge = id;
                            }
                            if (t < components[dest].cross) {
                                components[dest].cross = t;
                                components[dest].crossEdge = id;
                            }
                        }
                    }
                } else {
                    offer(edge.dst, id, p - 2 * edge.weight);
                    if (components[dest].color == Idle &&
                        p + potential(edge.dst) - 2 * edge.weight == time)
                        grow(u, edge.dst);
                }
            }
        }
        return false;
    }

    // 双対値を現在時刻で確定し、次の探索の時刻原点へ移す。
    void finishStage() {
        for (Int b = 1; b < Int(components.size()); ++b) {
            if (!components[b].live || components[b].parent)
                continue;
            touch(b);
            shift(b, slope(components[b].color) * time);
            components[b].color = Idle;
            components[b].changed = 0;
        }
        time = 0;
    }

    // 増加路を一つ求める。密版O(V^2)、疎版O(E log V)。
    bool stage() {
        queue.clear();
        queueHead = 0;
        heap = {};
        pending.clear();
        heapReady = false;
        for (Int v = 1; v <= n; ++v) {
            leaves[v].edge = -1;
            leaves[v].minimum = infinity;
            leaves[v].arg = 0;
        }
        for (Int b = 1; b < Int(components.size()); ++b) {
            components[b].best = infinity;
            components[b].bestVertex = 0;
            components[b].cross = infinity;
            components[b].crossEdge = -1;
            components[b].previous = {};
        }
        Int deadline = infinity;
        for (Int b = 1; b < Int(components.size()); ++b)
            if (components[b].color == Idle) {
                Int v = components[b].base;
                if (!mate[v]) {
                    deadline = std::min(deadline, potential(v));
                    label(b, Even);
                }
            }
        if (deadline == infinity)
            return false;
        while (true) {
            if (scan()) {
                finishStage();
                return true;
            }
            auto event = nextEvent();
            if (event.time >= deadline) {
                time = deadline;
                finishStage();
                return false;
            }
            assert(event.time >= time);
            time = event.time;
            if (event.kind == 0) {
                auto e = arcs[event.edge];
                grow(e.src, e.dst);
            } else if (event.kind == 1) {
                auto e = arcs[event.edge];
                Int common = ancestor(top(e.src), top(e.dst));
                if (!common) {
                    augment(e.src, e.dst);
                    finishStage();
                    return true;
                }
                contract(e.src, e.dst, common);
            } else if (event.kind == 2)
                expand(event.component);
            else
                assert(false);
        }
    }
};

// 整数重みの最大重みマッチングを求める。公開モジュール二種類から利用する共通実装。
template <bool Fast, UnDirectedGraph G> WeightedMatchingResult solve(const G &g) {
    static_assert(std::is_integral_v<typename G::cost_type> &&
                  std::is_signed_v<typename G::cost_type>);
    if constexpr (G::is_static)
        g.static_graph_initialized_check();
    std::vector<Int> indices(g.len), vertices{-1};
    Int largest = 0;
    for (auto e : g.edge_info) {
        if (e.src == e.dst || e.cost <= 0)
            continue;
        largest = std::max(largest, Int(e.cost));
        for (Int v : {e.src, e.dst})
            if (!indices[v]) {
                indices[v] = vertices.size();
                vertices.push_back(v);
            }
    }
    Int n = vertices.size() - 1;
    WeightedMatchingResult out;
    if (!n)
        return out;
    assert(largest <= infinity / 4 / n);
    Machine<Fast> s(n);
    for (Int b = 2 * n; b > n; --b)
        s.free.push_back(b);
    for (Int v = 1; v <= n; ++v) {
        auto &c = s.components[v];
        c.live = true;
        c.count = 1;
        c.first = c.base = c.root = v;
        c.color = Idle;
        auto &leaf = s.leaves[v];
        leaf.height = leaf.count = 1;
        leaf.owner = v;
        leaf.potential = largest;
        if constexpr (!Fast) {
            c.members = {v};
            s.owner[v] = v;
        }
    }
    if constexpr (Fast) {
        std::vector<std::vector<std::pair<Int, Int>>> raw(n + 1);
        for (auto e : g.edge_info) {
            if (e.src == e.dst || e.cost <= 0)
                continue;
            Int u = indices[e.src], v = indices[e.dst];
            raw[u].emplace_back(v, e.cost);
            raw[v].emplace_back(u, e.cost);
        }
        std::vector<Int> seen(n + 1), position(n + 1);
        for (Int u = 1; u <= n; ++u)
            for (auto [v, w] : raw[u]) {
                if (seen[v] != u) {
                    seen[v] = u;
                    position[v] = s.arcs.size();
                    s.adjacency[u].push_back(s.arcs.size());
                    s.arcs.push_back({u, v, w});
                } else
                    s.arcs[position[v]].weight = std::max(s.arcs[position[v]].weight, w);
            }
    } else {
        for (auto e : g.edge_info) {
            if (e.src == e.dst || e.cost <= 0)
                continue;
            Int a = indices[e.src], b = indices[e.dst];
            for (auto [u, v] : {std::pair{a, b}, std::pair{b, a}}) {
                Int id = s.between[u][v];
                if (id < 0) {
                    s.between[u][v] = s.arcs.size();
                    s.arcs.push_back({u, v, Int(e.cost)});
                } else
                    s.arcs[id].weight = std::max(s.arcs[id].weight, Int(e.cost));
            }
        }
        // 同じ頂点から出る辺を行先順に連続配置し、走査時の局所性を高める。
        std::vector<Edge> ordered;
        ordered.reserve(s.arcs.size());
        for (Int u = 1; u <= n; ++u)
            for (Int v = 1; v <= n; ++v) {
                Int id = s.between[u][v];
                if (id < 0)
                    continue;
                s.between[u][v] = ordered.size();
                s.adjacency[u].push_back(ordered.size());
                ordered.push_back(s.arcs[id]);
            }
        s.arcs = std::move(ordered);
    }
    // 最大重みの辺だけで作る初期マッチングは、その辺数に対して既に最適。
    for (auto e : s.arcs) {
        if (e.weight == largest && !s.mate[e.src] && !s.mate[e.dst]) {
            s.mate[e.src] = e.dst;
            s.mate[e.dst] = e.src;
        }
    }
    while (s.stage()) {
    }
    for (Int u = 1; u <= n; ++u) {
        Int v = s.mate[u];
        if (u >= v)
            continue;
        for (Int id : s.adjacency[u])
            if (s.arcs[id].dst == v) {
                out.weight += s.arcs[id].weight;
                break;
            }
        out.matching.emplace_back(std::min(vertices[u], vertices[v]),
                                  std::max(vertices[u], vertices[v]));
    }
    return out;
}
}
}
