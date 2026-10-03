#pragma once
#include <cplib/collections/private/retroactive_common.hpp>
#include <functional>
#include <variant>

namespace cplib::detail {
// Shared AVL engine: either ordered timestamps or implicit operation positions.
template <class K, class T, class S, bool Indexed> class RetroactiveAVL {
    using Aggregate = std::conditional_t<std::is_void_v<S>, std::monostate, S>;

    struct Node {
        Int left = 0, right = 0, parent = 0, height = 0, size = 0, pushes = 0, pops = 0,
            operationId = 0;
        K time{};
        T value{};
        QueueDebugKind operation = qdkNone;
        bool remains = false;
        Int balance = 0, minPrefix = 0, remaining = 0, removed = 0;
        Aggregate mapped{}, aggregate{};
    };

    std::vector<Node> nodes = std::vector<Node>(2);
    std::vector<Int> free, idNodes;
    Int root = 1, count_ = 0, dummyRemoved = 0;
    SortOrder order;
    T total{}, pushTotal{};
    std::function<Aggregate(Aggregate, Aggregate)> op;
    std::function<Aggregate(T)> lift;
    Aggregate identity{};

    Int nodeRank(Int i) const {
        Int out = nodes[nodes[i].left].size;
        while (nodes[i].parent) {
            Int p = nodes[i].parent;
            if (nodes[p].right == i)
                out += nodes[nodes[p].left].size + 1;
            i = p;
        }
        return out;
    }

    bool before(Int a, Int b) const {
        if (a == 1)
            return false;
        if (b == 1)
            return true;
        if (nodes[a].value < nodes[b].value)
            return order == Ascending;
        if (nodes[b].value < nodes[a].value)
            return order == Descending;
        if constexpr (Indexed)
            return nodeRank(a) < nodeRank(b);
        else
            return nodes[a].time < nodes[b].time;
    }

    Int choose(Int a, Int b, bool remains) const {
        if (!a)
            return b;
        if (!b)
            return a;
        bool first;
        if (a == 1)
            first = false;
        else if (b == 1)
            first = true;
        else if (nodes[a].value < nodes[b].value)
            first = order == Ascending;
        else if (nodes[b].value < nodes[a].value)
            first = order == Descending;
        else
            first = true;
        return first == remains ? a : b;
    }

    Int ownBalance(Int i) const {
        return i == 1                                               ? dummyRemoved
               : nodes[i].operation == qdkPop                       ? -1
               : nodes[i].operation == qdkPush && !nodes[i].remains ? 1
                                                                    : 0;
    }

    Int ownCandidate(Int i, bool remains) const {
        if (i == 1) {
            if (remains || dummyRemoved > 0)
                return 1;
        } else if (nodes[i].operation == qdkPush && nodes[i].remains == remains)
            return i;
        return 0;
    }

    void pull(Int i) {
        Int l = nodes[i].left, r = nodes[i].right;
        if (l)
            nodes[l].parent = i;
        if (r)
            nodes[r].parent = i;
        auto &a = nodes[i];
        a.pushes = nodes[l].pushes + nodes[r].pushes + (a.operation == qdkPush);
        a.pops = nodes[l].pops + nodes[r].pops + (a.operation == qdkPop);
        Int mid = nodes[l].balance + ownBalance(i);
        a.height = std::max(nodes[l].height, nodes[r].height) + 1;
        a.size = nodes[l].size + 1 + nodes[r].size;
        a.balance = mid + nodes[r].balance;
        a.minPrefix = mid;
        if (l)
            a.minPrefix = std::min(a.minPrefix, nodes[l].minPrefix);
        if (r)
            a.minPrefix = std::min(a.minPrefix, mid + nodes[r].minPrefix);
        a.remaining = choose(choose(nodes[l].remaining, ownCandidate(i, true), true),
                             nodes[r].remaining, true);
        a.removed = choose(choose(nodes[l].removed, ownCandidate(i, false), false),
                           nodes[r].removed, false);
        if constexpr (!std::is_void_v<S>)
            a.aggregate =
                op(op(nodes[l].aggregate,
                      i != 1 && a.operation == qdkPush && a.remains ? a.mapped : identity),
                   nodes[r].aggregate);
    }

    Int rotateLeft(Int i) {
        Int out = nodes[i].right;
        nodes[i].right = nodes[out].left;
        nodes[out].left = i;
        pull(i);
        pull(out);
        return out;
    }

    Int rotateRight(Int i) {
        Int out = nodes[i].left;
        nodes[i].left = nodes[out].right;
        nodes[out].right = i;
        pull(i);
        pull(out);
        return out;
    }

    Int rebalance(Int i) {
        pull(i);
        Int l = nodes[i].left, r = nodes[i].right;
        if (nodes[l].height > nodes[r].height + 1) {
            if (nodes[nodes[l].left].height < nodes[nodes[l].right].height)
                nodes[i].left = rotateLeft(l);
            return rotateRight(i);
        }
        if (nodes[r].height > nodes[l].height + 1) {
            if (nodes[nodes[r].right].height < nodes[nodes[r].left].height)
                nodes[i].right = rotateRight(r);
            return rotateLeft(i);
        }
        return i;
    }

    Int insertNode(Int r, Int i, Int position) {
        if (!r)
            return i;
        Int mid = nodes[nodes[r].left].size;
        if (position <= mid)
            nodes[r].left = insertNode(nodes[r].left, i, position);
        else
            nodes[r].right = insertNode(nodes[r].right, i, position - mid - 1);
        return rebalance(r);
    }

    Int detachMin(Int r, Int &minimum) {
        if (!nodes[r].left) {
            minimum = r;
            return nodes[r].right;
        }
        nodes[r].left = detachMin(nodes[r].left, minimum);
        return rebalance(r);
    }

    Int deleteNode(Int r, Int position) {
        Int mid = nodes[nodes[r].left].size;
        if (position == mid) {
            Int l = nodes[r].left, rr = nodes[r].right;
            if (!l)
                return rr;
            if (!rr)
                return l;
            Int successor = 0, rest = detachMin(rr, successor);
            nodes[successor].left = l;
            nodes[successor].right = rest;
            return rebalance(successor);
        }
        if (position < mid)
            nodes[r].left = deleteNode(nodes[r].left, position);
        else
            nodes[r].right = deleteNode(nodes[r].right, position - mid - 1);
        return rebalance(r);
    }

    void refresh(Int i) {
        while (i) {
            pull(i);
            i = nodes[i].parent;
        }
    }

    Int nodeAt(Int position) const {
        Int i = root;
        while (i) {
            Int mid = nodes[nodes[i].left].size;
            if (position < mid)
                i = nodes[i].left;
            else if (position == mid)
                return i;
            else {
                position -= mid + 1;
                i = nodes[i].right;
            }
        }
        return 0;
    }

    std::pair<Int, Int> locate(const K &t) const {
        Int i = root, rank = 0;
        while (i) {
            if (i == 1 || nodes[i].time < t) {
                rank += nodes[nodes[i].left].size + 1;
                i = nodes[i].right;
            } else if (t < nodes[i].time)
                i = nodes[i].left;
            else
                return {i, rank + nodes[nodes[i].left].size};
        }
        return {0, rank};
    }

    Int addEmpty(Int position, K t = {}) {
        Int i;
        if (free.empty()) {
            i = nodes.size();
            nodes.emplace_back();
        } else {
            i = free.back();
            free.pop_back();
        }
        nodes[i].time = t;
        if constexpr (Indexed) {
            nodes[i].operationId = idNodes.size();
            idNodes.push_back(i);
        }
        pull(i);
        root = insertNode(root, i, position);
        nodes[root].parent = 0;
        return i;
    }

    std::pair<Int, Int> ensure(const K &t) {
        auto [i, rank] = locate(t);
        if (!i)
            i = addEmpty(rank, t);
        return {i, rank};
    }

    Int findBridge(Int i, Int start, Int t, Int prefix, bool first) const {
        if (!i)
            return -1;
        Int finish = start + nodes[i].size;
        if (first ? finish < t : start >= t)
            return -1;
        if (prefix + nodes[i].minPrefix > 0)
            return -1;
        Int l = nodes[i].left, r = nodes[i].right, mid = start + nodes[l].size,
            after = prefix + nodes[l].balance + ownBalance(i), out;
        if (first) {
            out = findBridge(l, start, t, prefix, true);
            if (out >= 0)
                return out;
            if (mid + 1 >= t && after == 0)
                return mid + 1;
            return findBridge(r, mid + 1, t, after, true);
        }
        out = findBridge(r, mid + 1, t, after, false);
        if (out >= 0)
            return out;
        if (mid + 1 <= t && after == 0)
            return mid + 1;
        return findBridge(l, start, t, prefix, false);
    }

    Int candidate(Int i, Int start, Int l, Int r, bool remains) const {
        if (!i || r <= start || start + nodes[i].size <= l)
            return 0;
        if (l <= start && start + nodes[i].size <= r)
            return remains ? nodes[i].remaining : nodes[i].removed;
        Int mid = start + nodes[nodes[i].left].size;
        Int out = candidate(nodes[i].left, start, l, r, remains);
        if (l <= mid && mid < r)
            out = choose(out, ownCandidate(i, remains), remains);
        return choose(out, candidate(nodes[i].right, mid + 1, l, r, remains), remains);
    }

    void changeRemaining(Int i, bool remains, QueueDelta<K, T> &d) {
        assert(i);
        if (i == 1)
            dummyRemoved += remains ? -1 : 1;
        else {
            nodes[i].remains = remains;
            if (remains) {
                ++count_;
                if constexpr (std::is_void_v<S> && std::is_arithmetic_v<T>)
                    total += nodes[i].value;
            } else {
                --count_;
                if constexpr (std::is_void_v<S> && std::is_arithmetic_v<T>)
                    total -= nodes[i].value;
            }
            if constexpr (!Indexed) {
                QueueValue<K, T> entry{nodes[i].time, nodes[i].value};
                if (remains)
                    d.added.push_back(entry);
                else {
                    auto it = std::find_if(d.added.begin(), d.added.end(), [&](auto &a) {
                        return !(a.time < entry.time) && !(entry.time < a.time);
                    });
                    if (it != d.added.end())
                        d.added.erase(it);
                    else
                        d.removed.push_back(entry);
                }
            }
        }
        refresh(i);
    }

    void eraseOperation(Int i, Int rank, QueueDelta<K, T> &d) {
        auto kind = nodes[i].operation;
        if (kind == qdkNone)
            return;
        if (kind == qdkPush) {
            if constexpr (std::is_void_v<S> && std::is_arithmetic_v<T>)
                pushTotal = queue_sub(pushTotal, nodes[i].value);
            if (nodes[i].remains)
                changeRemaining(i, false, d);
            else
                changeRemaining(candidate(root, 0, 0, findBridge(root, 0, rank + 1, 0, true), true),
                                false, d);
        } else {
            Int bridge = std::max<Int>(0, findBridge(root, 0, rank, 0, false));
            changeRemaining(candidate(root, 0, bridge, nodes[root].size, false), true, d);
        }
        nodes[i].operation = qdkNone;
        nodes[i].value = T{};
        nodes[i].mapped = Aggregate{};
        nodes[i].remains = false;
        refresh(i);
    }

    QueueDelta<K, T> pushAt(Int i, Int rank, T value) {
        QueueDelta<K, T> d;
        eraseOperation(i, rank, d);
        Int bridge = std::max<Int>(0, findBridge(root, 0, rank, 0, false)),
            x = candidate(root, 0, bridge, nodes[root].size, false);
        nodes[i].operation = qdkPush;
        nodes[i].value = value;
        if constexpr (std::is_void_v<S> && std::is_arithmetic_v<T>)
            pushTotal = queue_add(pushTotal, value);
        if constexpr (!std::is_void_v<S>)
            nodes[i].mapped = lift(value);
        if (!x || before(x, i))
            changeRemaining(i, true, d);
        else {
            nodes[i].remains = false;
            refresh(i);
            changeRemaining(x, true, d);
        }
        return d;
    }

    QueueDelta<K, T> popAt(Int i, Int rank) {
        QueueDelta<K, T> d;
        if (nodes[i].operation == qdkPop)
            return d;
        eraseOperation(i, rank, d);
        Int bridge = findBridge(root, 0, rank, 0, true);
        changeRemaining(candidate(root, 0, 0, bridge, true), false, d);
        nodes[i].operation = qdkPop;
        refresh(i);
        return d;
    }

    QueueDelta<K, T> eraseAt(Int i, Int rank) {
        QueueDelta<K, T> d;
        eraseOperation(i, rank, d);
        root = deleteNode(root, rank);
        nodes[root].parent = 0;
        if constexpr (Indexed)
            idNodes[nodes[i].operationId] = 0;
        nodes[i] = {};
        free.push_back(i);
        return d;
    }

    Int operationIndex(Int k, bool push) const {
        assert(0 <= k && k < (push ? pushCount() : popCount()));
        Int i = root, start = 0;
        while (i) {
            Int l = nodes[i].left, c = push ? nodes[l].pushes : nodes[l].pops;
            if (k < c) {
                i = l;
                continue;
            }
            k -= c;
            Int mid = start + nodes[l].size;
            if (nodes[i].operation == (push ? qdkPush : qdkPop)) {
                if (!k)
                    return mid - 1;
                --k;
            }
            start = mid + 1;
            i = nodes[i].right;
        }
        return -1;
    }

public:
    // 空の操作列を生成します。Ascendingはpop min、Descendingはpop maxです。
    // KとTには一貫した<が必要です。
    explicit RetroactiveAVL(SortOrder ord = Ascending)
        requires std::is_void_v<S>
        : order(ord) {
        pull(1);
    }

    // モノイド集約付きの空の操作列を生成します。
    // operationは結合的、eは単位元、operation・mappingは副作用なしとしてください。
    // 残存pushを時刻順（添字版は現在の操作順）で集約します。
    template <class Op, class Lift>
    RetroactiveAVL(Op operation, Aggregate e, Lift mapping, SortOrder ord = Ascending)
        requires(!std::is_void_v<S>)
        : order(ord), op(operation), lift(mapping), identity(e) {
        nodes[0].aggregate = e;
        pull(1);
    }

    // 時刻版では時刻tにpushを挿入・上書きし、最終状態の差分を返します。
    // 添字版では既存のt番目の操作をpushで上書きします。
    QueueDelta<K, T> setPush(K t, T v) {
        if constexpr (Indexed) {
            assert(0 <= t && t < operationCount());
            return pushAt(nodeAt(t + 1), t + 1, v);
        } else {
            auto [i, rank] = ensure(t);
            return pushAt(i, rank, v);
        }
    }

    // 時刻版では時刻tにpopを挿入・上書きし、最終状態の差分を返します。
    // 添字版では既存のt番目の操作をpopで上書きします。
    QueueDelta<K, T> setPop(K t) {
        if constexpr (Indexed) {
            assert(0 <= t && t < operationCount());
            return popAt(nodeAt(t + 1), t + 1);
        } else {
            auto [i, rank] = ensure(t);
            return popAt(i, rank);
        }
    }

    // 時刻版では時刻tの操作を削除します。未登録なら何もしません。
    // 添字版では既存のt番目の操作を削除し、後続の位置は1つ前にずれます。
    QueueDelta<K, T> erase(K t) {
        if constexpr (Indexed) {
            assert(0 <= t && t < operationCount());
            return eraseAt(nodeAt(t + 1), t + 1);
        } else {
            auto [i, rank] = locate(t);
            return i ? eraseAt(i, rank) : QueueDelta<K, T>{};
        }
    }

    // 操作列の長さを返します。キュー内の要素数lenとは異なります。O(1)。
    Int operationCount() const {
        return nodes[root].size - 1;
    }

    // 操作列内のpushの個数を返します。取り出し済みのpushも含みます。O(1)。
    Int pushCount() const {
        return nodes[root].pushes;
    }

    // 操作列内のpopの個数を返します。空へのpopも含みます。O(1)。
    Int popCount() const {
        return nodes[root].pops;
    }

    // indexの直前にpushを挿入し、操作IDを返します。index=operationCountなら末尾です。償却O(log(M+2))。
    Int insertPush(Int index, T value)
        requires Indexed
    {
        assert(0 <= index && index <= operationCount());
        Int i = addEmpty(index + 1);
        pushAt(i, index + 1, value);
        return nodes[i].operationId;
    }

    // indexの直前にpopを挿入し、操作IDを返します。index=operationCountなら末尾です。償却O(log(M+2))。
    Int insertPop(Int index)
        requires Indexed
    {
        assert(0 <= index && index <= operationCount());
        Int i = addEmpty(index + 1);
        popAt(i, index + 1);
        return nodes[i].operationId;
    }

    // 現在の操作列でk番目のpushの位置を返します。kは0始まりです。O(log(M+2))。
    Int pushIndex(Int k) const
        requires Indexed
    {
        return operationIndex(k, true);
    }

    // 現在の操作列でk番目のpopの位置を返します。空へのpopも数えます。O(log(M+2))。
    Int popIndex(Int k) const
        requires Indexed
    {
        return operationIndex(k, false);
    }

    Int insertPushBeforePop(Int k, T v)
        requires Indexed
    {
        return insertPush(popIndex(k), v);
    }

    Int insertPopBeforePush(Int k)
        requires Indexed
    {
        return insertPop(pushIndex(k));
    }

    // 操作IDの現在位置を返します。削除済み・未発行IDなら-1です。O(log(M+2))。
    Int indexOf(Int id) const
        requires Indexed
    {
        return id < 0 || id >= Int(idNodes.size()) || !idNodes[id] ? -1 : nodeRank(idNodes[id]) - 1;
    }

    // 操作IDで削除します。削除済み・未発行IDなら何もしません。償却O(log(M+2))。
    void eraseById(Int id)
        requires Indexed
    {
        Int i = indexOf(id);
        if (i >= 0)
            erase(i);
    }

    // 全操作実行後に残る実要素数を返します。O(1)。
    Int len() const {
        return count_;
    }

    T sum() const
        requires(std::is_void_v<S> && std::is_arithmetic_v<T>)
    {
        return total;
    }

    // 現在の操作列でpopされる値の総和を返します。空へのpopは0として扱います。O(1)。
    // 整数は結果が型に収まる必要があります。内部の全push総和は桁あふれを許容します。
    // 浮動小数点は全push総和から残存総和を引くため、桁落ちが生じる場合があります。
    T poppedSum() const
        requires(std::is_void_v<S> && std::is_arithmetic_v<T>)
    {
        return queue_sub(pushTotal, total);
    }

    // 残存pushを時刻順（添字版は現在の操作順）で集約したモノイド積を返します。空なら単位元です。
    Aggregate fold() const
        requires(!std::is_void_v<S>)
    {
        return nodes[root].aggregate;
    }

    Aggregate get_all() const
        requires(!std::is_void_v<S>)
    {
        return fold();
    }

    // 最終状態の最優先要素を返します。空ならnoneです。O(1)。
    std::optional<T> peek() const {
        Int i = nodes[root].remaining;
        return i <= 1 ? std::optional<T>{} : nodes[i].value;
    }

    // 時刻tのpushが最後に残るかを返します。時刻版では未登録・popならfalseです。
    // 添字版ではtに既存の操作位置を指定してください。
    bool isRemaining(K t) const {
        Int i;
        if constexpr (Indexed) {
            assert(0 <= t && t < operationCount());
            i = nodeAt(t + 1);
        } else
            i = locate(t).first;
        return i && nodes[i].operation == qdkPush && nodes[i].remains;
    }

    // 登録中の操作を時刻順で返します。番兵・削除済み時刻を除き、pop結果は未計算です。O(M)。
    // 添字版では現在の操作位置をtimeとして返します。
    std::vector<QueueDebugEntry<K, T>> debugOperations() const {
        std::vector<QueueDebugEntry<K, T>> out;
        std::vector<Int> stack;
        Int i = root;
        while (i || !stack.empty()) {
            while (i) {
                stack.push_back(i);
                i = nodes[i].left;
            }
            i = stack.back();
            stack.pop_back();
            if (i != 1) {
                QueueDebugEntry<K, T> e;
                if constexpr (Indexed)
                    e.time = out.size();
                else
                    e.time = nodes[i].time;
                e.kind = nodes[i].operation;
                if (e.kind == qdkPush)
                    e.value = nodes[i].value;
                out.push_back(e);
            }
            i = nodes[i].right;
        }
        return out;
    }

    // 登録中の全操作と実際のpop結果を時刻順で返します。O(M log(M+2))時間・O(M)空間。
    // 添字版では現在の操作位置をtimeとして返します。
    auto debugTimeline() const {
        auto out = debugOperations();
        replayQueueDebug(out, order);
        return out;
    }

    std::string debugDump() const {
        return formatQueueDebug(debugTimeline());
    }
};
}
