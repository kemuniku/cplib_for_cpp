#pragma once
#include <cplib/tree/heavylightdecomposition.hpp>
#include <type_traits>
#include <cplib/utils/auto_rollback.hpp>

namespace cplib {
// 各頂点を含む部分木を集計してanswer(v)を呼ぶ。clear(v)は寄与を削除し、clear()なら一括削除する。
// 空の状態から開始し、終了・例外時は空に戻す。addは例外時に状態を変更せず、clearは例外を送出しないこと。
// 各コールバックがO(1)なら時間 O(N log N)、追加領域 O(1)。
template <class Add, class Answer, class Clear>
void dsuOnTree(const HeavyLightDecomposition &hld, Add add, Answer answer, Clear clear) {
    assert(hld.numVertices() > 0);
    Int activeLeft = 0, activeRight = 0;
    auto clearActive = [&]() {
        Int first = activeLeft, last = activeRight;
        activeLeft = activeRight;
        if (first == last)
            return;
        if constexpr (std::is_invocable_v<Clear>)
            clear();
        else
            for (Int j = first; j < last; ++j)
                clear(hld.toVtx(j));
    };
    try {
        for (Int i = hld.numVertices(); i-- > 0;) {
            Int v = hld.toVtx(i), heavy = hld.heavyChildOf(v);
            Int first = heavy < 0 ? i + 1 : hld.subtree(heavy).second;
            if (activeLeft == activeRight)
                activeLeft = activeRight = first;
            for (Int j = first; j < hld.subtree(v).second; ++j) {
                add(hld.toVtx(j));
                activeRight = j + 1;
            }
            add(v);
            activeLeft = i;
            answer(v);
            if (hld.heavyRootOf(v) == v)
                clearActive();
        }
    } catch (...) {
        clearActive();
        throw;
    }
    clearActive();
}
// 木を根付けてDSU on treeを実行する。時間 O(N log N)、HLD用領域 O(N)。
template <class G, class Add, class Answer, class Clear>
    requires (!std::is_same_v<G, HeavyLightDecomposition>)
void dsuOnTree(const G &g, Add add, Answer answer, Clear clear, Int root = 0) {
    auto hld = initHld(g, root);
    dsuOnTree(hld, add, answer, clear);
}
}

namespace cplib {
// 履歴に記録した更新を自動復元する。add(log,v)でlog.write/setを通して変更すること。
// answerの変更は履歴に記録しない。O(N log N)回の追加と記録・復元の費用、履歴領域は最大保存量に比例する。
template <class Add, class Answer>
    requires std::is_invocable_v<Add, RollbackLog &, Int>
void dsuOnTree(const HeavyLightDecomposition &hld, Add add, Answer answer) {
    RollbackLog log;
    dsuOnTree(hld, [&](Int v) { add(log, v); }, answer, [&]() { log.rollback(0); });
}
// 木を根付け、履歴付きadd(log,v)で部分木を集計する。時間 O(N log N)と記録・復元費用、HLD用領域 O(N)。
template <class G, class Add, class Answer>
    requires (!std::is_same_v<G, HeavyLightDecomposition> && std::is_invocable_v<Add, RollbackLog &, Int>)
void dsuOnTree(const G &g, Add add, Answer answer, Int root = 0) {
    auto hld = initHld(g, root);
    dsuOnTree(hld, add, answer);
}
}
