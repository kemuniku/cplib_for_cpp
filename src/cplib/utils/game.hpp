#pragma once
#include <cplib/common.hpp>
#include <functional>
#include <memory>
#include <stdexcept>
#include <unordered_map>
#include <unordered_set>

namespace cplib {
template <class T> using NextStates = std::function<std::vector<T>(T)>;
template <class T> using NextStatesByTurn = std::function<std::vector<T>(T, bool)>;

template <class T> struct OptimalPlay {
    bool is_win = false;
    std::vector<T> states;
    bool operator==(const OptimalPlay &) const = default;
};

namespace detail {
template <class F> struct GameFirstArgument : GameFirstArgument<decltype(&F::operator())> {};

template <class R, class A, class... B> struct GameFirstArgument<R (*)(A, B...)> {
    using type = std::remove_cvref_t<A>;
};

template <class C, class R, class A, class... B> struct GameFirstArgument<R (C::*)(A, B...) const> {
    using type = std::remove_cvref_t<A>;
};

template <class C, class R, class A, class... B> struct GameFirstArgument<R (C::*)(A, B...)> {
    using type = std::remove_cvref_t<A>;
};

template <class T, class Next> struct GameArgument {
    using type = T;
};

template <class Next> struct GameArgument<void, Next> : GameFirstArgument<Next> {};

template <class T> struct GameVisitGuard {
    std::unordered_set<T> &visiting;
    T state;

    ~GameVisitGuard() {
        visiting.erase(state);
    }
};

template <class T, class Next, int Mode> class GameSolver {
    struct Evaluation {
        bool winning = false;
        Int turns = 0;
        T next_state;
    };

    using Value =
        std::conditional_t<Mode == 0, Int, std::conditional_t<Mode == 1, bool, Evaluation>>;

    struct Engine {
        Next next;
        bool terminal;
        std::function<bool(const T &, const T &)> prefer;
        std::array<std::unordered_map<T, Value>, 2> memo;
        std::array<std::unordered_set<T>, 2> visiting;

        Engine(Next n, bool t, std::function<bool(const T &, const T &)> p)
            : next(std::move(n)), terminal(t), prefer(std::move(p)) {
        }

        static constexpr bool by_turn = std::is_invocable_v<Next, const T &, bool>;

        Value solve(const T &state, bool first = true) {
            std::size_t turn = by_turn ? std::size_t(first) : 0;
            auto &cache = memo[turn];
            if (auto it = cache.find(state); it != cache.end())
                return it->second;
            if (!visiting[turn].insert(state).second)
                throw std::invalid_argument("game graph contains a cycle");
            GameVisitGuard<T> guard{visiting[turn], state};
            auto next_states = [&]() {
                if constexpr (by_turn)
                    return next(state, first);
                else
                    return next(state);
            }();
            if constexpr (Mode == 0) {
                std::vector<bool> seen(next_states.size() + 1);
                for (const auto &s : next_states) {
                    Int value = solve(s, !first);
                    if (value < Int(seen.size()))
                        seen[value] = true;
                }
                Int value = 0;
                while (seen[value])
                    ++value;
                cache[state] = value;
                return value;
            } else if constexpr (Mode == 1) {
                bool winning = next_states.empty() ? terminal : false;
                for (const auto &s : next_states)
                    if (!solve(s, !first)) {
                        winning = true;
                        break;
                    }
                cache[state] = winning;
                return winning;
            } else {
                Evaluation result{terminal, 0, state};
                bool initial = true;
                for (const auto &s : next_states) {
                    auto child = solve(s, !first);
                    bool winning = !child.winning;
                    Int turns = child.turns + 1;
                    if (initial || (winning && !result.winning) ||
                        (winning == result.winning &&
                         ((winning && turns < result.turns) || (!winning && turns > result.turns) ||
                          (turns == result.turns && prefer && prefer(s, result.next_state)))))
                        result = {winning, turns, s};
                    initial = false;
                }
                cache.insert_or_assign(state, result);
                return result;
            }
        }
    };

    std::shared_ptr<Engine> engine_;

public:
    GameSolver(Next next, bool terminal = false,
               std::function<bool(const T &, const T &)> prefer = {})
        : engine_(std::make_shared<Engine>(std::move(next), terminal, std::move(prefer))) {
    }

    auto operator()(T state) const {
        if constexpr (Mode < 2)
            return engine_->solve(state);
        else {
            OptimalPlay<T> result;
            result.is_win = engine_->solve(state).winning;
            result.states.push_back(state);
            bool first = true;
            for (;;) {
                const auto &value =
                    engine_->memo[Engine::by_turn ? std::size_t(first) : 0].at(state);
                if (value.turns == 0)
                    break;
                state = value.next_state;
                first = !first;
                result.states.push_back(state);
            }
            return result;
        }
    }
};
}

// 同じ引数への遷移結果を固定して使用する。メモは返された関数のコピー間で共有。
// 新規状態 V、遷移 E に対して O(V+E)。最適手順の復元は別途 O(手順長)。
// 両者の合法手が同じで、合法手がない側が負ける有限かつ閉路なしのゲームを対象とする。
// 遷移先のGrundy数のmexを返し、終端状態では0。値0は手番側の負け、独立なゲームの和はxorで求める。
// 返された関数を使う間、nextの結果は同じ引数に対して変化しないこと。ルール変更時は関数を作り直す。
// 状態のハッシュ・比較・コピーはO(1)、nextは列挙数に比例する時間とする。
template <class T = void, class Next> auto init_grundy(Next next) {
    using State = typename detail::GameArgument<T, Next>::type;
    return detail::GameSolver<State, Next, 0>(std::move(next));
}

// 指定状態のGrundy数を、呼び出しごとに新しいメモで計算する。ゲーム条件・計算量はinit_grundyと同じ。
template <class T, class Next> Int grundy(T state, Next next) {
    return init_grundy<T>(std::move(next))(state);
}

// 手番側が必勝かを判定し、呼び出し間でメモを再利用する関数を返す。
// next(state)は一手の遷移先を返す。合法手がない状態ではterminalがtrueの場合に限り手番側の勝ち。
// ゲームグラフは有限かつ閉路なし。関数を使う間、nextの結果は同じ引数に対して固定する。
// next(state, is_first)の形では最初の手番がtrueで、一手ごとに反転し、状態と手番の組でメモ化する。
template <class T = void, class Next> auto init_can_win(Next next, bool terminal = false) {
    using State = typename detail::GameArgument<T, Next>::type;
    return detail::GameSolver<State, Next, 1>(std::move(next), terminal);
}

// 指定状態の勝敗を、呼び出しごとに新しいメモで判定する。ゲームグラフは有限かつ閉路なし。
template <class T, class Next> bool can_win(T state, Next next, bool terminal = false) {
    return init_can_win<T>(std::move(next), terminal)(state);
}

// 勝敗と最適な状態列を返し、呼び出し間でメモを再利用する関数を作る。
// 返り値のis_winは初期状態の手番側が必勝かを表し、statesは初期状態と終端状態を含む。
// 勝つ手数を最小化し、負ける手数を最大化する。勝敗と手数が同じならnextの順で選ぶ。
// 合法手がない状態はterminalがtrueの場合に限り勝ち。有限かつ閉路なしのゲームで、nextの結果は固定する。
// 手番別next(state, is_first)では各呼び出しを先手から始め、statesの偶数番目がtrue、奇数番目がfalse。
// 新規状態V、遷移E、復元列長Lに対してO(V+E+L)。ハッシュ・比較・コピーはO(1)、nextは列挙数に比例する時間とする。
template <class T = void, class Next> auto init_optimal_play(Next next, bool terminal = false) {
    using State = typename detail::GameArgument<T, Next>::type;
    return detail::GameSolver<State, Next, 2>(std::move(next), terminal);
}

// 勝敗と終局手数が同じ候補は、遷移先のevaluate(state)が大きい手を優先してメモ化する。
// 両プレイヤーとも評価値を最大化し、同値ならnextの順。評価値は状態列全体の合計ではない。
// 評価値は<で比較でき、関数を使う間、同じ引数へのevaluateの結果は固定する。
// その他の仕様は評価関数なしのinit_optimal_playと同じ。evaluateと比較がO(1)なら計算量も同じ。
template <class T = void, class Next, class Evaluate>
    requires(!std::is_same_v<std::remove_cvref_t<Evaluate>, bool>)
auto init_optimal_play(Next next, Evaluate evaluate, bool terminal = false) {
    using State = typename detail::GameArgument<T, Next>::type;
    return detail::GameSolver<State, Next, 2>(
        std::move(next), terminal,
        [evaluate = std::move(evaluate)](const State &a, const State &b) {
            return evaluate(b) < evaluate(a);
        });
}

// 呼び出しごとに新しいメモで、勝敗と初期状態から終端状態までの状態列を返す。
// 選択規則・前提・計算量はinit_optimal_playと同じで、nextの結果は呼び出し中一定。
template <class T, class Next> auto optimal_play(T state, Next next, bool terminal = false) {
    return init_optimal_play<T>(std::move(next), terminal)(state);
}

// 評価関数付きinit_optimal_playと同じ規則で、呼び出しごとに新しいメモを使って復元する。
// nextとevaluateの結果は呼び出し中一定。
template <class T, class Next, class Evaluate>
    requires(!std::is_same_v<std::remove_cvref_t<Evaluate>, bool>)
auto optimal_play(T state, Next next, Evaluate evaluate, bool terminal = false) {
    return init_optimal_play<T>(std::move(next), std::move(evaluate), terminal)(state);
}
}
