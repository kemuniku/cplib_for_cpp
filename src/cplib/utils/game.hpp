#pragma once
#include <cplib/common.hpp>
#include <functional>
#include <memory>
#include <stdexcept>
#include <unordered_map>
#include <unordered_set>
namespace cplib {
template<class T> using NextStates=std::function<std::vector<T>(T)>;
template<class T> using NextStatesByTurn=std::function<std::vector<T>(T,bool)>;
template<class T> struct OptimalPlay {bool is_win=false;std::vector<T> states;bool operator==(const OptimalPlay&)const=default;};
namespace detail {
template<class F> struct GameFirstArgument:GameFirstArgument<decltype(&F::operator())>{};
template<class R,class A,class... B> struct GameFirstArgument<R(*)(A,B...)>{using type=std::remove_cvref_t<A>;};
template<class C,class R,class A,class... B> struct GameFirstArgument<R(C::*)(A,B...)const>{using type=std::remove_cvref_t<A>;};
template<class C,class R,class A,class... B> struct GameFirstArgument<R(C::*)(A,B...)>{using type=std::remove_cvref_t<A>;};
template<class T,class Next> struct GameArgument {using type=T;};template<class Next> struct GameArgument<void,Next>:GameFirstArgument<Next>{};
template<class T> struct GameVisitGuard {std::unordered_set<T>& visiting;T state;~GameVisitGuard(){visiting.erase(state);}};
template<class T,class Next,int Mode> class GameSolver {
    struct Evaluation {bool winning=false;Int turns=0;T next_state;};
    using Value=std::conditional_t<Mode==0,Int,std::conditional_t<Mode==1,bool,Evaluation>>;
    struct Engine {
        Next next;bool terminal;std::function<bool(const T&,const T&)> prefer;
        std::array<std::unordered_map<T,Value>,2> memo;std::array<std::unordered_set<T>,2> visiting;
        Engine(Next n,bool t,std::function<bool(const T&,const T&)> p):next(std::move(n)),terminal(t),prefer(std::move(p)){}
        static constexpr bool by_turn=std::is_invocable_v<Next,const T&,bool>;
        Value solve(const T& state,bool first=true){std::size_t turn=by_turn?std::size_t(first):0;auto& cache=memo[turn];if(auto it=cache.find(state);it!=cache.end())return it->second;if(!visiting[turn].insert(state).second)throw std::invalid_argument("game graph contains a cycle");GameVisitGuard<T> guard{visiting[turn],state};auto next_states=[&](){if constexpr(by_turn)return next(state,first);else return next(state);}();
            if constexpr(Mode==0){std::vector<bool> seen(next_states.size()+1);for(const auto& s:next_states){Int value=solve(s,!first);if(value<Int(seen.size()))seen[value]=true;}Int value=0;while(seen[value])++value;cache[state]=value;return value;}
            else if constexpr(Mode==1){bool winning=next_states.empty()?terminal:false;for(const auto& s:next_states)if(!solve(s,!first)){winning=true;break;}cache[state]=winning;return winning;}
            else{Evaluation result{terminal,0,state};bool initial=true;for(const auto& s:next_states){auto child=solve(s,!first);bool winning=!child.winning;Int turns=child.turns+1;if(initial||(winning&&!result.winning)||(winning==result.winning&&((winning&&turns<result.turns)||(!winning&&turns>result.turns)||(turns==result.turns&&prefer&&prefer(s,result.next_state)))))result={winning,turns,s};initial=false;}cache.insert_or_assign(state,result);return result;}
        }
    };
    std::shared_ptr<Engine> engine_;
public:
    GameSolver(Next next,bool terminal=false,std::function<bool(const T&,const T&)> prefer={}):engine_(std::make_shared<Engine>(std::move(next),terminal,std::move(prefer))){}
    auto operator()(T state)const{if constexpr(Mode<2)return engine_->solve(state);else{OptimalPlay<T> result;result.is_win=engine_->solve(state).winning;result.states.push_back(state);bool first=true;for(;;){const auto& value=engine_->memo[Engine::by_turn?std::size_t(first):0].at(state);if(value.turns==0)break;state=value.next_state;first=!first;result.states.push_back(state);}return result;}}
};
}
// 同じ引数への遷移結果を固定して使用する。メモは返された関数のコピー間で共有。
// 新規状態 V、遷移 E に対して O(V+E)。最適手順の復元は別途 O(手順長)。
template<class T=void,class Next> auto init_grundy(Next next){using State=typename detail::GameArgument<T,Next>::type;return detail::GameSolver<State,Next,0>(std::move(next));}
template<class T,class Next> Int grundy(T state,Next next){return init_grundy<T>(std::move(next))(state);}
template<class T=void,class Next> auto init_can_win(Next next,bool terminal=false){using State=typename detail::GameArgument<T,Next>::type;return detail::GameSolver<State,Next,1>(std::move(next),terminal);}
template<class T,class Next> bool can_win(T state,Next next,bool terminal=false){return init_can_win<T>(std::move(next),terminal)(state);}
template<class T=void,class Next> auto init_optimal_play(Next next,bool terminal=false){using State=typename detail::GameArgument<T,Next>::type;return detail::GameSolver<State,Next,2>(std::move(next),terminal);}
template<class T=void,class Next,class Evaluate> requires (!std::is_same_v<std::remove_cvref_t<Evaluate>,bool>) auto init_optimal_play(Next next,Evaluate evaluate,bool terminal=false){using State=typename detail::GameArgument<T,Next>::type;return detail::GameSolver<State,Next,2>(std::move(next),terminal,[evaluate=std::move(evaluate)](const State& a,const State& b){return evaluate(b)<evaluate(a);});}
template<class T,class Next> auto optimal_play(T state,Next next,bool terminal=false){return init_optimal_play<T>(std::move(next),terminal)(state);}
template<class T,class Next,class Evaluate> requires (!std::is_same_v<std::remove_cvref_t<Evaluate>,bool>) auto optimal_play(T state,Next next,Evaluate evaluate,bool terminal=false){return init_optimal_play<T>(std::move(next),std::move(evaluate),terminal)(state);}
}
