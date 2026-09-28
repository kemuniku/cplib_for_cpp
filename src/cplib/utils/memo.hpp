#pragma once
#include <cplib/common.hpp>
#include <functional>
#include <memory>
#include <tuple>
#include <unordered_map>
namespace cplib {
namespace detail {
struct MemoTupleHash{template<class... Args>std::size_t operator()(const std::tuple<Args...>& key)const{std::size_t h=0;std::apply([&](const auto&... xs){((h^=std::hash<std::decay_t<decltype(xs)>>{}(xs)+0x9e3779b97f4a7c15ULL+(h<<6)+(h>>2)),...);},key);return h;}};
template<class Signature,class Body> class Memoized;
template<class Result,class... Args,class Body> class Memoized<Result(Args...),Body> {
    static_assert(sizeof...(Args)>0);using Key=std::tuple<std::decay_t<Args>...>;
    struct State{Body body;std::unordered_map<Key,Result,MemoTupleHash> cache;explicit State(Body f):body(std::move(f)){}};
    std::shared_ptr<State> state;
public:
    explicit Memoized(Body body):state(std::make_shared<State>(std::move(body))){}
    Result operator()(Args... args)const{Key key(args...);auto it=state->cache.find(key);if(it!=state->cache.end())return it->second;Result value;if constexpr(std::is_invocable_v<Body,const Memoized&,Args...>)value=state->body(*this,args...);else value=state->body(args...);state->cache.emplace(std::move(key),value);return value;}
};
}
// memo<R(Args...)>([](auto& self,Args...){...})。再帰はself(...)で呼び、版間でメモを共用。
template<class Signature,class Body> auto memo(Body body){return detail::Memoized<Signature,Body>(std::move(body));}
}
