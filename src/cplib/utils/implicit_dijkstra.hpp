#pragma once
#include <cplib/utils/constants.hpp>
#include <functional>
#include <queue>
#include <unordered_map>
namespace cplib {
template<class T> using ImplicitDijkstraAdjacent=std::function<std::vector<std::pair<T,Int>>(T)>;
template<class T> using ImplicitDijkstraFinish=std::function<bool(T)>;
template<class T> struct ImplicitDijkstraResult {std::unordered_map<T,Int> costs;std::unordered_map<T,T> prev;};
template<class T> struct ImplicitDijkstraPath {std::vector<T> path;Int cost=0;};
namespace detail {
// ヒープには発見順の整数 ID を置き、頂点型に大小比較を要求しない。
template<bool Restore,class T,class Adjacent,class Finish> auto implicit_dijkstra_run(T start,Adjacent adjacent,Finish finish,Int inf){
    ImplicitDijkstraResult<T> result;std::unordered_map<T,Int> ids;std::vector<T> vertices{start};std::priority_queue<std::pair<Int,Int>,std::vector<std::pair<Int,Int>>,std::greater<std::pair<Int,Int>>> queue;ids[start]=0;result.costs[start]=0;queue.emplace(0,0);Int answer=inf;
    while(!queue.empty()){auto [cost,id]=queue.top();queue.pop();T vertex=vertices[id];if(cost!=result.costs.at(vertex))continue;if(finish(vertex)){answer=cost;break;}for(auto [next,edgeCost]:adjacent(vertex)){Int nextCost=cost+edgeCost;auto it=result.costs.find(next);if(it==result.costs.end()||nextCost<it->second){if(!ids.contains(next)){ids[next]=vertices.size();vertices.push_back(next);}result.costs[next]=nextCost;if constexpr(Restore)result.prev[next]=vertex;queue.emplace(nextCost,ids.at(next));}}}
    return std::pair{std::move(result),answer};
}
}
// 到達可能な状態を全探索。非負辺を仮定し O((V+E) log(V+E))。
template<class T,class Adjacent> auto restore_implicit_dijkstra(T start,Adjacent adjacent){return detail::implicit_dijkstra_run<true>(start,adjacent,[](const T&){return false;},INF64).first;}
template<class T,class Adjacent> auto implicit_dijkstra(T start,Adjacent adjacent){return restore_implicit_dijkstra(start,adjacent).costs;}
template<class T,class Adjacent,class Finish> Int implicit_dijkstra_until(T start,Adjacent adjacent,Finish finish,Int inf=INF64){return detail::implicit_dijkstra_run<false>(start,adjacent,finish,inf).second;}
template<class T,class Adjacent> auto shortest_path_implicit_dijkstra(T start,T goal,Adjacent adjacent,Int inf=INF64){bool reached=false;auto [result,cost]=detail::implicit_dijkstra_run<true>(start,adjacent,[&](const T& v){return reached=(v==goal);},inf);ImplicitDijkstraPath<T> out{{},cost};if(reached){T current=goal;out.path.push_back(current);while(!(current==start)){current=result.prev.at(current);out.path.push_back(current);}std::reverse(out.path.begin(),out.path.end());}return out;}
}
