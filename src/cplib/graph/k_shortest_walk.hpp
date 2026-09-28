#pragma once
#include <cplib/graph/graph.hpp>
#include <cplib/utils/constants.hpp>
#include <queue>
namespace cplib {
// Eppstein法。非負整数辺、同長の別ウォークと空ウォークを数える。
// 時間O((V+E)log V+k log k)、領域O(E+V log V+k)。不足分はinfinity。
template<GraphTypes G,std::signed_integral T> std::vector<T> k_shortest_walk(const G& graph,Int s,Int t,Int k,T infinity){
 Int n=graph.len;assert(0<=s&&s<n&&0<=t&&t<n&&k>=0);if(!k)return {};std::vector<T> result(k,infinity);
 struct ReverseEdge{Int vertex;T cost;Int edge;};std::vector<std::vector<ReverseEdge>> reverse(n);
 for(Int u=0;u<n;++u){Int edge=0;for(auto [v,cost]:graph.to_and_cost(u)){assert(cost>=0);reverse[v].push_back({u,T(cost),edge++});}}
 std::vector<T> distance(n);std::vector<bool> reached(n),settled(n);std::vector<Int> parent(n,-1),tree_edge(n,-1),position(n,-1),heap,order;
 // decrease-key付きヒープで辺数に対するlog Eを避ける。
 auto sift_up=[&](Int vertex){Int i=position[vertex];if(i<0){i=heap.size();heap.push_back(vertex);}while(i){Int p=(i-1)/2;if(distance[heap[p]]<=distance[vertex])break;heap[i]=heap[p];position[heap[i]]=i;i=p;}heap[i]=vertex;position[vertex]=i;};
 auto pop_vertex=[&](){Int result=heap[0],last=heap.back();heap.pop_back();position[result]=-1;if(heap.empty())return result;Int i=0;while(2*i+1<Int(heap.size())){Int child=2*i+1;if(child+1<Int(heap.size())&&distance[heap[child+1]]<distance[heap[child]])++child;if(distance[last]<=distance[heap[child]])break;heap[i]=heap[child];position[heap[i]]=i;i=child;}heap[i]=last;position[last]=i;return result;};
 reached[t]=true;sift_up(t);
 while(!heap.empty()){Int v=pop_vertex();settled[v]=true;order.push_back(v);for(auto edge:reverse[v]){Int u=edge.vertex;if(settled[u])continue;T candidate=distance[v]+edge.cost;if(!reached[u]||candidate<distance[u]){reached[u]=true;distance[u]=candidate;parent[u]=v;tree_edge[u]=edge.edge;sift_up(u);}}}
 if(!reached[s])return result;
 result[0]=distance[s];if(k==1)return result;
 struct Sidetrack{T delta;Int dest;};std::vector<std::vector<Sidetrack>> local(n);
 for(Int u:order){Int edge=0;for(auto [v,cost]:graph.to_and_cost(u)){if(reached[v]&&edge!=tree_edge[u])local[u].push_back({T(cost+distance[v]-distance[u]),v});++edge;}
  for(Int start=Int(local[u].size())/2-1;start>=0;--start){Int i=start;auto value=local[u][i];while(2*i+1<Int(local[u].size())){Int child=2*i+1;if(child+1<Int(local[u].size())&&local[u][child+1].delta<local[u][child].delta)++child;if(value.delta<=local[u][child].delta)break;local[u][i]=local[u][child];i=child;}local[u][i]=value;}
 }
 struct Node{Int vertex=0,left=0,right=0,rank=0;};std::vector<Node> nodes(1);std::vector<Int> roots(n);
 // 各頂点の最小非木辺を持つ永続leftist heap。
 auto meld=[&](auto&& self,Int a,Int b)->Int{if(!a)return b;if(!b)return a;if(local[nodes[b].vertex][0].delta<local[nodes[a].vertex][0].delta)std::swap(a,b);auto node=nodes[a];node.right=self(self,node.right,b);if(nodes[node.left].rank<nodes[node.right].rank)std::swap(node.left,node.right);node.rank=nodes[node.right].rank+1;nodes.push_back(node);return nodes.size()-1;};
 for(Int u:order){if(parent[u]>=0)roots[u]=roots[parent[u]];if(!local[u].empty()){Int id=nodes.size();nodes.push_back({u,0,0,1});roots[u]=meld(meld,roots[u],id);}}
 using Candidate=std::tuple<T,Int,Int,Int>;std::priority_queue<Candidate,std::vector<Candidate>,std::greater<Candidate>> queue;
 auto push_root=[&](T cost,Int root){if(root){Int u=nodes[root].vertex;queue.emplace(cost+local[u][0].delta,root,u,0);}};
 push_root(distance[s],roots[s]);Int count=1;
 while(!queue.empty()&&count<k){auto [cost,node,u,i]=queue.top();queue.pop();result[count++]=cost;if(count==k)break;auto edge=local[u][i];T base=cost-edge.delta;if(node){push_root(base,nodes[node].left);push_root(base,nodes[node].right);}for(Int child:{2*i+1,2*i+2})if(child<Int(local[u].size()))queue.emplace(base+local[u][child].delta,0,u,child);push_root(cost,roots[edge.dest]);}
 return result;
}
template<GraphTypes G> requires std::signed_integral<typename G::cost_type> auto k_shortest_walk(const G& graph,Int s,Int t,Int k){using T=typename G::cost_type;T infinity;if constexpr(sizeof(T)>=8)infinity=T(INF64);else if constexpr(sizeof(T)>=4)infinity=T(INF32);else infinity=std::numeric_limits<T>::max();return k_shortest_walk(graph,s,t,k,infinity);}
}
