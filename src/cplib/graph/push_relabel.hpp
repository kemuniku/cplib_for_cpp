#pragma once
#include <cplib/graph/private/flow_graph.hpp>
namespace cplib {
template<class Cap> using PushRelabelEdge=FlowEdge<Cap>;
template<class Cap=Int> class PushRelabel:public detail::FlowGraph<Cap> {
    using Base=detail::FlowGraph<Cap>;using Base::graph;using Base::positions;
public:
    explicit PushRelabel(Int n):Base(n){}
    // Highest-label、global relabel、gap heuristicを保持。一般O(VE+V²√E)、追加領域O(V)。
    Cap flow(Int src,Int dst,Cap limit=std::numeric_limits<Cap>::max()){
        assert(0<=src&&src<Int(graph.size())&&0<=dst&&dst<Int(graph.size())&&src!=dst&&limit>=Cap(0));if(limit==Cap(0))return 0;
        Int source=graph.size(),sourceIndex=graph[src].size();graph.push_back({{src,sourceIndex,Cap(0)}});graph[src].push_back({source,0,limit});
        // 仮想始点を後で除く。例外が発生しても残余グラフの形を復元する。
        struct Restore{decltype(graph)& g;Int src,source,index;~Restore(){g[src].resize(index);g.resize(source);}} restore{graph,src,source,sourceIndex};
        Int n=graph.size(),unreachable=2*n;std::vector<Int> height(n),count(unreachable+1),iter(n),bucket(unreachable+1),next(n),prev(n),heightHead(n),heightNext(n),heightPrev(n),queue(n);std::vector<Cap> excess(n);std::vector<bool> active(n);Int highest=-1,highestLow=-1,highestFinite=-1,relabelWork=0,relabelWorkLimit=4*n+2*positions.size()+2;excess[src]=limit;
        auto activate=[&](Int v){if(v!=source&&v!=dst){Int h=height[v];next[v]=bucket[h];prev[v]=-1;if(next[v]>=0)prev[next[v]]=v;bucket[h]=v;active[v]=true;highest=std::max(highest,h);if(h<n)highestLow=std::max(highestLow,h);}};
        auto removeActive=[&](Int v){if(prev[v]<0)bucket[height[v]]=next[v];else next[prev[v]]=next[v];if(next[v]>=0)prev[next[v]]=prev[v];active[v]=false;};
        auto addHeight=[&](Int v){Int h=height[v];if(h<n){heightNext[v]=heightHead[h];heightPrev[v]=-1;if(heightNext[v]>=0)heightPrev[heightNext[v]]=v;heightHead[h]=v;highestFinite=std::max(highestFinite,h);}};
        auto removeHeight=[&](Int v){if(height[v]<n){if(heightPrev[v]<0)heightHead[height[v]]=heightNext[v];else heightNext[heightPrev[v]]=heightNext[v];if(heightNext[v]>=0)heightPrev[heightNext[v]]=heightPrev[v];}};
        auto rebuildActive=[&](){std::fill(bucket.begin(),bucket.end(),-1);highest=-1;highestLow=-1;for(Int v=0;v<n;++v){active[v]=false;if(excess[v]>Cap(0))activate(v);}};
        auto globalRelabel=[&](){std::fill(height.begin(),height.end(),unreachable);std::fill(iter.begin(),iter.end(),0);height[dst]=0;height[source]=n;for(Int root:{dst,source}){Int head=0,tail=1;queue[0]=root;while(head<tail){Int v=queue[head++];for(auto e:graph[v])if(height[e.dst]==unreachable&&graph[e.dst][e.rev].cap>Cap(0)){height[e.dst]=height[v]+1;queue[tail++]=e.dst;}}}std::fill(count.begin(),count.end(),0);std::fill(heightHead.begin(),heightHead.end(),-1);highestFinite=-1;for(Int v=0;v<n;++v){++count[height[v]];addHeight(v);}rebuildActive();relabelWork=0;};
        globalRelabel();while(highest>=0){if(bucket[highest]<0){if(highest==n)highest=highestLow;else{--highest;if(highest<n)highestLow=highest;}continue;}Int v=bucket[highest];removeActive(v);while(excess[v]>Cap(0)){
            if(iter[v]==Int(graph[v].size())){Int old=height[v],best=unreachable,bestIndex=0;for(Int i=0;i<Int(graph[v].size());++i){auto e=graph[v][i];if(e.cap>Cap(0)&&height[e.dst]<best){best=height[e.dst];bestIndex=i;}}relabelWork+=graph[v].size();removeHeight(v);--count[old];height[v]=best+1;++count[height[v]];addHeight(v);iter[v]=bestIndex;
                if(old<n&&count[old]==0){for(Int h=old+1;h<=highestFinite;++h)while(heightHead[h]>=0){Int u=heightHead[h];removeHeight(u);bool wasActive=active[u];if(wasActive)removeActive(u);--count[h];height[u]=n+1;++count[height[u]];iter[u]=0;if(wasActive)activate(u);}highestFinite=old;activate(v);break;}
            }else{Int i=iter[v];auto e=graph[v][i];if(e.cap>Cap(0)&&height[v]==height[e.dst]+1){Cap amount=std::min(excess[v],e.cap);if(excess[e.dst]==Cap(0))activate(e.dst);graph[v][i].cap-=amount;graph[e.dst][e.rev].cap+=amount;excess[v]-=amount;excess[e.dst]+=amount;}else ++iter[v];}
        }if(relabelWork>=relabelWorkLimit)globalRelabel();}return excess[dst];
    }
};
template<class Cap=Int> auto initPushRelabel(Int n,Cap=Cap(0)){return PushRelabel<Cap>(n);}
template<class Cap> Int add_edge(PushRelabel<Cap>& g,Int s,Int t,Cap c){return g.add_edge(s,t,c);}
template<class Cap> auto get_edge(const PushRelabel<Cap>& g,Int i){return g.get_edge(i);}
template<class Cap> auto get_edges(const PushRelabel<Cap>& g){return g.get_edges();}
template<class Cap> Cap flow(PushRelabel<Cap>& g,Int s,Int t,Cap limit=std::numeric_limits<Cap>::max()){return g.flow(s,t,limit);}
template<class Cap> auto min_cut(const PushRelabel<Cap>& g,Int s){return g.min_cut(s);}
}
