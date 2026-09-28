#pragma once
#include <cplib/common.hpp>
namespace cplib {
class HopcroftKarp {
    std::vector<std::pair<Int,Int>> edges;
    std::vector<Int> leftOffset,rightOffset,leftEdges,rightEdges,leftMatch,rightMatch;
    Int size=0;bool built=false;
    void build(){if(built)return;Int left=leftMatch.size(),right=rightMatch.size(),m=edges.size();leftOffset.assign(left+1,0);rightOffset.assign(right+1,0);for(auto [v,u]:edges){++leftOffset[v];++rightOffset[u];}for(Int v=1;v<left;++v)leftOffset[v]+=leftOffset[v-1];for(Int u=1;u<right;++u)rightOffset[u]+=rightOffset[u-1];leftOffset[left]=m;rightOffset[right]=m;leftEdges.resize(m);rightEdges.resize(m);for(auto [v,u]:edges){leftEdges[--leftOffset[v]]=u;rightEdges[--rightOffset[u]]=v;}built=true;}
public:
    HopcroftKarp(Int left,Int right):leftMatch(left,-1),rightMatch(right,-1){assert(left>=0&&right>=0);}
    void add_edge(Int left,Int right){assert(0<=left&&left<Int(leftMatch.size())&&0<=right&&right<Int(rightMatch.size()));edges.emplace_back(left,right);built=false;}
    // 両方向CSR・global relabel前処理・非再帰HK。O((V+E)√V)、領域O(V+E)。
    Int matching(bool useRelabel=true){
        Int n=leftMatch.size(),right=rightMatch.size(),limit=std::min(n,right);if(size==limit)return size;build();std::vector<Int> dist(n),queue(n);
        if(useRelabel){std::vector<Int> active(right);Int head=0,tail=0,count=0;for(Int u=0;u<right;++u)if(rightMatch[u]<0&&rightOffset[u]<rightOffset[u+1]){active[tail++]=u;++count;}if(tail==right)tail=0;Int period=n+right,workLimit=16*(n+right+edges.size()),work=0,steps=0;
            while(count>0&&size<limit&&work<workLimit){if(steps==0){Int first=0,last=0;for(Int v=0;v<n;++v)if(leftMatch[v]<0){dist[v]=0;queue[last++]=v;}else dist[v]=n;work+=n;while(first<last){Int v=queue[first++];for(Int i=leftOffset[v];i<leftOffset[v+1];++i){Int w=rightMatch[leftEdges[i]];if(w>=0&&dist[w]==n){dist[w]=dist[v]+1;queue[last++]=w;}}work+=leftOffset[v+1]-leftOffset[v];}steps=period;if(work>=workLimit)break;}
                Int u=active[head++];if(head==right)head=0;--count;--steps;Int best=-1,bestHeight=n;for(Int i=rightOffset[u];i<rightOffset[u+1];++i){++work;Int v=rightEdges[i];if(dist[v]<bestHeight){best=v;bestHeight=dist[v];if(bestHeight==0)break;}}if(best<0)continue;Int previous=leftMatch[best];if(previous<0)++size;else{rightMatch[previous]=-1;active[tail++]=previous;if(tail==right)tail=0;++count;}leftMatch[best]=u;rightMatch[u]=best;dist[best]=bestHeight+1;
            }if(size==limit)return size;
        }else for(Int v=0;v<n;++v)if(leftMatch[v]<0)for(Int i=leftOffset[v];i<leftOffset[v+1];++i){Int u=leftEdges[i];if(rightMatch[u]<0){leftMatch[v]=u;rightMatch[u]=v;++size;break;}}
        std::vector<Int> iter(n),freeLeft;freeLeft.reserve(n-size);for(Int v=0;v<n;++v)if(leftMatch[v]<0&&leftOffset[v]<leftOffset[v+1])freeLeft.push_back(v);
        while(size<limit){Int head=0,tail=0;for(Int v=0;v<n;++v){iter[v]=leftOffset[v];dist[v]=-1;}for(Int v:freeLeft){dist[v]=0;queue[tail++]=v;}Int shortest=n;while(head<tail){Int v=queue[head++];if(dist[v]>=shortest)break;for(Int i=leftOffset[v];i<leftOffset[v+1];++i){Int u=leftEdges[i],w=rightMatch[u];if(w<0)shortest=dist[v];else if(dist[w]<0&&dist[v]<shortest){dist[w]=dist[v]+1;queue[tail++]=w;}}}if(shortest==n)break;
            for(Int root:freeLeft){if(leftMatch[root]>=0||dist[root]!=0)continue;Int depth=0;queue[0]=root;while(depth>=0){Int v=queue[depth];bool descended=false;Int endpoint=-1;while(iter[v]<leftOffset[v+1]){Int u=leftEdges[iter[v]++],w=rightMatch[u];if(w<0){if(dist[v]==shortest){endpoint=u;break;}}else if(dist[v]<shortest&&dist[w]==dist[v]+1){queue[++depth]=w;descended=true;break;}}
                    if(endpoint>=0){while(depth>=0){Int x=queue[depth],previous=leftMatch[x];leftMatch[x]=endpoint;rightMatch[endpoint]=x;endpoint=previous;dist[x]=-1;--depth;}++size;break;}if(!descended){dist[v]=-1;--depth;}
                }
            }Int remaining=0;for(Int v:freeLeft)if(leftMatch[v]<0)freeLeft[remaining++]=v;freeLeft.resize(remaining);
        }return size;
    }
    std::vector<std::pair<Int,Int>> get_matching()const{std::vector<std::pair<Int,Int>> result;result.reserve(size);for(Int l=0;l<Int(leftMatch.size());++l)if(leftMatch[l]>=0)result.emplace_back(l,leftMatch[l]);return result;}
};
inline auto initHopcroftKarp(Int left,Int right){return HopcroftKarp(left,right);}
inline void add_edge(HopcroftKarp& g,Int left,Int right){g.add_edge(left,right);}
inline Int matching(HopcroftKarp& g,bool useRelabel=true){return g.matching(useRelabel);}
inline auto get_matching(const HopcroftKarp& g){return g.get_matching();}
}
