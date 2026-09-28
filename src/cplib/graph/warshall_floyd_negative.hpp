#pragma once
#include <cplib/common.hpp>
namespace cplib {
// 負閉路検出後に影響を受ける組を-infにし、残りの距離を完成する。O(V^3)、領域O(V^2)。
template<class T,class Run> void warshall_floyd_negative_finish(std::vector<std::vector<T>>& d,T zero,T inf,Run run){
    Int n=d.size();std::vector<bool> seen(n);std::vector<Int> order;
    for(Int root=0;root<n;++root){
        if(seen[root])continue;
        std::vector<std::pair<Int,Int>> stack{{root,0}};seen[root]=true;
        while(!stack.empty()){
            auto [v,j]=stack.back();
            if(j==n){order.push_back(v);stack.pop_back();}
            else {++stack.back().second;if(d[v][j]!=inf && !seen[j]){seen[j]=true;stack.emplace_back(j,0);}}
        }
    }
    std::vector<Int> component(n,-1);std::vector<std::vector<Int>> groups;
    for(auto it=order.rbegin();it!=order.rend();++it){
        Int root=*it;if(component[root]!=-1)continue;Int id=groups.size();std::vector<Int> vertices{root};component[root]=id;
        for(std::size_t head=0;head<vertices.size();++head){Int v=vertices[head];for(Int u=0;u<n;++u)if(d[u][v]!=inf && component[u]==-1){component[u]=id;vertices.push_back(u);}}
        groups.push_back(std::move(vertices));
    }
    std::vector<std::vector<bool>> affected(n,std::vector<bool>(n));
    for(const auto& vertices:groups){
        bool negative=groups.size()==1;
        for(Int v:vertices)if(d[v][v]<zero)negative=true;
        if(!negative && vertices.size()>1){
            std::vector<std::vector<T>> block(vertices.size(),std::vector<T>(vertices.size()));
            for(std::size_t i=0;i<vertices.size();++i)for(std::size_t j=0;j<vertices.size();++j)block[i][j]=d[vertices[i]][vertices[j]];
            negative=run(block,zero,inf);
        }
        if(!negative)continue;
        std::vector<bool> before(n),after(n);
        for(bool reverse:{false,true}){
            std::vector<bool> reached(n);std::vector<Int> queue{vertices[0]};reached[vertices[0]]=true;
            for(std::size_t head=0;head<queue.size();++head){Int v=queue[head];for(Int u=0;u<n;++u){T edge=reverse?d[u][v]:d[v][u];if(edge!=inf && !reached[u]){reached[u]=true;queue.push_back(u);}}}
            if(reverse)before=std::move(reached);else after=std::move(reached);
        }
        for(Int i=0;i<n;++i)if(before[i])for(Int j=0;j<n;++j)if(after[j])affected[i][j]=true;
    }
    for(Int i=0;i<n;++i)for(Int j=0;j<n;++j)if(affected[i][j])d[i][j]=inf;
    run(d,zero,inf);
    for(Int i=0;i<n;++i)for(Int j=0;j<n;++j)if(affected[i][j])d[i][j]=-inf;
}
}
