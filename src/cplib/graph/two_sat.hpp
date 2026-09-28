#pragma once
#include <cplib/common.hpp>
#include <memory>
#include <stdexcept>
namespace cplib {
namespace detail {
struct TwoSatState {
    std::vector<std::pair<Int,Int>> clauses;
    bool evaluated=false,solved=false;
    std::vector<bool> assignment;
};
}
struct Constraint2sat {
    std::shared_ptr<detail::TwoSatState> problem;
    std::array<std::pair<Int,Int>,2> clauses{};
    Int count=0;
};
struct Literal2sat {
    std::shared_ptr<detail::TwoSatState> problem;
    Int vertex=0;
    // 解におけるリテラルの値を取得する。O(1)。
    bool get()const{if(!problem || !problem->solved)throw std::invalid_argument("solve must succeed first");bool v=problem->assignment[vertex/2];return vertex&1?v:!v;}
    std::string to_string()const{if(!problem || !problem->solved)return "-";return get()?"1":"0";}
};
inline Literal2sat toLiteral(Literal2sat a){return a;}
inline Literal2sat operator!(Literal2sat a){a.vertex^=1;return a;}
// 2リテラルのOR制約を作る。O(1)。
inline Constraint2sat operator||(Literal2sat a,Literal2sat b){if(!a.problem || a.problem!=b.problem)throw std::invalid_argument("literals must belong to the same problem");Constraint2sat r;r.problem=a.problem;r.clauses[0]={a.vertex,b.vertex};r.count=1;return r;}
inline Constraint2sat operator|(Literal2sat a,Literal2sat b){return a||b;}
inline Constraint2sat nand(Literal2sat a,Literal2sat b){return !a||!b;}
inline Constraint2sat implies(Literal2sat a,Literal2sat b){return !a||b;}
inline Constraint2sat operator^(Literal2sat a,Literal2sat b){auto r=a||b;r.clauses[1]={a.vertex^1,b.vertex^1};r.count=2;return r;}
inline Constraint2sat operator==(Literal2sat a,Literal2sat b){return a^!b;}
inline Constraint2sat operator!=(Literal2sat a,Literal2sat b){return a^b;}
inline Constraint2sat operator==(Literal2sat a,bool b){return b?(a||a):(!a||!a);}
inline Constraint2sat operator==(bool a,Literal2sat b){return b==a;}
inline Constraint2sat operator!=(Literal2sat a,bool b){return a==!b;}
inline Constraint2sat operator!=(bool a,Literal2sat b){return b==!a;}
class Problem2sat {
    std::shared_ptr<detail::TwoSatState> p;
public:
    Problem2sat()=default;
    // n変数の問題を構築する。コピーは元実装と同じ問題を共有する。O(n)。
    explicit Problem2sat(Int n){if(n<0)throw std::invalid_argument("negative variable count");p=std::make_shared<detail::TwoSatState>();p->assignment.resize(n);}
    Literal2sat operator[](Int k)const{if(!p)throw std::invalid_argument("uninitialized problem");if(k<0 || k>=Int(p->assignment.size()))throw std::out_of_range("variable index");return {p,2*k+1};}
    Problem2sat& operator+=(const Constraint2sat& c){if(!p || c.problem!=p || !c.count)throw std::invalid_argument("constraint belongs to another problem");for(Int i=0;i<c.count;++i)p->clauses.push_back(c.clauses[i]);p->solved=p->evaluated=false;return *this;}
    void add_clause(Int i,bool f,Int j,bool g){auto a=(*this)[i],b=(*this)[j];*this+=(f?a:!a)||(g?b:!b);}
    // 含意グラフのCSRと対偶を用いて求解する。O(n+m)、再実行O(1)。
    bool solve(){
        if(!p)throw std::invalid_argument("uninitialized problem");
        if(p->evaluated)return p->solved;
        p->solved=false;Int n=2*p->assignment.size();std::vector<Int> offsets(n+1);
        for(auto [a,b]:p->clauses){++offsets[(a^1)+1];++offsets[(b^1)+1];}
        for(Int v=0;v<n;++v)offsets[v+1]+=offsets[v];
        auto cursor=offsets;std::vector<Int> edges(2*p->clauses.size());
        for(auto [a,b]:p->clauses){edges[cursor[a^1]++]=b;edges[cursor[b^1]++]=a;}
        cursor=offsets;std::vector<bool> used(n);std::vector<Int> order,stack;order.reserve(n);stack.reserve(n);
        for(Int root=0;root<n;++root){
            if(used[root])continue;
            used[root]=true;stack.push_back(root);
            while(!stack.empty()){Int v=stack.back();if(cursor[v]==offsets[v+1]){order.push_back(v);stack.pop_back();}else{Int dst=edges[cursor[v]++];if(!used[dst]){used[dst]=true;stack.push_back(dst);}}}
        }
        std::vector<Int> component(n,-1);Int id=0;
        for(auto it=order.rbegin();it!=order.rend();++it){
            Int root=*it;if(component[root]!=-1)continue;component[root]=id;stack.push_back(root);
            while(!stack.empty()){Int v=stack.back();stack.pop_back();Int opposite=v^1;for(Int e=offsets[opposite];e<offsets[opposite+1];++e){Int dst=edges[e]^1;if(component[dst]==-1){component[dst]=id;stack.push_back(dst);}}}
            ++id;
        }
        p->evaluated=true;
        for(Int i=0;i<Int(p->assignment.size());++i){if(component[2*i]==component[2*i+1])return false;p->assignment[i]=component[2*i]<component[2*i+1];}
        p->solved=true;return true;
    }
    bool satisfiable(){return solve();}
    std::vector<bool> answer()const{if(!p || !p->solved)throw std::invalid_argument("solve must succeed first");return p->assignment;}
};
inline Problem2sat initTwoSat(Int n){return Problem2sat(n);}
inline void add_clause(Problem2sat& p,Int i,bool f,Int j,bool g){p.add_clause(i,f,j,g);}
inline bool solve(Problem2sat& p){return p.solve();}
inline bool satisfiable(Problem2sat& p){return p.satisfiable();}
inline auto answer(const Problem2sat& p){return p.answer();}
inline bool get(const Literal2sat& a){return a.get();}
inline std::string to_string(const Literal2sat& a){return a.to_string();}
}
