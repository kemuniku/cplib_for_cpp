#pragma once
#include <cplib/graph/graph.hpp>
namespace cplib {
template<char First='a',char Last='z'> class AhoCorasick {
    static_assert(static_cast<unsigned char>(First)<=static_cast<unsigned char>(Last));
    static constexpr int alphabet=static_cast<unsigned char>(Last)-static_cast<unsigned char>(First)+1;
    struct Node{std::array<std::int32_t,alphabet> next{};std::int32_t parent=-1,failure=0,output=-1;char character=0;Int terminal=0,matched=0;};
    std::vector<Node> nodes;
    std::vector<std::int32_t> patterns;
    static int code(char c){return static_cast<unsigned char>(c)-static_cast<unsigned char>(First);}
    static bool valid(char c){return static_cast<unsigned char>(First)<=static_cast<unsigned char>(c)&&static_cast<unsigned char>(c)<=static_cast<unsigned char>(Last);}
public:
    AhoCorasick()=default;
    // 辞書とfailure遷移を構築する。O(Σ|s|+語数+σN)時間・O(σN+語数)領域。
    explicit AhoCorasick(std::span<const std::string> words):nodes(1){
        for(const auto& s:words){Int node=0;for(char c:s){assert(valid(c));if(!nodes[node].next[code(c)]){Int child=nodes.size();assert(child<=std::numeric_limits<std::int32_t>::max());nodes.emplace_back();nodes[child].parent=node;nodes[child].character=c;nodes[node].next[code(c)]=child;}node=nodes[node].next[code(c)];}++nodes[node].terminal;patterns.push_back(node);}
        nodes[0].matched=nodes[0].terminal;std::vector<Int> queue{0};
        for(std::size_t head=0;head<queue.size();++head){Int node=queue[head];for(int c=0;c<alphabet;++c){Int child=nodes[node].next[c];if(!child){if(node)nodes[node].next[c]=nodes[nodes[node].failure].next[c];continue;}Int failure=node?nodes[nodes[node].failure].next[c]:0;nodes[child].failure=failure;nodes[child].matched=nodes[child].terminal+nodes[failure].matched;nodes[child].output=nodes[failure].terminal>0?failure:nodes[failure].output;queue.push_back(child);}}
    }
    Int root()const{return 0;}
    Int nodeCount()const{return nodes.size();}
    Int patternNode(Int index)const{return patterns[index];}
    Int findNode(std::string_view s)const{if(nodes.empty())return -1;Int result=0;for(char c:s){if(!valid(c))return -1;Int child=nodes[result].next[code(c)];if(nodes[child].parent!=result)return -1;result=child;}return result;}
    Int next(Int node,char c)const{assert(node>=0&&node<Int(nodes.size()));return valid(c)?nodes[node].next[code(c)]:0;}
    Int next(Int node,std::string_view s)const{assert(node>=0&&node<Int(nodes.size()));for(char c:s)node=next(node,c);return node;}
    Int getParent(Int node)const{return nodes.at(node).parent;}
    Int failure(Int node)const{return nodes.at(node).failure;}
    Int terminal(Int node)const{return nodes.at(node).terminal;}
    Int matchCount(Int node)const{return nodes.at(node).matched;}
    struct MatchRange {
        const AhoCorasick* owner;Int first;
        struct iterator{const AhoCorasick* owner;Int id;Int operator*()const{return id;}iterator& operator++(){id=owner->nodes[id].output;return *this;}bool operator==(const iterator& other)const{return id==other.id;}};
        iterator begin()const{return {owner,first};}iterator end()const{return {owner,-1};}
    };
    // 一致語の節点を長い順に重複なく列挙する。O(列挙数+1)、追加領域O(1)。
    MatchRange matches(Int node)const{assert(node>=0&&node<Int(nodes.size()));return {this,nodes[node].terminal?node:nodes[node].output};}
    std::string restoreString(Int node)const{assert(node>=0&&node<Int(nodes.size()));std::string result;for(;node;node=nodes[node].parent)result+=nodes[node].character;std::reverse(result.begin(),result.end());return result;}
    auto toTrieGraph()const{assert(!nodes.empty());auto g=initWeightedDirectedGraph<char>(nodes.size(),nodes.size()-1);for(Int i=1;i<Int(nodes.size());++i)g.add_edge(nodes[i].parent,i,nodes[i].character);return g;}
    auto toFailureGraph()const{assert(!nodes.empty());auto g=initUnWeightedDirectedGraph(nodes.size(),nodes.size()-1);for(Int i=1;i<Int(nodes.size());++i)g.add_edge(i,nodes[i].failure);return g;}
};
template<char First='a',char Last='z'> class AhoCorasickPointer {
    AhoCorasick<First,Last>* owner=nullptr;Int id=0;
public:
    AhoCorasickPointer()=default;
    explicit AhoCorasickPointer(AhoCorasick<First,Last>& ac,Int node=0):owner(&ac),id(node){assert(node>=0&&node<ac.nodeCount());}
    Int nodeId()const{assert(owner);return id;}
    AhoCorasickPointer getParent()const{assert(owner&&id);return AhoCorasickPointer(*owner,owner->getParent(id));}
    AhoCorasickPointer failure()const{assert(owner);return AhoCorasickPointer(*owner,owner->failure(id));}
    template<class S> AhoCorasickPointer next(const S& s)const{assert(owner);return AhoCorasickPointer(*owner,owner->next(id,s));}
    template<class S> AhoCorasickPointer operator&(const S& s)const{return next(s);}
    template<class S> void add(const S& s){*this=next(s);}
    template<class S> AhoCorasickPointer& operator&=(const S& s){add(s);return *this;}
    Int terminal()const{assert(owner);return owner->terminal(id);}
    Int matchCount()const{assert(owner);return owner->matchCount(id);}
    auto matches()const{assert(owner);return owner->matches(id);}
    std::string restoreString()const{assert(owner);return owner->restoreString(id);}
};
template<char First='a',char Last='z'> auto initAhoCorasick(std::span<const std::string> words){return AhoCorasick<First,Last>(words);}
template<char First,char Last> auto initAhoCorasickPointer(AhoCorasick<First,Last>& ac,Int node=0){return AhoCorasickPointer<First,Last>(ac,node);}
template<char First,char Last> std::string to_string(const AhoCorasickPointer<First,Last>& p){return p.restoreString();}
}
