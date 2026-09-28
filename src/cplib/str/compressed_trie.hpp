#pragma once
#include <cplib/str/static_string.hpp>
#include <cplib/graph/graph.hpp>
namespace cplib {
struct CompressedTrieNode {
    std::weak_ptr<CompressedTrieNode> parent;std::array<std::shared_ptr<CompressedTrieNode>,26> child{};StaticString<char> s;std::int32_t cnt=0,subtree_sum=0;
};
using CompressedTriePtr=std::shared_ptr<CompressedTrieNode>;
inline CompressedTriePtr initCompressedTrie(std::vector<StaticString<char>> strings,bool sorted=false){if(!sorted)std::sort(strings.begin(),strings.end());auto root=std::make_shared<CompressedTrieNode>();if(strings.empty())return root;root->s=strings[0].substr(0,0);std::vector<CompressedTriePtr> stack{root};auto slot=[](char c){assert('a'<=c&&c<='z');return c-'a';};for(const auto& s:strings){while(!startsWith(s,stack.back()->s))stack.pop_back();Int l=lcp(stack.back()->s,s);if(l==s.len()){++stack.back()->cnt;continue;}auto parent=stack.back();auto& child=parent->child[slot(s[l])];if(!child){child=std::make_shared<CompressedTrieNode>();child->parent=parent;child->s=s;child->cnt=child->subtree_sum=1;stack.push_back(child);}else{Int x=lcp(s,child->s);auto tmp=std::make_shared<CompressedTrieNode>();tmp->parent=parent;tmp->s=s.substr(0,x);child->parent=tmp;tmp->child[slot(child->s[x])]=child;auto leaf=std::make_shared<CompressedTrieNode>();leaf->parent=tmp;leaf->s=s;leaf->cnt=leaf->subtree_sum=1;tmp->child[slot(s[x])]=leaf;child=tmp;stack.push_back(tmp);stack.push_back(leaf);}}
    std::vector<std::pair<CompressedTriePtr,bool>> dfs{{root,false}};while(!dfs.empty()){auto [node,visited]=dfs.back();dfs.pop_back();if(visited){node->subtree_sum=node->cnt;for(const auto& c:node->child)if(c)node->subtree_sum+=c->subtree_sum;}else{dfs.emplace_back(node,true);for(const auto& c:node->child)if(c)dfs.emplace_back(c,false);}}return root;
}
inline auto toGraph(const CompressedTriePtr& root){std::vector<std::pair<Int,StaticString<char>>> p{{-1,{}}};auto dfs=[&](auto&& self,const CompressedTriePtr& node,Int id)->void{for(const auto& child:node->child)if(child){Int l=lcp(node->s,child->s),next=p.size();p.emplace_back(id,child->s.substr(l,child->s.len()));self(self,child,next);}};dfs(dfs,root,0);auto g=initWeightedDirectedGraph<StaticString<char>>(p.size());for(Int i=1;i<Int(p.size());++i)g.add_edge(p[i].first,i,p[i].second);return g;}
struct VirtualTrieNode {
    CompressedTriePtr current_node;StaticString<char> now;
    bool has_child(char c)const{assert('a'<=c&&c<='z');if(!now.len())return bool(current_node->child[c-'a']);return current_node->child[now[0]-'a']->s[current_node->s.len()+now.len()]==c;}
    VirtualTrieNode get_child(char c)const{assert(has_child(c));auto child=current_node->child[(now.len()?now[0]:c)-'a'];Int n=current_node->s.len(),length=now.len()+1;if(child->s.len()==n+length)return {child,now.substr(0,0)};return {current_node,child->s.substr(n,n+length)};}
    Int get_cnt()const{return now.len()?0:current_node->cnt;}Int get_subtree_sum()const{return now.len()?current_node->child[now[0]-'a']->subtree_sum:current_node->subtree_sum;}
};
inline VirtualTrieNode get_virtualnode(const CompressedTriePtr& n){return {n,n->s.substr(0,0)};}
inline auto get_child(const VirtualTrieNode& n,char c){return n.get_child(c);}inline bool has_child(const VirtualTrieNode& n,char c){return n.has_child(c);}inline Int get_cnt(const VirtualTrieNode& n){return n.get_cnt();}inline Int get_subtree_sum(const VirtualTrieNode& n){return n.get_subtree_sum();}
}
