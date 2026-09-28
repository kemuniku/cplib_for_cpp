#include <cplib/str/zalgorithm.hpp>
#include <cplib/str/manacher.hpp>
#include <cplib/str/lcp_naive.hpp>
#include <cplib/str/lcs.hpp>
#include <cplib/str/suffix_array.hpp>
#include <cplib/str/edit_distance.hpp>
#include <cplib/str/trie.hpp>
#include <cplib/str/aho_corasick.hpp>
#include <random>
#include <set>
using namespace cplib;
int main(){
    std::mt19937 rng(51);
    for(int trial=0;trial<800;++trial){
        Int n=rng()%70;std::string s;for(Int i=0;i<n;++i)s+=char('a'+rng()%4);
        auto z=zalgorithm(s);auto pal=get_palindromes(s);auto odd=manacher(s);
        for(Int i=0;i<n;++i){assert(z[i]==lcp_naive(s,s.substr(i)));Int radius=0;while(i-radius>=0&&i+radius<n&&s[i-radius]==s[i+radius])++radius;assert(odd[i]==radius);}
        for(Int c=0;c<2*n-1;++c){Int l=c/2,r=(c+1)/2;while(l>=0&&r<n&&s[l]==s[r]){--l;++r;}std::pair<Int,Int> want=l+1==r?std::pair<Int,Int>{-1,-1}:std::pair{l+1,r};assert(pal[c]==want);}
        std::vector<Int> expected(n);std::iota(expected.begin(),expected.end(),0);std::sort(expected.begin(),expected.end(),[&](Int a,Int b){return s.substr(a)<s.substr(b);});auto sa=suffix_array(s);assert(sa==expected);auto lcps=lcp_array(s,sa);for(Int i=0;i<n-1;++i)assert(lcps[i]==lcp_naive(s.substr(sa[i]),s.substr(sa[i+1])));
        std::vector<Int> v(s.begin(),s.end());assert(suffix_array(v)==sa);assert(suffix_array(v,255)==sa);
        std::string t;Int m=rng()%50;for(Int i=0;i<m;++i)t+=char('a'+rng()%4);
        std::vector<std::vector<Int>> dp(n+1,std::vector<Int>(m+1)),ed=dp;for(Int i=0;i<=n;++i)ed[i][0]=i;for(Int j=0;j<=m;++j)ed[0][j]=j;
        for(Int i=1;i<=n;++i)for(Int j=1;j<=m;++j){dp[i][j]=s[i-1]==t[j-1]?dp[i-1][j-1]+1:std::max(dp[i-1][j],dp[i][j-1]);ed[i][j]=std::min({ed[i-1][j]+1,ed[i][j-1]+1,ed[i-1][j-1]+(s[i-1]!=t[j-1])});}
        assert(LCS(s,t)==dp[n][m]);auto subseq=restoreLCS(s,t);assert(Int(subseq.size())==dp[n][m]);for(const auto& text:{s,t}){std::size_t pos=0;for(char c:text)if(pos<subseq.size()&&c==subseq[pos])++pos;assert(pos==subseq.size());}
        for(Int k:{Int(0),Int(1),Int(3),ed[n][m],std::max(n,m)})assert(editDistance(s,t,k)==(ed[n][m]<=k?ed[n][m]:-1));
    }
    std::string bytes;for(int i=255;i>=0;--i)bytes+=char(i);auto sa=suffix_array(bytes);for(Int i=0;i<256;++i)assert(sa[i]==255-i);
    std::vector<Int> large(10000);for(Int& x:large)x=Int(rng()%4000000000ULL)-2000000000LL;auto rmq=initRMQ(large);
    std::vector<std::int32_t> small(large.begin(),large.end());auto rmq32=initRMQ(small);std::vector<double> real(large.begin(),large.end());auto rmq_real=initRMQ(real);
    for(int trial=0;trial<10000;++trial){Int l=rng()%large.size(),r=l+1+rng()%(large.size()-l);Int want=*std::min_element(large.begin()+l,large.begin()+r);assert(rmq.query(l,r)==want);assert(rmq32.query(l,r)==want);assert(rmq_real.query(l,r)==want);}
    auto trie=initTrie<'a','c'>();std::multiset<std::string> words;
    for(int trial=0;trial<500;++trial){
        std::string word;Int length=rng()%7;for(Int i=0;i<length;++i)word+=char('a'+rng()%3);
        if(rng()%2){Int count=rng()%4;trie.incl(word,count);for(Int i=0;i<count;++i)words.insert(word);}
        else {trie.excl(word,2);for(int i=0;i<2;++i){auto it=words.find(word);if(it!=words.end())words.erase(it);}}
        assert(trie.len()==Int(words.size()));assert(trie.count(word)==Int(words.count(word)));assert(trie.lowerBound(word)==std::distance(words.begin(),words.lower_bound(word)));assert(trie.upperBound(word)==std::distance(words.begin(),words.upper_bound(word)));
        Int prefixes=0;for(auto& w:words)prefixes+=w.starts_with(word);assert(trie.countPrefix(word)==prefixes);
    }
    auto tp=initTriePointer(trie);tp&='a';auto old=tp;tp&='b';assert(tp.restoreString()=="ab" && old.restoreString()=="a");assert(tp.pop()=='b');assert(trie.toGraph().len==Int(trie.nodes.size()));
    std::vector<std::string> patterns{"","a","ab","bab","aba","b","ab"};auto ac=initAhoCorasick<'a','c'>(patterns);auto p=initAhoCorasickPointer(ac);std::string text;
    for(int trial=0;trial<500;++trial){
        if(trial){char c=char('a'+rng()%4);text+=c;p&=c;}
        Int count=0;std::set<Int> nodes;for(std::size_t i=0;i<patterns.size();++i)if(text.ends_with(patterns[i])){++count;nodes.insert(ac.patternNode(i));}
        assert(p.matchCount()==count);std::set<Int> actual;for(Int node:p.matches())actual.insert(node);assert(actual==nodes);
        auto prefix=p.restoreString();assert(text.ends_with(prefix));assert(ac.findNode(prefix)==p.nodeId());
    }
    assert(ac.toTrieGraph().edge_count()==ac.nodeCount()-1);assert(ac.toFailureGraph().edge_count()==ac.nodeCount()-1);
}
