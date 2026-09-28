#include <cplib/str/lcs_bitset.hpp>
#include <cplib/str/edit_distance_bitset.hpp>
#include <random>
using namespace cplib;
template<class A,class B> Int lcs_naive(const A& a,const B& b){std::vector<Int> dp(b.size()+1);for(auto x:a){Int prev=0;for(std::size_t j=0;j<b.size();++j){Int save=dp[j+1];dp[j+1]=x==b[j]?prev+1:std::max(dp[j],dp[j+1]);prev=save;}}return dp.back();}
template<class A,class B> bool subsequence(const A& a,const B& b){std::size_t j=0;for(auto x:b)if(j<a.size()&&a[j]==x)++j;return j==a.size();}
template<class A> void check_lcs(const A& a,const A& b){auto n=lcs_naive(a,b);assert(cplib::LCS(a,b)==n);auto answer=cplib::restoreLCS(a,b);assert(Int(answer.size())==n&&subsequence(answer,a)&&subsequence(answer,b));}
Int edit_naive(std::string_view a,std::string_view b){std::vector<Int> dp(b.size()+1);std::iota(dp.begin(),dp.end(),0);for(char x:a){Int prev=dp[0]++;for(std::size_t j=0;j<b.size();++j){Int save=dp[j+1];dp[j+1]=std::min({prev+(x!=b[j]),dp[j]+1,dp[j+1]+1});prev=save;}}return dp.back();}
struct EqualityOnly {int value;bool operator==(const EqualityOnly&)const=default;};
int main(){
 std::mt19937 rng(8546);
 for(int trial=0;trial<300;++trial){int n=rng()%150,m=rng()%150,k=1+rng()%8;std::string a(n,'\0'),b(m,'\0');for(char& x:a)x=char(rng()%k+(trial%2?248:0));for(char& x:b)x=char(rng()%k+(trial%2?248:0));check_lcs(a,b);assert(editDistance_bitset(a,b)==edit_naive(a,b));}
 for(int n:{0,1,63,64,65,511,512,513,1025}){std::string a(n,'a'),b(n,'b');assert(editDistance_bitset(a,b)==n&&editDistance_bitset(a,a)==0);check_lcs(a,a);}
 for(int trial=0;trial<50;++trial){std::vector<int> a(200),b(160);for(auto& x:a)x=rng()%1000;for(auto& x:b)x=rng()%1000;check_lcs(a,b);}
 check_lcs(std::vector<EqualityOnly>{{1},{2},{1},{3}},std::vector<EqualityOnly>{{2},{1},{4},{3}});
 check_lcs(std::vector<bool>{true,false,true,true,false},std::vector<bool>{false,true,false});
 check_lcs(std::vector<std::string>{"abc","d","abc","f"},std::vector<std::string>{"d","f"});
}
