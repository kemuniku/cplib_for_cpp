#include <cplib/str/rolling_hash.hpp>
#include <random>
using namespace cplib;
int main(){std::mt19937 rng(452);for(Int n:{0,1,2,5,64,500,501001}){std::string text(n,'a');for(char& c:text)c=char(rng()%256);auto h=initRollingHash(text);for(int q=0;q<200;++q){Int l=rng()%(n+1),r=rng()%(n+1);if(l>r)std::swap(l,r);auto sub=text.substr(l,r-l);if(!sub.empty()){auto ref=initRollingHash(sub);assert(h.query(closed_slice(l,r-1))==ref.query(closed_slice(Int(0),Int(sub.size())-1)));}}}std::vector<Int> numbers{5,7,1,5,7,1};auto h=initRollingHash(numbers);assert(h.query(closed_slice(0,2))==h.query(closed_slice(3,5)));}
