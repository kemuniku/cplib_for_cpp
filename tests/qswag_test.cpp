#include <cplib/collections/QSWAG.hpp>
#include <random>
#include <deque>
using namespace cplib;
int main(){std::mt19937 rng(781);auto op=[](const std::string& a,const std::string& b){return a+b;};auto queue=initSWAG(op,std::string());std::deque<std::string> qq;for(int trial=0;trial<5000;++trial){std::string x(1,char('a'+rng()%26)),want;        if(qq.empty()||rng()%2){queue.push(x);qq.push_back(x);}else{assert(queue.pop()==qq.front());qq.pop_front();}want.clear();for(auto& x:qq)want+=x;assert(queue.fold()==want);
        if(!qq.empty())assert(queue[from_end(1)]==qq.back());
}
    auto sum=[](Int x,Int y){return x+y;};std::vector<Int> input{1,2,3,4,5,6};auto ok=[](Int x){return x<=7;};std::vector<Int> want;for(Int l=0;l<Int(input.size());++l){Int r=l,s=0;while(r<Int(input.size())&&s+input[r]<=7)s+=input[r++];want.push_back(r);}assert(get_maxrights(input,sum,Int(0),ok)==want);assert(cplib::get_maxrights(input,sum,Int(0),ok)==want);
}
