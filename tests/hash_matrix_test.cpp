#include <cplib/matrix/rolling_hash_2d.hpp>
#include <random>
using namespace cplib;
int main(){std::mt19937 rng(7163);for(int trial=0;trial<200;++trial){Int h=1+rng()%25,w=1+rng()%25;std::vector<std::vector<Int>> a(h,std::vector<Int>(w));for(auto& row:a)for(auto& x:row)x=rng()%256;auto hm=initHashMatrix(a);for(int query=0;query<50;++query){Int top=rng()%(h+1),bottom=rng()%(h+1),left=rng()%(w+1),right=rng()%(w+1);if(top>bottom)std::swap(top,bottom);if(left>right)std::swap(left,right);auto sub=hm.get(closed_slice(top,bottom-1),closed_slice(left,right-1));std::vector<std::vector<Int>> copy;for(Int i=top;i<bottom;++i)copy.emplace_back(a[i].begin()+left,a[i].begin()+right);assert(sub==initHashMatrix(copy));for(Int i=0;i<bottom-top;++i)for(Int j=0;j<right-left;++j)assert(sub[i][j]==a[top+i][left+j]);if(top<bottom&&left<right){copy[0][0]+=1;assert(!(sub==initHashMatrix(copy)));auto row=sub[0];sub=hm;assert(row[0]==a[top][left]);}}}
 auto zeros=initHashMatrix(std::vector<std::vector<Int>>(3,std::vector<Int>(5)));assert(zeros==initHashMatrix(std::vector<std::vector<Int>>{}));
 using namespace cplib;for(Int n:{0,1,500000,500001,1000000000}){assert(calc_mod(mul(base_powi(n),inv_base_powi(n)))==1);assert(calc_mod(mul(base_powj(n),inv_base_powj(n)))==1);}
 std::vector<std::vector<char>> big(1,std::vector<char>(500005,'a'));big[0][500003]=char(255);auto hm=initHashMatrix(big);assert(hm.get(closed_slice(0,0),closed_slice(500002,500004))==initHashMatrix(std::vector<std::vector<char>>{{'a',char(255),'a'}}));
}
