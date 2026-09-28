#include <cplib/collections/bitset_avx512.hpp>
#include <random>
using namespace cplib;
using namespace cplib;
int main(){
 std::mt19937_64 rng(39429);
 for(Int n:{0,1,63,64,65,255,256,511,512,513,1024,4099})for(int trial=0;trial<80;++trial){
  BitSet a(n),b(n),c(n),d(n);for(Int j=0;j<n;++j){a[j]=int(rng()&1);b[j]=int(rng()&1);c[j]=int(rng()&1);d[j]=int(rng()&1);}
  BitSet out(n);
  fuse(out,[](auto x,auto y,auto z){return (x&y)|(~x&z);},a,b,c);assert(out==((a&b)|(~a&c)));
  fuse(out,[](auto x,auto y,auto z,auto w){return ~((x|y)^z)&w;},a,b,c,d);assert(out==(~((a|b)^c)&d));
  fuse(out,[](auto x,auto y,auto z){return ((x+y)^z)+(x&y);},a,b,c);assert(out==(((a+b)^c)+(a&b)));
  for(Int shift:{1,3,63,64,65,511,1025}){
   fuse(out,[&](auto x,auto y){return (x<<shift)^~(y>>shift);},a,b);assert(out==((a<<shift)^~(b>>shift)));
   auto dst=a;fuse(dst,[&](auto x,auto y){return (x<<shift)^~(y>>shift);},dst,b);assert(dst==out);
   dst=a;fuse(dst,[&](auto x){return x|(x<<shift);},dst);assert(dst==(a|(a<<shift)));
   dst=a;fuse(dst,[&](auto x){return x|(x>>shift);},dst);assert(dst==(a|(a>>shift)));
  }
  fuse(out,[](auto x,auto y){return ((x&y)<<3)<<5;},a,b);assert(out==(((a&b)<<3)<<5));
  auto x=a,y=b;bool flag=false;
  auto t=(a&b)+c;auto xx=(t^b)<<3;auto yy=(t|d)+xx;bool expected=yy.lastBit();if(n)xx[0]=true;
  fuse([&](auto& f){auto p=f.var(x),q=f.var(y);auto r=f.read(c),s=f.read(d);auto tmp=(p&q)+r;p=(tmp^q)<<3;q=(tmp|s)+p;f.lastBit(flag,q);if(n)p.low(true);});
  assert(x==xx&&y==yy&&flag==expected);
  // 同一変数を複数の名前で参照しても文の順序を保持する。
  x=a;fuse([&](auto& f){auto p=f.var(x),q=f.var(x);p^=f.read(b);q|=f.read(c);});assert(x==((a^b)|c));
  // 共有中間値・同一出力値・死んだ式・大きい内部シフトのフォールバック。
  x=a;y=b;fuse([&](auto& f){auto p=f.var(x),q=f.var(y);auto tmp=(p&q)>>65;auto unused=p+q;(void)unused;p=tmp;q=tmp;});assert(x==((a&b)>>65)&&y==x);
  x=BitSet(n+3);fuse(x,[](auto p,auto q){return p^q;},a,b);assert(x==(a^b));
  if(n){x=a;fuse([&](auto& f){auto p=f.var(x);p<<=1;p.low(false);f.lastBit(flag,p);});assert(x==(a<<1)&&flag==x.lastBit());}
 }
 for(Int n:{513,1024,10003}){BitSet a(n),one(n),out(n);a.fill();one[0]=true;fuse(out,[](auto a,auto b){return a+b;},a,one);assert(!out.any());a[n-1]=false;fuse([&](auto& f){auto p=f.var(a);auto tmp=p+f.read(one);p=tmp;});assert(a.popcount()==1&&a[n-1]);}
 bool threw=false;try{BitSet a(5),b(6),out(5);fuse(out,[](auto x,auto y){return x&y;},a,b);}catch(const std::invalid_argument&){threw=true;}assert(threw);
}
