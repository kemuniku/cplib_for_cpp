#pragma once
#include <cplib/collections/bitset_avx512.hpp>
namespace cplib {
// Myers法。バイト単位、時間O(n ceil(m/64))、領域O(256 ceil(m/64))語。
// AVX512F環境では融合した主要ループを8語ずつ処理する。
inline Int editDistance_bitset(std::string_view s,std::string_view t){
 if(s.size()<t.size())return editDistance_bitset(t,s);
 Int m=t.size();if(!m)return s.size();
 std::array<BitSetAvx512,256> matches;
 for(Int i=0;i<m;++i){auto& mask=matches[static_cast<unsigned char>(t[i])];if(!mask.len())mask=BitSetAvx512(m);mask[i]=true;}
 BitSetAvx512 positive(m),negative(m);positive.fill();Int result=m;
 for(unsigned char c:s){auto& mask=matches[c];if(!mask.len())mask=BitSetAvx512(m);bool positiveLast=false,negativeLast=false;
  cplib::fuse([&](auto& f){
   auto p=f.var(positive),n=f.var(negative);auto equal=f.read(mask);
   auto vertical=equal|n;
   auto horizontal=(((equal&p)+p)^p)|equal;
   auto positiveHorizontal=~(horizontal|p)|n;
   auto negativeHorizontal=p&horizontal;
   auto shiftedPositive=f.low(positiveHorizontal<<1,true);
   p=~(vertical|shiftedPositive)|(negativeHorizontal<<1);
   n=shiftedPositive&vertical;
   f.lastBit(positiveLast,positiveHorizontal);f.lastBit(negativeLast,negativeHorizontal);
  });
  result+=Int(positiveLast)-Int(negativeLast);
 }
 return result;
}
}
