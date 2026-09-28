#define CPLIB_DEBUG 1
#include "template_test_support.hpp"
#include <cplib/tmpl/combined.hpp>
#include <cplib/tmpl/optimize.hpp>
#include <cplib/modint/modint.hpp>
namespace ct=cplib;namespace ts=cplib::detail::template_support;
Int combined_read_other();
void old_input(const std::string& text,bool stream){file_input(text);cplib::detail::combined_input_native::cursor=cplib::detail::combined_input_native::length=0;cplib::detail::combined_input_native::mapped=nullptr;cplib::detail::combined_input_native::initialized=stream;}
void release(){if(cplib::detail::combined_input_native::mapped){munmap(const_cast<char*>(cplib::detail::combined_input_native::mapped),cplib::detail::combined_input_native::length);cplib::detail::combined_input_native::mapped=nullptr;}}
void inputs(){
 for(bool stream:{true,false}){old_input("  -9223372036854775808 9223372036854775807 12345678 123456789 4294967295 -1 abcd tail",stream);assert(ct::ii()==std::numeric_limits<Int>::min());assert(combined_read_other()==std::numeric_limits<Int>::max());assert(ct::lii(2)==std::vector<Int>({12345678,123456789}));assert(ct::lii2(2)==std::vector<std::uint32_t>({4294967295u,4294967295u}));assert(ct::si()=="abcd"&&ct::si()=="tail"&&ct::si().empty());assert(ct::ii()==0);assert(bool(cplib::detail::combined_input_native::mapped)==!stream);release();}
 // 64 MiB の境界を跨ぐ負整数と文字列。旧版特有のバッファ経路。
 old_input(std::string(cplib::detail::combined_input_native::buffer_size-7,' ')+"-9223372036854775808 hello",true);assert(ct::ii()==std::numeric_limits<Int>::min());assert(ct::si()=="hello");
}
template<class T> void joins(){std::mt19937_64 rng(912);std::vector<T> a{std::numeric_limits<T>::min(),std::numeric_limits<T>::max(),0};for(int i=0;i<1000;++i)a.push_back(T(rng()));std::string expected;for(std::size_t i=0;i<a.size();++i){if(i)expected+="::";expected+=std::to_string(a[i]);}assert(ct::join(a,"::")==expected);}
int main(){inputs();joins<int8_t>();joins<int16_t>();joins<int32_t>();joins<int64_t>();joins<uint8_t>();joins<uint16_t>();joins<uint32_t>();joins<uint64_t>();
 using Mint=cplib::StaticMontgomeryModint<998244353>;assert(ct::join(std::vector<Mint>{-1,2},",")=="998244352,2");assert(ct::join(std::vector<std::string>{"","x",""},",")==",x,");
 assert(ts::stringValue(std::vector<std::string>{"a","b"})=="@[\"a\", \"b\"]");assert(ts::stringValue(std::array<int,2>{1,2})=="[1, 2]");assert(ts::stringValue(std::tuple{1,"x"})=="(1, \"x\")");assert(ts::stringValue(1.0)=="1.0");
 auto out=capture([&]{ct::print(1,2);ct::print(ct::sep("::"),3,"end");ct::print(ct::splat(std::vector<Int>{4,5}));ct::print(ct::sep(","),ct::splat(std::vector<Int>{6,7}));ct::print(ts::PrintOptions{stdout,"-","!",true},8,9);ct::yes(false);});
 assert(out=="1 2\n3::end\n4 5\n6,7\n8-9!No\n");
 Int value=9;assert(capture([&]{CPLIB_DEBUG_PRINT(value,"text");},stderr)=="value: 9, \"text\"\n");CPLIB_DBLOCK(value=10;);assert(value==10);assert(std::string(CPLIB_SYMBOL_NAME(value))=="value");
 std::vector<Int> r;for(Int x:ts::range(5,-3,-2))r.push_back(x);assert(r==std::vector<Int>({5,3,1,-1}));for(Int x:ts::range(1,9,0)){(void)x;assert(false);}auto sq=ts::mapIt(cplib::closed_slice(Int(1),Int(4)),[](Int x){return x*x;});assert(sq==std::vector<Int>({1,4,9,16}));assert(ts::filterIt(cplib::closed_slice(Int(1),Int(4)),[](Int x){return x%2;})==std::vector<Int>({1,3}));assert(ts::sum(cplib::closed_slice(Int(-2),Int(4)))==7);assert(ts::floor_div(-7,3)==-3&&ts::floor_mod(-7,3)==2&&ts::floor_div(7,-3)==-3&&ts::floor_mod(7,-3)==-2);assert(ts::power(3,10)==59049);ts::set_bit(value,63,true);assert(ts::bit(value,63));assert(cplib::debug_build&&!cplib::optimized_build);
}
