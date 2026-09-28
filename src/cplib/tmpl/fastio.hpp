#pragma once
#include <cplib/tmpl/private/template_support.hpp>
#include <cplib/common.hpp>
#include <cplib/tmpl/private/fastio_native.hpp>
#include <span>
#include <sstream>
#include <array>
namespace cplib {
struct PrintOptions {FILE* f=stdout;std::string sepc=" ",endc="\n";bool flush=false;};
// 元の$相当の通常変換。浮動小数点や独自型はストリーム出力を利用する。
template<class T> std::string stringify(const T& value){return detail::template_support::stringValue(value);}
// 型と符号をコンパイル時に選ぶ。uint32単独/配列は元のSIMD専用経路を使用する。
template<class T=Int> __attribute__((target("avx2"))) T input(){if constexpr(std::is_same_v<T,std::string>){std::size_t length;auto p=detail::fastio_native::cplib_fio_read_token(&length);return std::string(p,length);}else{static_assert(std::integral<T>&&!std::is_same_v<T,bool>&&!std::is_same_v<T,char>);if constexpr(std::is_signed_v<T>)return T(detail::fastio_native::cplib_fio_read_int());else if constexpr(std::is_same_v<T,std::uint32_t>)return detail::fastio_native::cplib_fio_read_u32();else return T(detail::fastio_native::cplib_fio_read_uint());}}
template<class T=Int> __attribute__((target("avx2"))) std::vector<T> input(Int n){assert(n>=0);std::vector<T> out(n);if(n){if constexpr(std::is_same_v<T,std::string>){for(auto& x:out)x=input<T>();}else{static_assert(std::integral<T>&&!std::is_same_v<T,bool>&&!std::is_same_v<T,char>);if constexpr(std::is_signed_v<T>){if constexpr(sizeof(T)==8)detail::fastio_native::cplib_fio_read_i64_array(out.data(),n);else if constexpr(sizeof(T)==4)detail::fastio_native::cplib_fio_read_i32_array(out.data(),n);else if constexpr(sizeof(T)==2)detail::fastio_native::cplib_fio_read_i16_array(out.data(),n);else detail::fastio_native::cplib_fio_read_i8_array(out.data(),n);}else{if constexpr(sizeof(T)==8)detail::fastio_native::cplib_fio_read_u64_array(out.data(),n);else if constexpr(sizeof(T)==4)detail::fastio_native::cplib_fio_read_u32_array(reinterpret_cast<std::uint32_t*>(out.data()),n);else if constexpr(sizeof(T)==2)detail::fastio_native::cplib_fio_read_u16_array(out.data(),n);else detail::fastio_native::cplib_fio_read_u8_array(out.data(),n);}}}return out;}
inline Int ii(){return input<Int>();}inline std::vector<Int> lii(Int n){return input<Int>(n);}inline std::string si(){return input<std::string>();}
inline int fastioGetChar(){return detail::fastio_native::cplib_fio_get_char();}
inline void print_internal(const PrintOptions& prop,std::span<const std::string> args){for(std::size_t i=0;i<args.size();++i){std::fwrite(args[i].data(),1,args[i].size(),prop.f);const auto& sep=i+1==args.size()?prop.endc:prop.sepc;std::fwrite(sep.data(),1,sep.size(),prop.f);}if(prop.flush)std::fflush(prop.f);}
template<class... T> void print(const PrintOptions& prop,const T&... args){std::array<std::string,sizeof...(T)> values{stringify(args)...};print_internal(prop,values);}
template<class T> std::string join(std::span<const T> values,std::string_view sep=""){
 if(values.empty())return {};
 if constexpr(std::integral<T>&&(sizeof(T)==4||sizeof(T)==8)){constexpr std::size_t digits=sizeof(T)==8?20:std::is_signed_v<T>?11:10;std::string out(values.size()*digits+(values.size()-1)*sep.size(),'\0');std::size_t written;
  if constexpr(std::is_signed_v<T>){if constexpr(sizeof(T)==8)written=detail::fastio_native::cplib_fio_join_i64(values.data(),values.size(),out.data(),sep.data(),sep.size());else written=detail::fastio_native::cplib_fio_join_i32(values.data(),values.size(),out.data(),sep.data(),sep.size());}else{if constexpr(sizeof(T)==8)written=detail::fastio_native::cplib_fio_join_u64(values.data(),values.size(),out.data(),sep.data(),sep.size());else written=detail::fastio_native::cplib_fio_join_u32(values.data(),values.size(),out.data(),sep.data(),sep.size());}out.resize(written);return out;
 }else if constexpr(requires(T x){T::umod();x.val();}){std::vector<std::uint32_t> canonical(values.size());for(std::size_t i=0;i<values.size();++i)canonical[i]=values[i].val();return join<std::uint32_t>(canonical,sep);}else{std::string out;for(std::size_t i=0;i<values.size();++i){if(i)out+=sep;out+=stringify(values[i]);}return out;}
}
template<class T,class A> std::string join(const std::vector<T,A>& values,std::string_view sep=""){if constexpr(std::is_same_v<T,bool>){std::string out;for(std::size_t i=0;i<values.size();++i){if(i)out+=sep;out+=values[i]?"true":"false";}return out;}else return join<T>(std::span<const T>(values.data(),values.size()),sep);}
template<class T,std::size_t N> std::string join(const std::array<T,N>& values,std::string_view sep=""){return join<T>(std::span<const T>(values),sep);}
// 配列出力は固定長バッファを使用する。巨大な区切り文字も分割して書き出す。
template<class T> void print_array(std::span<const T> values,std::string_view sep=" ",FILE* file=stdout){if(values.empty()){std::fwrite("\n",1,1,file);return;}if constexpr(std::integral<T>&&(sizeof(T)==4||sizeof(T)==8)){if constexpr(std::is_signed_v<T>){if constexpr(sizeof(T)==8)detail::fastio_native::cplib_fio_print_array_i64(file,values.data(),values.size(),sep.data(),sep.size());else detail::fastio_native::cplib_fio_print_array_i32(file,values.data(),values.size(),sep.data(),sep.size());}else{if constexpr(sizeof(T)==8)detail::fastio_native::cplib_fio_print_array_u64(file,values.data(),values.size(),sep.data(),sep.size());else detail::fastio_native::cplib_fio_print_array_u32(file,values.data(),values.size(),sep.data(),sep.size());}}else if constexpr(requires(T x){T::umod();x.val();}){std::vector<std::uint32_t> canonical(values.size());for(std::size_t i=0;i<values.size();++i)canonical[i]=values[i].val();print_array<std::uint32_t>(canonical,sep,file);}else{auto text=join<T>(values,sep);std::fwrite(text.data(),1,text.size(),file);std::fwrite("\n",1,1,file);}}
template<class T,class A> void print_array(const std::vector<T,A>& values,std::string_view sep=" ",FILE* file=stdout){if constexpr(std::is_same_v<T,bool>){auto text=join(values,sep);std::fwrite(text.data(),1,text.size(),file);std::fwrite("\n",1,1,file);}else print_array<T>(std::span<const T>(values.data(),values.size()),sep,file);}
template<class T,std::size_t N> void print_array(const std::array<T,N>& values,std::string_view sep=" ",FILE* file=stdout){print_array<T>(std::span<const T>(values),sep,file);}
// 1値の整数出力はstdio自身のバッファへ直接書き込む経路を維持する。
template<class T> void print_one(const T& value,FILE* file=stdout){if constexpr(std::integral<T>&&(sizeof(T)==4||sizeof(T)==8)){if constexpr(std::is_signed_v<T>){if constexpr(sizeof(T)==8)detail::fastio_native::cplib_fio_print_one_i64(file,value);else detail::fastio_native::cplib_fio_print_one_i32(file,value);}else{if constexpr(sizeof(T)==8)detail::fastio_native::cplib_fio_print_one_u64(file,value);else detail::fastio_native::cplib_fio_print_one_u32(file,value);}}else if constexpr(requires{T::umod();value.val();})detail::fastio_native::cplib_fio_print_one_u32(file,value.val());else print(PrintOptions{file},value);}
struct Separator {std::string value;};inline Separator sep(std::string value){return {std::move(value)};}inline Separator sep(char c){return {std::string(1,c)};}
template<class R> struct Splat {const R& values;};template<class R> Splat<R> splat(const R& values){return {values};}
template<class R> void print(Splat<R> values){print_array(values.values);}template<class R> void print(Separator separator,Splat<R> values){print_array(values.values,separator.value);}
template<class T> requires(!std::is_same_v<T,Separator>&&!std::is_same_v<T,PrintOptions>) void print(const T& value){print_one(value);}
template<class... T> void print(Separator separator,const T&... values){if constexpr(sizeof...(T)>0&&(std::integral<T>&&...)){using V=std::common_type_t<T...>;if constexpr((std::is_same_v<V,T>&&...)){std::array<V,sizeof...(T)> a{values...};print_array(a,separator.value);}else print(PrintOptions{stdout,separator.value,"\n",false},values...);}else print(PrintOptions{stdout,separator.value,"\n",false},values...);}

inline void print(){} // 元の0引数版は改行も出力しない。
template<class T,class U,class... V> requires(!std::is_same_v<T,Separator>&&!std::is_same_v<T,PrintOptions>) void print(const T& a,const U& b,const V&... rest){print(sep(" "),a,b,rest...);}
struct InputSource {template<class T> T input(){return cplib::input<T>();}};
}
