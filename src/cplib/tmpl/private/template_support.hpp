#pragma once
#include <cplib/common.hpp>
#include <cplib/utils/constants.hpp>
#include <cplib/utils/backwards_index.hpp>
#include <array>
#include <string_view>
#include <sstream>
#include <iomanip>
#include <tuple>
#include <queue>
#include <set>
#include <unordered_set>
#include <unordered_map>
#include <map>
#include <deque>
#include <charconv>
#include <cstdio>
#include <cmath>
#include <functional>
namespace cplib::detail::template_support {
struct PrintOptions {FILE* f=stdout;std::string sepc=" ",endc="\n";bool flush=false;};
inline std::string quoted(std::string_view text){std::string out="\"";for(unsigned char c:text){if(c=='\n')out+="\\n";else if(c=='\r')out+="\\r";else if(c=='\t')out+="\\t";else if(c=='\\'||c=='\"'){out+='\\';out+=char(c);}else out+=char(c);}return out+'"';}
template<class T> inline constexpr bool IsArray=false;
template<class T,std::size_t N> inline constexpr bool IsArray<std::array<T,N>> = true;
template<class T> std::string stringValue(const T&);
template<class T> std::string repr(const T& value){if constexpr(std::is_convertible_v<T,std::string_view>)return quoted(std::string_view(value));else if constexpr(std::is_same_v<T,char>)return std::string("'")+value+"'";else return stringValue(value);}
template<class T> std::string stringValue(const T& value){if constexpr(std::is_same_v<T,bool>)return value?"true":"false";else if constexpr(std::is_convertible_v<T,std::string_view>)return std::string(std::string_view(value));else if constexpr(std::is_same_v<T,char>)return std::string(1,value);else if constexpr(std::integral<T>)return std::to_string(value);else if constexpr(std::floating_point<T>){char buffer[128];auto r=std::to_chars(buffer,buffer+sizeof(buffer),value);std::string out(buffer,r.ptr);if(out.find_first_of(".eE") ==std::string::npos&&std::isfinite(value))out+=".0";return out;}else if constexpr(requires{value.str();})return value.str();else if constexpr(requires{value.val();})return std::to_string(value.val());else if constexpr(IsArray<T>){std::string out="[";bool first=true;for(const auto& x:value){if(!first)out+=", ";first=false;out+=repr(x);}return out+"]";}else if constexpr(requires{std::tuple_size<T>::value;}){std::string out="(";std::size_t i=0;std::apply([&](const auto&... x){((out+=(i++?", ":"")+repr(x)),...);},value);if constexpr(std::tuple_size_v<T> == 1)out+=',';return out+")";}else if constexpr(requires{value.begin();value.end();}){std::string out="@[";bool first=true;for(const auto& x:value){if(!first)out+=", ";first=false;out+=repr(x);}return out+"]";}else{std::ostringstream os;os<<value;return os.str();}}
inline void write(const PrintOptions& prop,std::span<const std::string> values){for(std::size_t i=0;i<values.size();++i){std::fwrite(values[i].data(),1,values[i].size(),prop.f);const auto& sep=i+1==values.size()?prop.endc:prop.sepc;std::fwrite(sep.data(),1,sep.size(),prop.f);}if(prop.flush)std::fflush(prop.f);}
template<class... T> void print(const PrintOptions& prop,const T&... values){std::array<std::string,sizeof...(T)> text{stringValue(values)...};write(prop,text);}
template<class... T> void print(const T&... values){print(PrintOptions{},values...);}
// 組み込み整数の演算子は再定義できないため、床除算・余りは名前付き関数で公開する。
inline Int floor_mod(Int x,Int y){Int result=x%y;if(y>0&&result<0)result+=y;if(y<0&&result>0)result+=y;return result;}
inline Int floor_div(Int x,Int y){Int result=x/y,remainder=x%y;if(remainder&&((x<0)!=(y<0)))--result;return result;}
inline void floor_mod_assign(Int& x,Int y){x=floor_mod(x,y);}inline void floor_div_assign(Int& x,Int y){x=floor_div(x,y);}
template<class T> bool chmin(T& x,const T& y){if(x>y){x=y;return true;}return false;}template<class T> bool chmax(T& x,const T& y){if(x<y){x=y;return true;}return false;}template<class T,class U> void min_assign(T& x,U y){if(x>y)x=y;}template<class T,class U> void max_assign(T& x,U y){if(x<y)x=y;}
inline bool bit(Int x,Int n){assert(n>=0&&n<64);return (UInt(x)>>n)&1;}inline void set_bit(Int& x,Int n,bool value){assert(n>=0&&n<64);UInt v=UInt(x),mask=UInt(1)<<n;v=value?v|mask:v&~mask;x=std::bit_cast<Int>(v);}inline Int at(char x,char a='0'){return Int(static_cast<unsigned char>(x))-Int(static_cast<unsigned char>(a));}
inline Int power(Int x,Int n){assert(n>=0);Int result=1;while(n){if(n&1)result*=x;n>>=1;if(n)x*=x;}return result;}
struct Range {Int first,last,step;struct Iterator {Int value,last,step;Int operator*()const{return value;}Iterator& operator++(){value+=step;return *this;}bool operator!=(std::default_sentinel_t)const{return step>0?value<last:step<0?value>last:false;}};Iterator begin()const{return {first,last,step};}std::default_sentinel_t end()const{return {};}};
inline Range range(Int first,Int last,Int step){return {first,last,step};}inline Range range(Int first,Int last){return {first,last,1};}inline Range range(Int last){return {0,last,1};}
template<class T,class F> auto mapIt(ClosedSlice<T,T> slice,F operation){using Out=std::remove_cvref_t<std::invoke_result_t<F,T>>;std::vector<Out> out;for(T value=slice.a;value<=slice.b;++value)out.push_back(operation(value));return out;}
template<class T,class F> auto filterIt(ClosedSlice<T,T> slice,F predicate){std::vector<T> out;for(T value=slice.a;value<=slice.b;++value)if(predicate(value))out.push_back(value);return out;}
inline Int sum(ClosedSlice<Int,Int> slice){return floor_div((slice.a+slice.b)*std::max(Int(0),slice.b-slice.a+1),2);}
template<class T> void dump(const std::vector<std::vector<T>>& values){for(const auto& row:values)print(row);}
// debugマクロの引数名を、文字列リテラルと括弧の入れ子を避けて分割する。
inline std::vector<std::string_view> debugNames(std::string_view names){std::vector<std::string_view> out;std::size_t start=0;int depth=0;char quote=0;bool escape=false;for(std::size_t i=0;i<names.size();++i){char c=names[i];if(quote){if(escape)escape=false;else if(c=='\\')escape=true;else if(c==quote)quote=0;}else if(c=='\"'||c=='\'')quote=c;else if(c=='('||c=='['||c=='{')++depth;else if(c==')'||c==']'||c=='}')--depth;else if(c==','&&depth==0){out.push_back(names.substr(start,i-start));start=i+1;}}if(start<names.size())out.push_back(names.substr(start));for(auto& s:out){while(!s.empty()&&s.front()==' ')s.remove_prefix(1);while(!s.empty()&&s.back()==' ')s.remove_suffix(1);}return out;}
template<class... T> void debugPrintValues(std::string_view names,const T&... values){auto n=debugNames(names);std::size_t i=0;auto one=[&](const auto& v){if(i)std::fwrite(", ",1,2,stderr);auto name=n.at(i++);std::string text;if(!name.empty()&&name.front()=='"')text=detail::template_support::quoted(stringValue(v));else text=std::string(name)+": "+stringValue(v);std::fwrite(text.data(),1,text.size(),stderr);};(one(values),...);if(sizeof...(T))std::fwrite("\n",1,1,stderr);std::fflush(stderr);}
}
#if defined(CPLIB_DEBUG) || defined(debug)
#define CPLIB_DEBUG_PRINT(...) ::cplib::detail::template_support::debugPrintValues(#__VA_ARGS__, __VA_ARGS__)
#define CPLIB_DBLOCK(...) do { __VA_ARGS__ } while(false)
#else
#define CPLIB_DEBUG_PRINT(...) ((void)0)
#define CPLIB_DBLOCK(...) ((void)0)
#endif
#define CPLIB_SYMBOL_NAME(...) #__VA_ARGS__
