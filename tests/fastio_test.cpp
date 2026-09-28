#include <cplib/tmpl/fastio.hpp>
#include <cplib/tmpl/replayable_input.hpp>
#include <cplib/modint/modint.hpp>
#include <unistd.h>
#include <random>
using namespace cplib;
namespace fio=cplib;
Int fastio_read_other();
bool owned_mapping=false;
void clear_input(){auto* s=cplib::detail::fastio_native::cplib_fio_input_state();if(owned_mapping&&s->mapped)munmap(const_cast<char*>(s->mapped),s->length);owned_mapping=false;*s={};}
void mapped(const std::string& text){clear_input();auto* s=cplib::detail::fastio_native::cplib_fio_input_state();s->mapped=text.data();s->length=text.size();s->initialized=true;}
void file_input(const std::string& text){clear_input();FILE* f=std::tmpfile();assert(f);assert(std::fwrite(text.data(),1,text.size(),f)==text.size());std::fflush(f);std::rewind(f);assert(dup2(fileno(f),fileno(stdin))>=0);std::fclose(f);std::clearerr(stdin);assert(std::fseek(stdin,0,SEEK_SET)==0);owned_mapping=true;}
std::string contents(FILE* f){std::fflush(f);long size=std::ftell(f);assert(size>=0);std::rewind(f);std::string s(size,'\0');assert(std::fread(s.data(),1,s.size(),f)==s.size());return s;}
template<class F> std::string capture(F fn){std::fflush(stdout);int saved=dup(fileno(stdout));FILE* f=std::tmpfile();assert(f&&saved>=0);assert(dup2(fileno(f),fileno(stdout))>=0);fn();std::fflush(stdout);assert(dup2(saved,fileno(stdout))>=0);close(saved);auto out=contents(f);std::fclose(f);return out;}
template<class T> void integers(){std::mt19937_64 rng(5381);std::vector<T> values{std::numeric_limits<T>::min(),std::numeric_limits<T>::max(),0,1};for(int i=0;i<1300;++i)values.push_back(T(rng()));std::string text,expected;for(std::size_t i=0;i<values.size();++i){auto s=std::to_string(values[i]);text+=s+(i%3==0?" \n\t":" ");if(i)expected+="|";expected+=s;}for(bool stream:{false,true}){if(stream)file_input(text);else mapped(text);auto actual=fio::input<T>(values.size());assert(actual==values);assert(fio::input<T>()==0);}assert(fio::join(values,"|")==expected);FILE* f=std::tmpfile();assert(f);fio::print_array(values,"|",f);assert(contents(f)==expected+"\n");std::fclose(f);f=std::tmpfile();for(T v:values)fio::print_one(v,f);std::string one;for(T v:values)one+=std::to_string(v)+'\n';assert(contents(f)==one);std::fclose(f);}
struct Source {Int* position;template<class T> T input(){if constexpr(std::is_same_v<T,std::string>)return "s"+std::to_string((*position)++);else return T((*position)++);}};
struct Large {Int a,b,c;bool operator==(const Large&)const=default;};
void replay(){Int position=0;replayableInput(Source{&position},[&](auto& input){assert(input.ii()==0);input.peekInput([&]{assert(input.ii()==1);assert(input.si()=="s2");input.peekInput([&]{assert(input.ii()==3);assert(input.si()=="s4");});assert(input.ii()==3);assert(input.si()=="s4");});assert(position==5);assert(input.ii()==1);assert(input.si()=="s2");assert(input.ii()==3);assert(input.si()=="s4");assert(input.ii()==5);try{input.peekInput([&]{assert(input.si()=="s6");throw std::runtime_error("test");});assert(false);}catch(const std::runtime_error&){}assert(input.si()=="s6");input.peekInput([&]{assert(input.ii()==7);});bool mismatch=false;try{input.si();}catch(const std::invalid_argument&){mismatch=true;}assert(mismatch);assert(input.ii()==7);});assert(position==8);
 ReplayInputState<> state;Int reads=0;auto next=[&]{++reads;return Large{1,2,3};};state.peekInput([&]{assert(state.read<Large>(next)==(Large{1,2,3}));state.peekInput([&]{assert(state.read<std::string>([]{return "abc";})=="abc");});});assert(state.read<Large>(next)==(Large{1,2,3})&&reads==1);assert(state.read<std::string>([]{return "wrong";})=="abc");
 position=0;replayableInput<true>(Source{&position},[&](auto& input){for(int trial=0;trial<1000;++trial){Int x=position;input.peekInput([&]{assert(input.lii(3)==std::vector<Int>({x,x+1,x+2}));});assert(input.lii(3)==std::vector<Int>({x,x+1,x+2}));}});assert(position==3000);
}
int main(){integers<std::int8_t>();integers<std::int16_t>();integers<std::int32_t>();integers<std::int64_t>();integers<std::uint8_t>();integers<std::uint16_t>();integers<std::uint32_t>();integers<std::uint64_t>();
 std::string input="11 22";mapped(input);assert(fio::ii()==11&&fastio_read_other()==22);input.clear();std::vector<std::uint32_t> expected;for(int i=0;i<500;++i){expected.push_back(123456789+i);input+=std::to_string(expected.back())+' ';}for(bool stream:{false,true}){if(stream)file_input(input);else mapped(input);assert(fio::input<std::uint32_t>(500)==expected);}
 input=std::string((1<<20)-8,' ')+"-9223372036854775808 "+std::string((1<<20)+3,'0')+"4294967295 "+std::string((1<<21)+9,'x')+" tail";for(bool stream:{false,true}){if(stream)file_input(input);else mapped(input);assert(fio::ii()==std::numeric_limits<Int>::min());assert(fio::input<std::uint32_t>()==4294967295u);assert(fio::si()==std::string((1<<21)+9,'x'));assert(fio::si()=="tail");assert(fio::si().empty());}
 input="  +12 -1 +18446744073709551615 +00000042 z";file_input(input);assert(fio::ii()==12);assert(fio::input<UInt>()==UInt(-1));assert(fio::input<UInt>()==UInt(-1));assert(fio::input<std::uint32_t>()==42);assert(fio::si()=="z");
#ifdef CPLIB_FASTIO_MMAP
 assert(cplib::detail::fastio_native::cplib_fio_input_state()->mapped);
#else
 assert(!cplib::detail::fastio_native::cplib_fio_input_state()->mapped);
#endif
 assert(fio::join(std::array<int,3>{1,2,3},",")=="1,2,3");assert(fio::join(std::vector<bool>{true,false},",")=="true,false");using Mint=StaticMontgomeryModint<998244353>;std::vector<Mint> mods{Mint(-1),2,3};assert(fio::join(mods,",")=="998244352,2,3");std::string separator(70000,'s');FILE* f=std::tmpfile();fio::print_array(std::array<Int,3>{1,2,3},separator,f);assert(contents(f)=="1"+separator+"2"+separator+"3\n");std::fclose(f);
 auto output=capture([&]{fio::print(1);fio::print(2,3);fio::print(fio::sep(','),4,5,6);fio::print(fio::sep("::"),7,"last string");fio::print(fio::splat(mods));fio::print(fio::sep(','),fio::splat(std::vector<Int>{9,10}));fio::print(fio::PrintOptions{stdout,"-","!",true},"a","b");fio::print();fio::print(fio::sep(","));std::fwrite("end",1,3,stdout);});assert(output=="1\n2 3\n4,5,6\n7::last string\n998244352 2 3\n9,10\na-b!end");replay();clear_input();}
