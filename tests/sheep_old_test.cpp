#include "template_test_support.hpp"
#include <cplib/tmpl/sheep_old.hpp>
template<class T> void joins(){std::mt19937_64 rng(912);std::vector<T> a{std::numeric_limits<T>::min(),std::numeric_limits<T>::max(),0};for(int i=0;i<1000;++i)a.push_back(T(rng()));std::string expected;for(std::size_t i=0;i<a.size();++i){if(i)expected+="::";expected+=std::to_string(a[i]);}assert(cplib::join(a,"::")==expected);}
int main(){ file_input(" 12\n 34 abc def\n");assert(cplib::lii(2)==std::vector<Int>({12,34}));assert(cplib::si()=="abc");assert(cplib::si()=="def");file_input(" word");assert(cplib::si().empty());assert(cplib::si()=="word");
joins<int8_t>();joins<int16_t>();joins<int32_t>();joins<int64_t>();joins<uint64_t>();assert(capture([]{cplib::takahashi();cplib::aoki();cplib::yes();})=="Takahashi\nAoki\nYes\n");}
