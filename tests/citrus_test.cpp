#include "template_test_support.hpp"
#include <cplib/tmpl/citrus.hpp>
int main(){ std::istringstream source("12 1.25\nabc def 42\n");cplib::Input in(source);assert(in.input<Int>()==12);assert(in.input<double>()==1.25);assert(in.input<char>()=='a');assert(in.input<std::string>()=="b");assert(in.input<char>()=='c');assert(in.input<std::string>()=="def");assert(in.input<Int>()==42);assert(in.input<std::string>().empty());
 std::istringstream nested("1 2 3 4 5 6 7 x 8 y");cplib::reader=cplib::Input(nested);assert(cplib::input<Int>(2,3)==std::vector<std::vector<Int>>({{1,2,3},{4,5,6}}));assert((cplib::input<Int,std::string>(2)==std::vector<std::tuple<Int,std::string>>({{7,"x"},{8,"y"}})));
assert(cplib::pow(3,10,100)==49);assert(capture([]{cplib::print(1.25,std::vector<int>{1,2});cplib::print(cplib::PrintOptions{stdout,":","!",false},1,2);})=="1.2500000000000000 1 2\n1:2!");}
