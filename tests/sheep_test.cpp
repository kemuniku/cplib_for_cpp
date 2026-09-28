#include "template_test_support.hpp"
#include <cplib/tmpl/sheep.hpp>
int main(){assert(capture([]{cplib::yes();cplib::print(std::vector<int>{1,2});})=="Yes\n@[1, 2]\n");auto a=cplib::vv<int>(2,3,7);assert(a[1][2]==7);assert(cplib::powll(2,10)==1024);}
