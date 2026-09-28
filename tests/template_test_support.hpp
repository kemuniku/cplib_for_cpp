#pragma once
#include <cplib/tmpl/private/template_support.hpp>
#include <unistd.h>
#include <stdio_ext.h>
#include <random>
using cplib::Int;using cplib::UInt;
std::string contents(FILE* f){std::fflush(f);long size=std::ftell(f);assert(size>=0);std::rewind(f);std::string s(size,'\0');assert(std::fread(s.data(),1,s.size(),f)==s.size());return s;}
template<class F> std::string capture(F fn,FILE* target=stdout){std::fflush(target);int saved=dup(fileno(target));FILE* f=std::tmpfile();assert(f&&saved>=0);assert(dup2(fileno(f),fileno(target))>=0);fn();std::fflush(target);assert(dup2(saved,fileno(target))>=0);close(saved);auto out=contents(f);std::fclose(f);return out;}
void file_input(const std::string& text){__fpurge(stdin);FILE* f=std::tmpfile();assert(f);assert(std::fwrite(text.data(),1,text.size(),f)==text.size());std::fflush(f);std::rewind(f);std::string path="/proc/self/fd/"+std::to_string(fileno(f));assert(std::freopen(path.c_str(),"rb",stdin));std::fclose(f);}
