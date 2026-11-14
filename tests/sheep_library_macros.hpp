#pragma once
#include <cplib/tmpl/sheep.hpp>

// 元のsheep.hppの短いマクロでもライブラリを使えることを検証する。
// 手元のテンプレートで関数化されていても、テストでは元の呼び出しを再現する。
#ifndef len
#define len(value) ((::cplib::ll)(value).size())
#endif
#ifndef out
#define out(...)                                                                                  \
    overload_inner_output(__VA_ARGS__, inner_output8, inner_output7, inner_output6, inner_output5, \
                          inner_output4, inner_output3, inner_output2, inner_output1)(__VA_ARGS__)
#endif
