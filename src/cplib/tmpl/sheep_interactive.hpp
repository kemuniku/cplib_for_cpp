#pragma once
// 対話用入出力を有効にする。このヘッダーを他のテンプレートより先にincludeする。
#define CPLIB_FASTIO_INTERACTIVE
#include <cplib/tmpl/fastio.hpp>
#include <cplib/tmpl/sheep.hpp>

namespace cplib::detail {
// sheepのcout出力も各出力後に確定する。
inline const bool interactiveOutputInitialized = [] {
    std::cout << std::unitbuf;
    return true;
}();
}
