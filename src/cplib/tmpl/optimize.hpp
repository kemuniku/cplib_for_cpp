#pragma once

// C++はコンパイル時に外部コンパイラを呼べないため、tools/optimize.py が
// 元の optimize / optimizeCpp のビルド処理を担う。このヘッダは段階の判定用。
namespace cplib {
#if defined(CPLIB_SECOND_COMPILE)
inline constexpr bool second_compile = true;
#else
inline constexpr bool second_compile = false;
#endif
#if defined(CPLIB_DEBUG) || defined(debug)
inline constexpr bool debug_build = true;
#else
inline constexpr bool debug_build = false;
#endif
inline constexpr bool optimized_build = second_compile && !debug_build;
}
