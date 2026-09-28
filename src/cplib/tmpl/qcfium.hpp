#pragma once
// 元と同じGCCの翻訳単位設定。Clangではビルド時に-mavx2 -O3 -funroll-loopsを指定する。
#if defined(__GNUC__) && !defined(__clang__)
#pragma GCC target("avx2")
#pragma GCC optimize("O3")
#pragma GCC optimize("unroll-loops")
#endif
