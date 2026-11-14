#pragma once
#include <cplib/collections/private/bitset_avx512_fuse_block.hpp>

namespace cplib {
// Nimの構文マクロをC++のラムダとSSA式で表す。
// fuse([&](auto& f){auto a=f.var(x),b=f.var(y);auto t=a&b;a=t+b;a<<=1;a.low(true);});
// bodyで構築した式をまとめて評価する。式節点と登録操作の合計をm、ビット数をnとしてO(m³+m(1+n/64))とbodyの計算量。
template <class Function> void fuse(Function &&body) {
    FuseBlock block;
    std::invoke(std::forward<Function>(body), block);
    block.run();
}

// 単一代入では加算がなければAVX2も利用する。
// fuse(dst, [](auto a,auto b){return (a<<3)^~b;}, x,y);
// 入力から作った式を評価しdstに代入する。式節点と登録操作の合計をm、ビット数をnとしてO(m³+m(1+n/64))とexpressionの計算量。
template <class Function, class... Inputs>
void fuse(SimdBitSet<512> &dst, Function &&expression, const Inputs &...inputs) {
    FuseBlock block;
    block.assign(dst, std::invoke(std::forward<Function>(expression), block.readValues(inputs)...));
    block.run(true);
}
}
