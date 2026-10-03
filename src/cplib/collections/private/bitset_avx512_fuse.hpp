#pragma once
#include <cplib/collections/private/bitset_avx512_fuse_block.hpp>

namespace cplib {
// Nimの構文マクロをC++のラムダとSSA式で表す。
// fuse([&](auto& f){auto a=f.var(x),b=f.var(y);auto t=a&b;a=t+b;a<<=1;a.low(true);});
template <class Function> void fuse(Function &&body) {
    FuseBlock block;
    std::invoke(std::forward<Function>(body), block);
    block.run();
}

// 単一代入では加算がなければAVX2も利用する。
// fuse(dst, [](auto a,auto b){return (a<<3)^~b;}, x,y);
template <class Function, class... Inputs>
void fuse(SimdBitSet<512> &dst, Function &&expression, const Inputs &...inputs) {
    FuseBlock block;
    block.assign(dst, std::invoke(std::forward<Function>(expression), block.read(inputs)...));
    block.run(true);
}
}
