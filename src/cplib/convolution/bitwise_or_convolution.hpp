#pragma once
#include <cplib/convolution/set_transform.hpp>

namespace cplib {
// c[k]=Σ(i|j=k) a[i]b[j]を返す。同じ2冪長、両方空なら空。時間 O(N log N)、領域 O(N)。
template <class T>
std::vector<T> bitwiseOrConvolution(std::vector<T> a, std::vector<T> b) {
    assert(a.size() == b.size());
    subsetZetaTransform(a);
    subsetZetaTransform(b);
    for (std::size_t i = 0; i < a.size(); ++i)
        a[i] *= b[i];
    subsetMobiusTransform(a);
    return a;
}
}
