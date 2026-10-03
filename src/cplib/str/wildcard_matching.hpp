#pragma once
#include <cplib/convolution/convolution.hpp>
#include <string_view>

namespace cplib {
// 両側のワイルドカードに対応。3回の整数畳み込みで衝突なく判定。O((N+M)log(N+M))。
// sの各開始位置でtが一致するかを返す。空間O(s.size()+t.size())。
// 両文字列のwildは任意の1バイトに一致し、文字列はバイト単位で扱う。
// tが空ならs.size()+1個のtrue、t.size() > s.size()なら空列を返す。
// 0 < t.size() <= s.size()の場合、s.size()+t.size()-1 <= 2^24が必要。
inline std::vector<bool> wildcard_match(std::string_view s, std::string_view t, char wild = '?') {
    if (t.size() > s.size())
        return {};
    std::vector<bool> result(s.size() - t.size() + 1);
    if (t.empty()) {
        std::fill(result.begin(), result.end(), true);
        return result;
    }
    assert(s.size() <= (std::size_t(1) << 24) - t.size() + 1);
    std::array<std::vector<Int>, 3> a, b;
    for (Int k = 0; k < 3; ++k) {
        a[k].resize(s.size());
        b[k].resize(t.size());
    }
    for (std::size_t i = 0; i < s.size(); ++i)
        if (s[i] != wild) {
            Int x = static_cast<unsigned char>(s[i]) + 1;
            a[0][i] = x;
            a[1][i] = x * x;
            a[2][i] = x * x * x;
        }
    for (std::size_t i = 0; i < t.size(); ++i)
        if (t[i] != wild) {
            Int x = static_cast<unsigned char>(t[i]) + 1;
            auto j = t.size() - 1 - i;
            b[0][j] = x;
            b[1][j] = x * x;
            b[2][j] = x * x * x;
        }
    auto x = convolution_ll(a[0], b[2]), y = convolution_ll(a[1], b[1]),
         z = convolution_ll(a[2], b[0]);
    for (std::size_t i = 0; i < result.size(); ++i) {
        auto j = i + t.size() - 1;
        result[i] = x[j] + z[j] - 2 * y[j] == 0;
    }
    return result;
}
}
