#pragma once
#include <cplib/convolution/ntt.hpp>
#include <cplib/convolution/private/convolution_common.hpp>

namespace cplib {
using mint754974721 = StaticMontgomeryModint<754974721>;
using mint167772161 = StaticMontgomeryModint<167772161>;
using mint469762049 = StaticMontgomeryModint<469762049>;

template <class T> auto convolution_naive(const std::vector<T> &f, const std::vector<T> &g) {
    return detail::convolution_schoolbook(f, g);
}

// 旧版のスカラーNTT経路。長さ60以下だけ愚直積へ切り替える。O(N log N)。
template <class T>
    requires(BarrettModint<T> || MontgomeryModint<T>)
std::vector<T> convolution(std::vector<T> f, std::vector<T> g) {
    if (f.empty() || g.empty())
        return {};
    if (std::min(f.size(), g.size()) <= 60)
        return detail::convolution_schoolbook(f, g);
    std::size_t deg = f.size() + g.size() - 1, l = std::bit_ceil(deg);
    f.resize(l);
    g.resize(l);
    ntt(f);
    ntt(g);
    for (std::size_t i = 0; i < l; ++i)
        f[i] *= g[i];
    intt(f);
    f.resize(deg);
    return f;
}

inline auto convolution_ll(const std::vector<Int> &f, const std::vector<Int> &g) {
    auto run = [&]<class M>() {
        std::vector<M> a(f.begin(), f.end()), b(g.begin(), g.end());
        auto c = cplib::convolution(std::move(a), std::move(b));
        std::vector<std::uint32_t> out(c.size());
        for (std::size_t i = 0; i < c.size(); ++i)
            out[i] = c[i].val();
        return out;
    };
    auto c1 = run.template operator()<mint754974721>();
    auto c2 = run.template operator()<mint167772161>();
    auto c3 = run.template operator()<mint469762049>();
    return detail::convolution_crt_ll(c1, c2, c3);
}
}
