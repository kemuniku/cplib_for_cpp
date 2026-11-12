#pragma once
#include <cplib/common.hpp>

namespace cplib {
class Fenwick2D {
    Int n;
    std::vector<std::vector<Int>> ys, bit;

public:
    // 各x座標の更新候補を登録する。座標の総数をKとしてO(N+K log(N+1) log(K+2))時間。
    explicit Fenwick2D(const std::vector<std::vector<Int>> &posVals)
        : n(posVals.size()), ys(n + 1), bit(n + 1) {
        for (Int p = 0; p < n; ++p)
            for (Int i = p + 1; i <= n; i += i & -i)
                ys[i].insert(ys[i].end(), posVals[p].begin(), posVals[p].end());
        for (Int i = 1; i <= n; ++i) {
            std::sort(ys[i].begin(), ys[i].end());
            ys[i].erase(std::unique(ys[i].begin(), ys[i].end()), ys[i].end());
            bit[i].resize(ys[i].size() + 1);
        }
    }

    // 登録座標に加算する。O(log N log K)。
    void add(Int p, Int y, Int delta) {
        for (Int i = p + 1; i <= n; i += i & -i)
            for (Int k = std::lower_bound(ys[i].begin(), ys[i].end(), y) - ys[i].begin() + 1;
                 k < Int(bit[i].size()); k += k & -k)
                bit[i][k] += delta;
    }

    // x∈[0,r)、y<yUpperの和を返す。登録座標の総数をKとしてO(log(N+1) log(K+2))。
    Int prefix(Int r, Int yUpper) const {
        Int result = 0;
        for (Int i = r; i > 0; i -= i & -i)
            for (Int k = std::lower_bound(ys[i].begin(), ys[i].end(), yUpper) - ys[i].begin();
                 k > 0; k -= k & -k)
                result += bit[i][k];
        return result;
    }

    // x∈[l,r)、y<yUpperの和を返す。登録座標の総数をKとしてO(log(N+1) log(K+2))。
    Int getLess(Int l, Int r, Int yUpper) const {
        return prefix(r, yUpper) - prefix(l, yUpper);
    }

    // 登録座標の総数をKとしてO(log(N+1) log(K+2))。
    // 長方形領域 x in [l,r), y in [yLower,yUpper) の和を返します。
    Int get(Int l, Int r, Int yLower, Int yUpper) const {
        return getLess(l, r, yUpper) - getLess(l, r, yLower);
    }
};

// 更新候補座標から構築する。座標の総数をKとしてO(N+K log(N+1) log(K+2))。
inline Fenwick2D initFenwick2D(const std::vector<std::vector<Int>> &positions) {
    return Fenwick2D(positions);
}

// 登録済み座標(p,y)に加算する。座標総数をKとしてO(log(N+1) log(K+2))。
inline void add(Fenwick2D &f, Int p, Int y, Int delta) {
    f.add(p, y, delta);
}

// x∈[0,r)、座標がy未満の和を返す。O(log(N+1) log(K+2))。
inline Int prefix(const Fenwick2D &f, Int r, Int y) {
    return f.prefix(r, y);
}

// x∈[l,r)、座標がy未満の和を返す。O(log(N+1) log(K+2))。
inline Int getLess(const Fenwick2D &f, Int l, Int r, Int y) {
    return f.getLess(l, r, y);
}

// 長方形[l,r)×[yl,yr)の和を返す。O(log(N+1) log(K+2))。
inline Int get(const Fenwick2D &f, Int l, Int r, Int yl, Int yr) {
    return f.get(l, r, yl, yr);
}
}
