#pragma once
#include <cplib/utils/constants.hpp>
#include <bit>

namespace cplib {
// 座標圧縮 Li Chao 木。直線追加・点取得 O(log N)、線分追加 O(log² N)。
class LiChaoTree {
    std::vector<Int> x_, a_, b_;
    std::size_t last_ = 1;

    void query(Int a, Int b, std::size_t now) {
        auto depth = std::bit_width(now) - 1;
        std::size_t length = last_ >> depth, i = now ^ (std::size_t(1) << depth), l = i * length,
                    r = (i + 1) * length;
        for (;;) {
            auto m = (l + r) / 2;
            bool left = x_[l] * a + b < x_[l] * a_[now] + b_[now],
                 right = x_[r - 1] * a + b < x_[r - 1] * a_[now] + b_[now],
                 mid = x_[m] * a + b < x_[m] * a_[now] + b_[now];
            if (left && right) {
                a_[now] = a;
                b_[now] = b;
                return;
            }
            if (!left && !right)
                return;
            if (mid) {
                std::swap(a, a_[now]);
                std::swap(b, b_[now]);
            }
            if (left != mid) {
                now *= 2;
                r = m;
            } else {
                now = now * 2 + 1;
                l = m;
            }
        }
    }

public:
    // LiChaoTreeを初期化します
    // get_minで用いる可能性のあるXを配列で与えてください。
    explicit LiChaoTree(std::vector<Int> x) {
        std::sort(x.begin(), x.end());
        x.erase(std::unique(x.begin(), x.end()), x.end());
        while (last_ < x.size())
            last_ *= 2;
        x_.assign(last_ + 1, INF32);
        std::copy(x.begin(), x.end(), x_.begin());
        a_.assign(2 * last_, 0);
        b_.assign(2 * last_, INF64);
    }

    // 直線ax+bを追加します。
    // aは32bit整数に収まるようにしてください。
    void add_line(Int a, Int b) {
        query(a, b, 1);
    }

    // 線分ax+b (l<=x<r)を追加します。
    // aは32bit整数に収まるようにしてください。
    void add_segment(Int a, Int b, Int l, Int r) {
        std::size_t ql = std::lower_bound(x_.begin(), x_.end(), l) - x_.begin() + last_,
                    qr = std::lower_bound(x_.begin(), x_.end(), r) - x_.begin() + last_;
        while (ql < qr) {
            if (ql & 1)
                query(a, b, ql++);
            if (qr & 1)
                query(a, b, --qr);
            ql /= 2;
            qr /= 2;
        }
    }

    // xにおける最小値を返します。
    // xは32bit整数に収まるようにしてください。
    // 線分が存在しない場合、INF64が返ります。
    // xは初期化時に与える必要があります。
    Int get_min(Int x) const {
        auto it = std::lower_bound(x_.begin(), x_.end(), x);
        assert(it != x_.end() && *it == x);
        std::size_t now = it - x_.begin() + last_;
        Int result = INF64;
        while (now) {
            result = std::min(result, a_[now] * x + b_[now]);
            now /= 2;
        }
        return result;
    }
};

inline LiChaoTree initLiChaoTree(const std::vector<Int> &x) {
    return LiChaoTree(x);
}
}
