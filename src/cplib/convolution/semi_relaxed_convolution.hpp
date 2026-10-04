#pragma once
#include <cplib/convolution/convolution.hpp>

namespace cplib {
template <class T> class SemiRelaxedConvolution {
    std::vector<T> fixed, online, product;

public:
    // 左側をfixedに固定した逐次畳み込みを初期化する。
    // NTTが利用できる場合、N項追加する計算量はO(N log^2 N)となる。
    explicit SemiRelaxedConvolution(std::vector<T> fixed) : fixed(std::move(fixed)) {
    }

    Int len() const {
        return online.size();
    }

    Int size() const {
        return len();
    }

    // 現在までに確定した積の係数を返す。
    std::vector<T> coefficients() const {
        std::vector<T> out(online.size());
        std::copy_n(product.begin(), std::min(out.size(), product.size()), out.begin());
        return out;
    }

    // 低位ビットごとのブロックを加算。N項の追加はO(N log²N)。
    // オンライン側へ次の係数を追加し、その次数の積の係数を返す。
    T add(T value) {
        online.push_back(value);
        Int count = online.size();
        for (Int blockSize = 1;; blockSize *= 2) {
            Int fixedLeft = blockSize - 1,
                fixedRight = std::min(fixedLeft + blockSize, Int(fixed.size()));
            if (fixedLeft < fixedRight) {
                Int onlineLeft = count - blockSize, first = fixedLeft + onlineLeft;
                if (blockSize <= 16) {
                    Int required = fixedRight + count - 1;
                    if (Int(product.size()) < required)
                        product.resize(required);
                    for (Int i = fixedLeft; i < fixedRight; ++i)
                        for (Int j = onlineLeft; j < count; ++j)
                            product[i + j] += fixed[i] * online[j];
                } else {
                    auto values = convolution(
                        std::vector<T>(fixed.begin() + fixedLeft, fixed.begin() + fixedRight),
                        std::vector<T>(online.begin() + onlineLeft, online.begin() + count));
                    if (product.size() < first + values.size())
                        product.resize(first + values.size());
                    for (std::size_t i = 0; i < values.size(); ++i)
                        product[first + i] += values[i];
                }
            }
            if (blockSize == (count & -count))
                break;
        }
        return count - 1 < Int(product.size()) ? product[count - 1] : T{};
    }

    T append(T value) {
        return add(value);
    }

    T get(T value) {
        return add(value);
    }
};

template <class T> auto initSemiRelaxedConvolution(std::vector<T> fixed) {
    return SemiRelaxedConvolution<T>(std::move(fixed));
}

template <class T> Int len(const SemiRelaxedConvolution<T> &s) {
    return s.len();
}

template <class T> auto coefficients(const SemiRelaxedConvolution<T> &s) {
    return s.coefficients();
}

template <class T> T add(SemiRelaxedConvolution<T> &s, T value) {
    return s.add(value);
}

template <class T> T append(SemiRelaxedConvolution<T> &s, T value) {
    return s.add(value);
}

template <class T> T get(SemiRelaxedConvolution<T> &s, T value) {
    return s.add(value);
}
}
