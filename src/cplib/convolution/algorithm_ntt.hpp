#pragma once
#include <cplib/convolution/convolution.hpp>
#include <memory>

namespace cplib {
// 一連の畳み込みで正逆NTTの変換表を共有する。スペクトルはMontgomery表現。
class AlgorithmNtt {
    using Plan = detail::avx2_ntt::TransformPlan;
    std::unique_ptr<Plan> forward_, inverse_;
    std::uint32_t modulus_, root_;
    Int size_;
    void activate() const {
        assert(forward_ && inverse_);
        detail::avx2_ntt::modulus = modulus_;
        detail::avx2_ntt::primitive_root = root_;
    }
public:
    // 64以上の2冪sizeに対応するNTT素数で初期化する。時間・領域 O(size)。
    AlgorithmNtt(std::uint32_t modulus, Int size) : modulus_(modulus), size_(size) {
        assert(size >= 64 && std::has_single_bit(UInt(size)) && modulus > 1 &&
               UInt(size) <= modulus - 1 && detail::isNttFriendlyModulus(modulus, size));
        detail::avx2_ntt::modulus = modulus;
        root_ = detail::avx2_ntt::find_primitive_root(modulus);
        detail::avx2_ntt::primitive_root = root_;
        forward_ = std::make_unique<Plan>(size);
        inverse_ = std::make_unique<Plan>(size);
        inverse_->prepare_inverse();
    }
    // 変換長を返す。計算量 O(1)。
    Int size() const {
        return size_;
    }
    // 所有する変換表を解放する。繰り返し可能。計算量 O(size)。
    void close() {
        forward_.reset();
        inverse_.reset();
    }
    // fを変更せずMontgomery表現のスペクトルを返す。時間 O(size log size)、領域 O(size)。
    template <class T> std::vector<std::uint32_t> spectrum(const std::vector<T> &f) const {
        assert(f.size() <= std::size_t(size_) && T::umod() == modulus_);
        activate();
        const auto &mont = forward_->montgomery();
        std::vector<std::uint32_t> result(size_);
        for (std::size_t i = 0; i < f.size(); ++i)
            result[i] = mont.to_montgomery(f[i].val());
        forward_->forward(result.data());
        return result;
    }
    // 同じ法・長さのスペクトルの積をoutputに加算する。時間 O(size)、追加領域 O(1)。
    void addSpectrumProduct(std::vector<std::uint32_t> &output,
                            const std::vector<std::uint32_t> &left,
                            const std::vector<std::uint32_t> &right) const {
        assert(output.size() == std::size_t(size_) && left.size() == output.size() && right.size() == output.size());
        activate();
        const auto &mont = forward_->montgomery();
        for (Int i = 0; i < size_; ++i)
            output[i] = detail::avx2_ntt::add_mod(output[i], mont.multiply(left[i], right[i]));
    }
    // 同じコンテキストのスペクトルの積を返す。時間・領域 O(size)。
    auto spectrumProduct(const std::vector<std::uint32_t> &left,
                         const std::vector<std::uint32_t> &right) const {
        std::vector<std::uint32_t> result(size_);
        addSpectrumProduct(result, left, right);
        return result;
    }
    // valuesを消費し[start,start+count)の係数を返す。時間 O(size log size)、出力領域 O(count)。
    template <class T>
    std::vector<T> coefficients(std::vector<std::uint32_t> &values, Int start, Int count) const {
        assert(values.size() == std::size_t(size_) && 0 <= start && 0 <= count && start <= size_ - count);
        assert(T::umod() == modulus_);
        activate();
        inverse_->inverse(values.data());
        const auto &mont = inverse_->montgomery();
        auto inverseSize = detail::avx2_ntt::power_mod(size_ % modulus_, modulus_ - 2);
        std::vector<T> result(count);
        for (Int i = 0; i < count; ++i)
            result[i] = T(UInt(mont.multiply(values[start + i], 1)) * inverseSize % modulus_);
        return result;
    }
};
// 共有変換表を構築する。時間・領域 O(size)。
inline auto initAlgorithmNtt(std::uint32_t modulus, Int size) {
    return AlgorithmNtt(modulus, size);
}
}
