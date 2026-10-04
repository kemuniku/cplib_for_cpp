#pragma once
#include <cplib/matrix/matrix_avx2_kernel.hpp>
#include <cplib/matrix/field_matrix_ops.hpp>
#include <cplib/modint/modint.hpp>
#include <span>
#include <cstring>

namespace cplib::detail::matrix_avx2 {
template <class T>
concept Element = MontgomeryModint<T> || BarrettModint<T>;

template <Element T> std::uint32_t modulus() {
    static_assert(sizeof(T) == 4 && alignof(T) == alignof(std::uint32_t) &&
                  std::is_trivially_copyable_v<T>);
    auto p = T::umod();
    assert(p > 0 && p < (1u << 30) && (p & 1));
    return p;
}

inline Int size(Int h, Int w) {
    assert(h >= 0 && w >= 0 && h <= INT32_MAX && w <= INT32_MAX);
    assert(!h || w <= (std::numeric_limits<Int>::max() / 4) / h);
    return h * w;
}

template <class T> const std::uint32_t *raw(const T *p) {
    return reinterpret_cast<const std::uint32_t *>(p);
}

template <class T> std::uint32_t *raw(T *p) {
    return reinterpret_cast<std::uint32_t *>(p);
}

template <Element T, class R> std::optional<LinearSystemSolution<T>> solution(R &r, Int h, Int w) {
    for (Int i = r.rank; i < h; ++i)
        if (r.values[i * r.width + w])
            return std::nullopt;
    fieldRestoreKernel(r.values.data(), h * r.width, modulus<T>(), MontgomeryModint<T>);
    LinearSystemSolution<T> out;
    out.particular.resize(w);
    std::vector<bool> pivot(w);
    for (Int i = 0; i < r.rank; ++i) {
        Int col = r.pivots[i];
        pivot[col] = true;
        out.particular[col] = std::bit_cast<T>(r.values[i * r.width + w]);
    }
    for (Int col = 0; col < w; ++col)
        if (!pivot[col]) {
            std::vector<T> v(w);
            v[col] = 1;
            for (Int i = 0; i < r.rank; ++i)
                v[r.pivots[i]] = -std::bit_cast<T>(r.values[i * r.width + col]);
            out.basis.push_back(std::move(v));
        }
    return out;
}

// 所有権を持つ公開値vectorの領域を再利用する。vector自体の型・内部レイアウトは変換しない。
// byte配列の寿命開始により暗黙的寿命型Tの配列領域を用意し、各要素を同じ位置で構築する。
template <Element T> struct AdoptAllocator {
    using value_type = T;
    std::shared_ptr<std::vector<std::uint32_t>> owner;
    AdoptAllocator() = default;

    explicit AdoptAllocator(std::shared_ptr<std::vector<std::uint32_t>> x) : owner(std::move(x)) {
    }

    template <class U> AdoptAllocator(const AdoptAllocator<U> &x) : owner(x.owner) {
    }

    T *allocate(std::size_t n) {
        assert(owner && n == owner->size());
        void *p = owner->data();
        ::new (p) std::byte[n * sizeof(T)];
        return reinterpret_cast<T *>(p);
    }

    void deallocate(T *, std::size_t) noexcept {
    }

    void construct(T *p) {
        std::uint32_t value;
        std::memcpy(&value, p, sizeof(value));
        std::construct_at(p, T(value));
    }

    void destroy(T *p) noexcept {
        std::destroy_at(p);
        ::new (static_cast<void *>(p)) std::uint32_t;
    }

    template <class U> bool operator==(const AdoptAllocator<U> &b) const noexcept {
        return owner == b.owner;
    }
};

template <Element T> struct Storage {
    std::shared_ptr<void> owner;
    T *pointer = nullptr;
    Int count = 0;

    explicit Storage(std::vector<T> values) {
        auto v = std::make_shared<std::vector<T>>(std::move(values));
        pointer = v->data();
        count = v->size();
        owner = std::move(v);
    }

    explicit Storage(std::vector<std::uint32_t> values) {
        auto words = std::make_shared<std::vector<std::uint32_t>>(std::move(values));
        auto v = std::make_shared<std::vector<T, AdoptAllocator<T>>>(words->size(),
                                                                     AdoptAllocator<T>(words));
        pointer = v->data();
        count = v->size();
        owner = std::move(v);
    }
};
}
