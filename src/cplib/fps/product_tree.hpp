#pragma once
#include <cplib/fps/formal_power_series.hpp>

namespace cplib {
template <Modint T> struct PolynomialProductTree {
    Int pointCount, leafCount = 1;
    std::vector<std::vector<T>> nodes;

    explicit PolynomialProductTree(const std::vector<T> &xs) : pointCount(xs.size()) {
        while (leafCount < pointCount)
            leafCount *= 2;
        nodes.resize(leafCount * 2);
        for (Int i = 0; i < leafCount; ++i)
            nodes[leafCount + i] =
                i < pointCount ? std::vector<T>{-xs[i], T(1)} : std::vector<T>{T(1)};
        for (Int i = leafCount - 1; i >= 1; --i)
            nodes[i] = nodes[i * 2] * nodes[i * 2 + 1];
    }

    // 剰余木を64点のブロックまで降り、末端は直接評価。O(M(N)log N)。
    // 木の構築に使ったすべての点でfを評価する。
    std::vector<T> evaluate(const std::vector<T> &f) const {
        std::vector<T> out(pointCount);
        if (!pointCount)
            return out;
        Int directLevel = leafCount, pointsPerNode = 1;
        while (pointsPerNode < 64 && directLevel > 1) {
            directLevel /= 2;
            pointsPerNode *= 2;
        }
        std::vector<std::vector<T>> rem(nodes.size());
        rem[1] = f.size() < nodes[1].size() ? f : f % nodes[1];
        for (Int i = 1; i < directLevel; ++i) {
            rem[i * 2] = rem[i] % nodes[i * 2];
            rem[i * 2 + 1] = rem[i] % nodes[i * 2 + 1];
        }
        for (Int node = directLevel; node < directLevel * 2; ++node) {
            Int first = (node - directLevel) * pointsPerNode,
                last = std::min(first + pointsPerNode, pointCount);
            for (Int i = first; i < last; ++i)
                out[i] = eval(rem[node], -nodes[leafCount + i][0]);
        }
        return out;
    }
};

template <Modint T> auto initPolynomialProductTree(const std::vector<T> &xs) {
    return PolynomialProductTree<T>(xs);
}

template <Modint T> auto evaluate(const PolynomialProductTree<T> &tree, const std::vector<T> &f) {
    return tree.evaluate(f);
}

// 反転積木・中間積版。16点ブロックで直接評価する独立した実装。O(M(N)log N)。
// 反転積木と中間積を用いてfをすべての点で評価する。
template <Modint T>
std::vector<T> multipointEvaluation(const std::vector<T> &f, const std::vector<T> &xs) {
    if (xs.empty())
        return {};
    Int leafCount = 1;
    while (leafCount < Int(xs.size()))
        leafCount *= 2;
    Int blockSize = std::min(Int(16), leafCount), blockCount = leafCount / blockSize;
    std::vector<std::vector<T>> reversedProducts(blockCount * 2);
    for (Int block = 0; block < blockCount; ++block) {
        std::vector<T> product = {T(1)};
        for (Int offset = 0; offset < blockSize; ++offset) {
            Int point = block * blockSize + offset;
            T x = point < Int(xs.size()) ? xs[point] : T(0);
            product.push_back(T(0));
            for (Int i = Int(product.size()) - 1; i >= 1; --i)
                product[i] -= x * product[i - 1];
        }
        reversedProducts[blockCount + block] = std::move(product);
    }
    for (Int node = blockCount - 1; node >= 1; --node)
        reversedProducts[node] = reversedProducts[node * 2] * reversedProducts[node * 2 + 1];
    auto polynomial = f;
    if (Int(polynomial.size()) > leafCount) {
        auto root = reversedProducts[1];
        std::reverse(root.begin(), root.end());
        polynomial = polynomial % root;
    }
    std::vector<T> reversedPolynomial(leafCount);
    for (Int i = 0; i < Int(polynomial.size()); ++i)
        reversedPolynomial[leafCount - 1 - i] = polynomial[i];
    std::vector<std::vector<T>> transformed(blockCount * 2);
    transformed[1] = prefix(reversedPolynomial * inv(reversedProducts[1], leafCount), leafCount);
    for (Int node = 1; node < blockCount; ++node) {
        Int childSize = transformed[node].size() / 2;
        auto left = convolutionCyclicPowerOfTwo(transformed[node], reversedProducts[node * 2 + 1],
                                                transformed[node].size()),
             right = convolutionCyclicPowerOfTwo(transformed[node], reversedProducts[node * 2],
                                                 transformed[node].size());
        transformed[node * 2] = {left.begin() + childSize, left.begin() + childSize * 2};
        transformed[node * 2 + 1] = {right.begin() + childSize, right.begin() + childSize * 2};
    }
    std::vector<T> out(xs.size());
    for (Int block = 0; block < blockCount; ++block) {
        const auto &tb = transformed[blockCount + block];
        const auto &rb = reversedProducts[blockCount + block];
        std::vector<T> rem(blockSize);
        for (Int i = 0; i < blockSize; ++i)
            for (Int j = 0; j <= i; ++j)
                rem[i] += tb[j] * rb[i - j];
        std::reverse(rem.begin(), rem.end());
        Int first = block * blockSize, last = std::min(first + blockSize, Int(xs.size()));
        for (Int i = first; i < last; ++i)
            out[i] = eval(rem, xs[i]);
    }
    return out;
}
}
