#pragma once
#include <cplib/common.hpp>
#include <memory>
#include <sstream>

namespace cplib::detail {
struct BinaryTrieNode {
    std::shared_ptr<BinaryTrieNode> zero, one;
    Int value = 0;
};

class BinaryTrieBase {
protected:
    std::shared_ptr<BinaryTrieNode> root;
    Int h;

    explicit BinaryTrieBase(Int h) : root(std::make_shared<BinaryTrieNode>()), h(h) {
        assert(0 <= h && h <= 64);
    }

    static auto clone(const std::shared_ptr<BinaryTrieNode> &node)
        -> std::shared_ptr<BinaryTrieNode> {
        if (!node)
            return {};
        auto out = std::make_shared<BinaryTrieNode>(*node);
        out->zero = clone(node->zero);
        out->one = clone(node->one);
        return out;
    }

    Int bound(Int x, Int xor_value, bool inclusive) const {
        assert(x >= 0);
        auto now = root;
        Int out = 0;
        for (Int i = h - 1; i >= 0; --i) {
            auto zero = ((UInt(xor_value) >> i) & 1) ? now->one : now->zero,
                 one = ((UInt(xor_value) >> i) & 1) ? now->zero : now->one;
            if (!((UInt(x) >> i) & 1)) {
                if (!zero)
                    return out;
                now = zero;
            } else {
                if (zero)
                    out += zero->value;
                if (!one)
                    return out;
                now = one;
            }
        }
        return out + (inclusive ? now->value : 0);
    }

public:
    Int len() const {
        return root->value;
    }

    Int count(Int x) const {
        assert(x >= 0);
        auto now = root;
        for (Int i = h - 1; i >= 0; --i) {
            now = ((UInt(x) >> i) & 1) ? now->one : now->zero;
            if (!now)
                return 0;
        }
        return now->value;
    }

    bool contains(Int x) const {
        return count(x) != 0;
    }

    // XOR後の値で順序を付け、返す値は元の値。存在しなければ-1。O(H)。
    Int get_kth(Int k, Int xor_value = 0) const {
        assert(k >= 0);
        auto now = root;
        if (now->value <= k)
            return -1;
        Int cnt = 0;
        UInt out = 0;
        for (Int i = h - 1; i >= 0; --i) {
            if (!now->zero) {
                now = now->one;
                out = (out << 1) | 1;
            } else if (!now->one) {
                now = now->zero;
                out <<= 1;
            } else if (!((UInt(xor_value) >> i) & 1)) {
                if (cnt + now->zero->value > k) {
                    now = now->zero;
                    out <<= 1;
                } else {
                    cnt += now->zero->value;
                    now = now->one;
                    out = (out << 1) | 1;
                }
            } else {
                if (cnt + now->one->value > k) {
                    now = now->one;
                    out = (out << 1) | 1;
                } else {
                    cnt += now->one->value;
                    now = now->zero;
                    out <<= 1;
                }
            }
        }
        return std::bit_cast<Int>(out);
    }

    Int operator[](Int k) const {
        return get_kth(k);
    }

    std::vector<std::pair<Int, Int>> RLE() const {
        std::vector<std::pair<Int, Int>> out;
        auto dfs = [&](auto &&self, const std::shared_ptr<BinaryTrieNode> &node,
                       UInt value) -> void {
            if (node->zero)
                self(self, node->zero, value << 1);
            if (node->one)
                self(self, node->one, (value << 1) | 1);
            if (!node->zero && !node->one && node->value != 0)
                out.emplace_back(std::bit_cast<Int>(value), node->value);
        };
        dfs(dfs, root, 0);
        return out;
    }

    std::string to_string() const {
        std::ostringstream out;
        out << "@[";
        bool first = true;
        for (auto [value, count] : RLE()) {
            if (!first)
                out << ", ";
            first = false;
            out << '(' << value << ", " << count << ')';
        }
        return out << ']', out.str();
    }
};
}
