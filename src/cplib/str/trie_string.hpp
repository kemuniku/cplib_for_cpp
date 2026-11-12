#pragma once
#include <cplib/common.hpp>
#include <memory>
#include <stdexcept>
#include <unordered_map>
#include <ostream>

namespace cplib {
namespace detail {
template <class T, unsigned char Low, unsigned char High, bool Array> struct TrieStringStorage {
    using Children = std::conditional_t<Array, std::array<Int, High - Low + 1>, std::unordered_map<T, Int>>;
    struct Entry {
        Int parent = 0, jump = 0, depth = 0;
        T character{};
        Children children{};
    };
    std::vector<Entry> nodes{1};
    Int ancestor(Int node, Int depth) const {
        while (nodes[node].depth > depth) {
            Int jump = nodes[node].jump;
            node = nodes[jump].depth >= depth ? jump : nodes[node].parent;
        }
        return node;
    }
    Int lca(Int a, Int b) const {
        Int depth = std::min(nodes[a].depth, nodes[b].depth);
        a = ancestor(a, depth);
        b = ancestor(b, depth);
        if (a == b)
            return a;
        while (nodes[a].parent != nodes[b].parent) {
            Int ja = nodes[a].jump, jb = nodes[b].jump;
            if (ja != jb) {
                a = ja;
                b = jb;
            } else {
                a = nodes[a].parent;
                b = nodes[b].parent;
            }
        }
        return nodes[a].parent;
    }
};
}
template <class T, unsigned char Low, unsigned char High, bool Array> class BasicTrieString;
template <class T, unsigned char Low, unsigned char High, bool Array> class BasicTrieStringNode {
    using Storage = detail::TrieStringStorage<T, Low, High, Array>;
    std::shared_ptr<Storage> owner_;
    Int id_ = 0;
    friend class BasicTrieString<T, Low, High, Array>;
    BasicTrieStringNode(std::shared_ptr<Storage> owner, Int id) : owner_(std::move(owner)), id_(id) {
    }
    void check() const {
        if (!owner_)
            throw std::invalid_argument("未初期化のTrieString節点です");
    }
    void compatible(const BasicTrieStringNode &other) const {
        check();
        other.check();
        if (owner_ != other.owner_)
            throw std::invalid_argument("異なるTrieStringの節点です");
    }
public:
    // 未初期化のハンドルを作る。使用前にarena.root()等を代入する。計算量 O(1)。
    BasicTrieStringNode() = default;
    // arena内の不変な節点番号を返す。根は0。計算量 O(1)。
    Int id() const {
        check();
        return id_;
    }
    // 表す列の長さを返す。計算量 O(1)。
    Int size() const {
        check();
        return owner_->nodes[id_].depth;
    }
    // 表す列の長さを返す。計算量 O(1)。
    Int length() const {
        return size();
    }
    // 一要素追加した履歴を返す。元は変更しない。Table版は期待償却 O(1)、Array版は償却 O(σ)。
    BasicTrieStringNode add(T c) const {
        check();
        Int existing = 0;
        if constexpr (Array) {
            unsigned char key = c;
            if (key < Low || key > High)
                throw std::invalid_argument("文字がアルファベット範囲外です");
            existing = owner_->nodes[id_].children[key - Low];
        } else {
            auto it = owner_->nodes[id_].children.find(c);
            if (it != owner_->nodes[id_].children.end())
                existing = it->second;
        }
        if (existing != 0)
            return {owner_, existing};
        Int j = owner_->nodes[id_].jump, jj = owner_->nodes[j].jump;
        Int jump = id_;
        if (owner_->nodes[id_].depth - owner_->nodes[j].depth == owner_->nodes[j].depth - owner_->nodes[jj].depth)
            jump = jj;
        Int index = owner_->nodes.size();
        typename Storage::Entry node;
        node.parent = id_;
        node.jump = jump;
        node.depth = size() + 1;
        node.character = c;
        owner_->nodes.push_back(std::move(node));
        if constexpr (Array)
            owner_->nodes[id_].children[static_cast<unsigned char>(c) - Low] = index;
        else
            owner_->nodes[id_].children[c] = index;
        return {owner_, index};
    }
    // 末尾を除いた履歴を返す。根では根を返す。計算量 O(1)。
    BasicTrieStringNode pop() const {
        check();
        return {owner_, owner_->nodes[id_].parent};
    }
    // 先頭length要素の履歴を返す。不正な長さはinvalid_argument。時間 O(log(N+1))。
    BasicTrieStringNode prefix(Int length) const {
        if (length < 0 || length > size())
            throw std::invalid_argument("接頭辞の長さが範囲外です");
        return {owner_, owner_->ancestor(id_, length)};
    }
    // 指定位置の要素を返す。範囲外はinvalid_argument。時間 O(log(N+1))。
    T operator[](Int index) const {
        if (index < 0 || index >= size())
            throw std::invalid_argument("添字が範囲外です");
        return owner_->nodes[owner_->ancestor(id_, index + 1)].character;
    }
    // 最長共通接頭辞の履歴を返す。同じarenaであること。時間 O(log(N+1))。
    BasicTrieStringNode lcpNode(const BasicTrieStringNode &other) const {
        compatible(other);
        return {owner_, owner_->lca(id_, other.id_)};
    }
    // 最長共通接頭辞の長さを返す。時間 O(log(N+1))。
    Int lcp(const BasicTrieStringNode &other) const {
        return lcpNode(other).size();
    }
    // 辞書順で-1,0,1を返す。ハッシュを使わず厳密に比較する。時間 O(log(N+1))。
    Int cmp(const BasicTrieStringNode &other) const {
        compatible(other);
        if (id_ == other.id_)
            return 0;
        Int common = owner_->nodes[owner_->lca(id_, other.id_)].depth;
        if (common == std::min(size(), other.size()))
            return size() < other.size() ? -1 : 1;
        T a = (*this)[common], b = other[common];
        if constexpr (std::is_same_v<T, char>)
            return static_cast<unsigned char>(a) < static_cast<unsigned char>(b) ? -1 : 1;
        else
            return a < b ? -1 : 1;
    }
    // 同一arenaの同一履歴かを返す。異なるarenaはinvalid_argument。計算量 O(1)。
    bool operator==(const BasicTrieStringNode &other) const {
        compatible(other);
        return id_ == other.id_;
    }
    // 辞書順の小なりを返す。時間 O(log(N+1))。
    bool operator<(const BasicTrieStringNode &other) const {
        return cmp(other) < 0;
    }
    // 辞書順の大なりを返す。時間 O(log(N+1))。
    bool operator>(const BasicTrieStringNode &other) const {
        return cmp(other) > 0;
    }
    // 辞書順の以下を返す。時間 O(log(N+1))。
    bool operator<=(const BasicTrieStringNode &other) const {
        return cmp(other) <= 0;
    }
    // 辞書順の以上を返す。時間 O(log(N+1))。
    bool operator>=(const BasicTrieStringNode &other) const {
        return cmp(other) >= 0;
    }
    // 列を復元する。時間・出力領域 O(L+1)。
    std::vector<T> toSeq() const {
        std::vector<T> result(size());
        Int v = id_;
        for (std::size_t i = result.size(); i-- > 0;) {
            result[i] = owner_->nodes[v].character;
            v = owner_->nodes[v].parent;
        }
        return result;
    }
    // byte列を文字列に復元する。時間・出力領域 O(L+1)。
    std::string to_string() const requires std::is_same_v<T, char> {
        auto values = toSeq();
        return {values.begin(), values.end()};
    }
    // byte列を出力する。時間・一時領域 O(L+1)。
    friend std::ostream &operator<<(std::ostream &out, const BasicTrieStringNode &node)
        requires std::is_same_v<T, char> {
        return out << node.to_string();
    }
};
template <class T, unsigned char Low, unsigned char High, bool Array> class BasicTrieString {
    static_assert(!Array || Low <= High);
    using Storage = detail::TrieStringStorage<T, Low, High, Array>;
    std::shared_ptr<Storage> storage_ = std::make_shared<Storage>();
public:
    using Node = BasicTrieStringNode<T, Low, High, Array>;
    // 空列の根を持つarenaを作る。Table版 O(1)、Array版 O(σ)。
    BasicTrieString() = default;
    // 根のハンドルを返す。ハンドルはarenaの寿命を共有する。計算量 O(1)。
    Node root() const {
        return {storage_, 0};
    }
    // 根を含む作成済み節点数を返す。計算量 O(1)。
    Int nodeCount() const {
        return storage_->nodes.size();
    }
};
template <class T> using TrieStringTable = BasicTrieString<T, 0, 0, false>;
template <class T> using TrieStringTableNode = BasicTrieStringNode<T, 0, 0, false>;
template <unsigned char Low = 'a', unsigned char High = 'z'>
using TrieStringArray = BasicTrieString<char, Low, High, true>;
template <unsigned char Low = 'a', unsigned char High = 'z'>
using TrieStringArrayNode = BasicTrieStringNode<char, Low, High, true>;
using TrieString = TrieStringTable<char>;
// 空のTable版arenaを作る。時間・領域 O(1)。
template <class T = char> auto initTrieStringTable() {
    return TrieStringTable<T>();
}
// 指定範囲のArray版arenaを作る。時間・領域 O(σ)。
template <unsigned char Low = 'a', unsigned char High = 'z'> auto initTrieStringArray() {
    return TrieStringArray<Low, High>();
}
// byte列用のarenaを作る。時間・領域 O(1)。
inline auto initTrieString() {
    return TrieString();
}
// byte列ハンドルを文字列に復元する。時間・出力領域 O(L+1)。
template <unsigned char L, unsigned char H, bool A>
auto to_string(const BasicTrieStringNode<char, L, H, A> &node) {
    return node.to_string();
}
}
