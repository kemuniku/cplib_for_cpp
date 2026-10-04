#pragma once
#include <cplib/graph/graph.hpp>

namespace cplib {
template <char First = 'a', char Last = 'z'> struct TrieNode {
    static_assert(static_cast<unsigned char>(First) <= static_cast<unsigned char>(Last));
    std::array<std::int32_t,
               static_cast<unsigned char>(Last) - static_cast<unsigned char>(First) + 1>
        child{};
    std::int32_t terminal = 0, subtree = 0, parent = -1;
    char character = 0;
};

template <char First = 'a', char Last = 'z'> class Trie {
    static int code(char c) {
        return static_cast<unsigned char>(c) - static_cast<unsigned char>(First);
    }

    static bool valid(char c) {
        return static_cast<unsigned char>(First) <= static_cast<unsigned char>(c) &&
               static_cast<unsigned char>(c) <= static_cast<unsigned char>(Last);
    }

public:
    std::vector<TrieNode<First, Last>> nodes;

    // 指定した文字範囲の空の多重集合を作る。文字種数をσとして O(σ)。
    Trie() : nodes(1) {
    }

    // 文字列列から重複を保って構築する。O(σ(1 + Σ|s|) + words.size())。
    explicit Trie(std::span<const std::string> words) : Trie() {
        for (const auto &s : words)
            incl(s);
    }

    // 根の節点ID（0）を返す。O(1)。
    Int root() const {
        return 0;
    }

    // 必要なら子を作る。登録個数は変えない。償却O(文字種数)。
    Int getChild(Int node, char c) {
        assert(valid(c));
        if (nodes.empty()) {
            assert(node == 0);
            nodes.resize(1);
        }
        assert(node >= 0 && node < Int(nodes.size()));
        Int result = nodes[node].child[code(c)];
        if (!result) {
            result = nodes.size();
            assert(result <= std::numeric_limits<std::int32_t>::max());
            nodes.emplace_back();
            nodes[result].parent = node;
            nodes[result].character = c;
            nodes[node].child[code(c)] = result;
        }
        return result;
    }

    // 親の節点IDを返す。根の親は -1。O(1)。
    Int getParent(Int node) const {
        assert(node >= 0 && node < Int(nodes.size()));
        return nodes[node].parent;
    }

    // 節点が示す文字列を復元する。根なら空文字列。O(深さ + 1)。
    std::string restoreString(Int node) const {
        assert(node >= 0 && node < Int(nodes.size()));
        std::string result;
        for (; node; node = nodes[node].parent)
            result += nodes[node].character;
        std::reverse(result.begin(), result.end());
        return result;
    }

    // 重複を含む文字列の総数を返す。O(1)。
    Int len() const {
        return nodes.empty() ? 0 : nodes[0].subtree;
    }

    Int size() const {
        return len();
    }

    // sに対応する節点IDを個数が0でも返す。存在しなければ -1。O(|s|)。
    Int findNode(std::string_view s) const {
        if (nodes.empty())
            return -1;
        Int node = 0;
        for (char c : s) {
            if (!valid(c))
                return -1;
            node = nodes[node].child[code(c)];
            if (!node)
                return -1;
        }
        return node;
    }

    // sと完全一致する文字列の個数を返す。O(|s|)。
    Int count(std::string_view s) const {
        Int node = findNode(s);
        return node < 0 ? 0 : nodes[node].terminal;
    }

    // sが1個以上含まれるかを返す。O(|s|)。
    bool contains(std::string_view s) const {
        return count(s) > 0;
    }

    // sで始まる文字列の個数を返す。空のsなら総数。O(|s|)。
    Int countPrefix(std::string_view s) const {
        Int node = findNode(s);
        return node < 0 ? 0 : nodes[node].subtree;
    }

    // 文字列をv個追加する。O(σ|s|+1)。
    void incl(std::string_view s, Int v = 1) {
        assert(v >= 0);
        if (!v)
            return;
        for (char c : s) {
            assert(valid(c));
            (void)c;
        }
        assert(v <= std::numeric_limits<std::int32_t>::max() - len());
        if (nodes.empty())
            nodes.resize(1);
        Int node = 0;
        nodes[0].subtree += v;
        for (char c : s) {
            node = getChild(node, c);
            nodes[node].subtree += v;
        }
        nodes[node].terminal += v;
    }

    // 最大v個削除する。節点は保持する。O(|s|)。
    void excl(std::string_view s, Int v = 1) {
        assert(v >= 0);
        Int removed = std::min(v, count(s));
        if (!removed)
            return;
        Int node = 0;
        nodes[0].subtree -= removed;
        for (char c : s) {
            node = nodes[node].child[code(c)];
            nodes[node].subtree -= removed;
        }
        nodes[node].terminal -= removed;
    }

    // 辞書順でs未満の登録個数を返す。O(σ|s|+1)。
    Int lowerBound(std::string_view s) const {
        if (nodes.empty())
            return 0;
        Int result = 0, node = 0;
        for (char c : s) {
            result += nodes[node].terminal;
            for (int smaller = static_cast<unsigned char>(First);
                 smaller <= static_cast<unsigned char>(Last) &&
                 smaller < static_cast<unsigned char>(c);
                 ++smaller) {
                Int next = nodes[node].child[smaller - static_cast<unsigned char>(First)];
                if (next)
                    result += nodes[next].subtree;
            }
            if (!valid(c))
                return result;
            node = nodes[node].child[code(c)];
            if (!node)
                return result;
        }
        return result;
    }

    // 辞書順でs以下の文字列の個数を返す。O(σ|s| + 1)。
    Int upperBound(std::string_view s) const {
        return lowerBound(s) + count(s);
    }

    // 節点IDを頂点番号とし、親から子へ文字を重みとする辺を張る。全節点を含め O(N + 1)。
    auto toGraph() const {
        auto g = initWeightedDirectedGraph<char>(std::max<std::size_t>(1, nodes.size()));
        for (Int node = 1; node < Int(nodes.size()); ++node)
            g.add_edge(nodes[node].parent, node, nodes[node].character);
        return g;
    }
};

template <char First = 'a', char Last = 'z'> class TriePointer {
    Trie<First, Last> *owner = nullptr;
    Int id = 0;

public:
    TriePointer() = default;

    // 指定した節点を指すポインタを作る。省略時は根。σ固定で O(1)。
    explicit TriePointer(Trie<First, Last> &trie, Int node = 0) : owner(&trie), id(node) {
        if (trie.nodes.empty()) {
            assert(node == 0);
            trie = Trie<First, Last>();
        }
        assert(node >= 0 && node < Int(trie.nodes.size()));
    }

    // 指している節点IDを返す。O(1)。
    Int nodeId() const {
        assert(owner);
        return id;
    }

    // 指している文字列をprefixとする登録個数を返す。O(1)。
    Int subtree() const {
        assert(owner);
        return owner->nodes[id].subtree;
    }

    // 指している文字列と完全一致する登録個数を返す。O(1)。
    Int terminal() const {
        assert(owner);
        return owner->nodes[id].terminal;
    }

    // 指している節点の文字列を復元する。O(深さ + 1)。
    std::string restoreString() const {
        assert(owner);
        return owner->restoreString(id);
    }

    // cを末尾に追加した位置を返し、なければ節点を作る。σ固定で償却 O(1)。
    TriePointer getChild(char c) const {
        assert(owner);
        return TriePointer(*owner, owner->getChild(id, c));
    }

    // 末尾を1文字削除した位置を返す。根では呼べない。O(1)。
    TriePointer getParent() const {
        assert(owner && id != 0);
        return TriePointer(*owner, owner->getParent(id));
    }

    // 末尾にcを追加した位置を返す。σ固定で償却 O(1)。
    TriePointer operator&(char c) const {
        return getChild(c);
    }

    // 末尾にcを追加した位置へ移動する。登録個数は変えない。σ固定で償却 O(1)。
    void add(char c) {
        *this = getChild(c);
    }

    // 末尾にcを追加した位置へ移動する。σ固定で償却 O(1)。
    TriePointer &operator&=(char c) {
        add(c);
        return *this;
    }

    // 末尾の文字を返して親へ移動する。空文字列では呼べない。O(1)。
    char pop() {
        assert(owner && id != 0);
        char c = owner->nodes[id].character;
        *this = getParent();
        return c;
    }
};

template <char First = 'a', char Last = 'z'> Trie<First, Last> initTrie() {
    return {};
}

template <char First = 'a', char Last = 'z'>
Trie<First, Last> initTrie(std::span<const std::string> words) {
    return Trie<First, Last>(words);
}

template <char First, char Last> auto initTriePointer(Trie<First, Last> &trie, Int node = 0) {
    return TriePointer<First, Last>(trie, node);
}

template <char First, char Last> std::string to_string(const TriePointer<First, Last> &p) {
    return p.restoreString();
}
}
