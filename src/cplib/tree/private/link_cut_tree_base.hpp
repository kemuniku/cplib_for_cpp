#pragma once
#include <cplib/utils/backwards_index.hpp>
#include <memory>
#include <functional>

namespace cplib::detail {
// 集約・遅延情報の処理を派生型へ委譲する共通splay操作。各公開操作は償却O(log N)。
template <class Derived, class S> struct LinkCutTreeOperations {
    Derived &self() {
        return static_cast<Derived &>(*this);
    }

    const Derived &self() const {
        return static_cast<const Derived &>(*this);
    }

    bool isAuxRoot(Int v) const {
        auto &n = self().nodes;
        Int p = n[v].parent;
        return !p || (n[p].left != v && n[p].right != v);
    }

    void toggle(Int v) {
        if (!v)
            return;
        auto &x = self().nodes[v];
        std::swap(x.left, x.right);
        std::swap(x.prod, x.rprod);
        x.rev = !x.rev;
    }

    void rotate(Int v) {
        auto &s = self();
        auto &n = s.nodes;
        Int p = n[v].parent, g = n[p].parent;
        bool right = n[p].right == v;
        Int middle = right ? n[v].left : n[v].right;
        if (!isAuxRoot(p)) {
            if (n[g].left == p)
                n[g].left = v;
            else
                n[g].right = v;
        }
        s.setParent(middle, p);
        s.setParent(v, g);
        s.setParent(p, v);
        if (right) {
            n[p].right = middle;
            n[v].left = p;
        } else {
            n[p].left = middle;
            n[v].right = p;
        }
        s.pull(p);
        s.pull(v);
    }

    void splay(Int v) {
        auto &s = self();
        auto &n = s.nodes;
        s.stack.clear();
        Int x = v;
        s.stack.push_back(x);
        while (!isAuxRoot(x)) {
            x = n[x].parent;
            s.stack.push_back(x);
        }
        for (auto it = s.stack.rbegin(); it != s.stack.rend(); ++it)
            s.push(*it);
        while (!isAuxRoot(v)) {
            Int p = n[v].parent, g = n[p].parent;
            if (!isAuxRoot(p)) {
                if ((n[p].left == v) == (n[g].left == p))
                    rotate(p);
                else
                    rotate(v);
            }
            rotate(v);
        }
    }

    void accessNode(Int v) {
        auto &s = self();
        Int last = 0, x = v;
        while (x) {
            splay(x);
            s.addVirtual(x, s.nodes[x].right);
            s.removeVirtual(x, last);
            s.nodes[x].right = last;
            s.pull(x);
            last = x;
            x = s.nodes[x].parent;
        }
        splay(v);
    }

    Int len() const {
        return Int(self().nodes.size()) - 1;
    }

    void check(Int v) const {
        assert(v >= 0 && v < len());
        (void)v;
    }

    // vを所属する木の根にする。償却O(log N)。
    void makeRoot(Int v) {
        check(v);
        accessNode(v + 1);
        toggle(v + 1);
    }

    // vが所属する木の現在の根を返す。償却O(log N)。
    Int findRoot(Int v) {
        check(v);
        auto &s = self();
        Int x = v + 1;
        accessNode(x);
        s.push(x);
        while (s.nodes[x].left) {
            x = s.nodes[x].left;
            s.push(x);
        }
        splay(x);
        return x - 1;
    }

    // uとvが同じ木に属するかを返す。償却O(log N)。
    bool connected(Int u, Int v) {
        check(u);
        check(v);
        return u == v || findRoot(u) == findRoot(v);
    }

    // 結合後の根は、結合前のv側の根になる。
    // 異なる木の頂点u, vを辺で結ぶ。償却O(log N)。
    void link(Int u, Int v) {
        check(u);
        check(v);
        auto &s = self();
        makeRoot(u);
        assert(findRoot(v) != u);
        accessNode(v + 1);
        s.setParent(u + 1, v + 1);
        s.addVirtual(v + 1, u + 1);
        s.pull(v + 1);
    }

    // 切断後の根はそれぞれuとvになる。
    // 存在する辺(u, v)を削除する。償却O(log N)。
    void cut(Int u, Int v) {
        check(u);
        check(v);
        auto &s = self();
        makeRoot(u);
        accessNode(v + 1);
        s.push(u + 1);
        assert(s.nodes[v + 1].left == u + 1 && s.nodes[u + 1].right == 0);
        s.setParent(u + 1, 0);
        s.nodes[v + 1].left = 0;
        s.pull(v + 1);
    }

    void update(Int v, const S &value) {
        check(v);
        accessNode(v + 1);
        self().nodes[v + 1].value = value;
        self().pull(v + 1);
    }

    S point(Int v) {
        check(v);
        accessNode(v + 1);
        return self().nodes[v + 1].value;
    }

    // 同じ木のuからvへのパスを両端込みで順に集約する。償却O(log N)。根をuに変更する。
    S pathProd(Int u, Int v) {
        check(u);
        check(v);
        makeRoot(u);
        accessNode(v + 1);
        return self().nodes[v + 1].prod;
    }

    // vを含む木全体を集約する。可換群と逆元の指定が必要。償却O(log N)。
    S componentProd(Int v) {
        check(v);
        assert(self().inverse);
        accessNode(v + 1);
        return self().nodes[v + 1].all;
    }

    // 辺(v, parent)のv側を集約する。可換群と逆元の指定が必要。償却O(log N)。根をparentに変更する。
    S subtreeProd(Int v, Int parent) {
        check(v);
        check(parent);
        auto &s = self();
        assert(s.inverse);
        makeRoot(parent);
        accessNode(v + 1);
        s.push(parent + 1);
        assert(s.nodes[v + 1].left == parent + 1 && s.nodes[parent + 1].right == 0);
        return s.merge(s.nodes[v + 1].value, s.nodes[v + 1].virtualValue);
    }
};

// Nimのref objectと同様に、コピーしたハンドルは森を共有する。
template <class Impl, class S> class LinkCutTreeHandle {
protected:
    std::shared_ptr<Impl> impl_;

    explicit LinkCutTreeHandle(std::shared_ptr<Impl> p) : impl_(std::move(p)) {
    }

public:
    using value_type = S;

    Int len() const {
        return impl_->len();
    }

    void makeRoot(Int v) const {
        impl_->makeRoot(v);
    }

    Int findRoot(Int v) const {
        return impl_->findRoot(v);
    }

    bool connected(Int u, Int v) const {
        return impl_->connected(u, v);
    }

    void link(Int u, Int v) const {
        impl_->link(u, v);
    }

    void cut(Int u, Int v) const {
        impl_->cut(u, v);
    }

    void update(Int v, const S &value) const {
        impl_->update(v, value);
    }

    S get(Int v) const {
        return impl_->point(v);
    }

    S get(Int u, Int v) const {
        return impl_->pathProd(u, v);
    }

    S pathProd(Int u, Int v) const {
        return impl_->pathProd(u, v);
    }

    S componentProd(Int v) const {
        return impl_->componentProd(v);
    }

    S subtreeProd(Int v, Int parent) const {
        return impl_->subtreeProd(v, parent);
    }

    struct Reference {
        const LinkCutTreeHandle *owner;
        Int vertex;

        operator S() const {
            return owner->get(vertex);
        }

        Reference &operator=(const S &value) {
            owner->update(vertex, value);
            return *this;
        }

        Reference &operator=(const Reference &b) {
            return *this = S(b);
        }
    };

    Reference operator[](Int v) {
        impl_->check(v);
        return {this, v};
    }

    S operator[](Int v) const {
        return get(v);
    }

    CPLIB_BACKWARDS_INDEX_OVERLOADS
};
}
