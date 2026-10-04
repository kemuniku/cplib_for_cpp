#pragma once
#include <cplib/common.hpp>
#include <memory>
#include <deque>

namespace cplib::detail {
template <class Node> struct PersistentArena {
    struct Owner {
        std::shared_ptr<std::deque<Node>> arena = std::make_shared<std::deque<Node>>();
        std::shared_ptr<Owner> parent;
        std::vector<std::shared_ptr<std::deque<Node>>> arenas{arena};
    };

    std::shared_ptr<Owner> owner = std::make_shared<Owner>();

    Node *make(Node node) const {
        owner->arena->push_back(std::move(node));
        return &owner->arena->back();
    }

    static std::shared_ptr<Owner> root(std::shared_ptr<Owner> p) {
        auto r = p;
        while (r->parent)
            r = r->parent;
        while (p->parent) {
            auto next = p->parent;
            p->parent = r;
            p = next;
        }
        return r;
    }

    void share(const PersistentArena &other) const {
        auto a = root(owner), b = root(other.owner);
        if (a == b)
            return;
        if (a->arenas.size() < b->arenas.size())
            std::swap(a, b);
        a->arenas.insert(a->arenas.end(), b->arenas.begin(), b->arenas.end());
        b->arenas.clear();
        b->parent = a;
    }
};
}
