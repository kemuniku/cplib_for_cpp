#pragma once
#include <cplib/common.hpp>
#include <cplib/utils/backwards_index.hpp>
#include <functional>
#include <stdexcept>
#include <sstream>

namespace cplib::detail {
template <class T> class SwagBase {
protected:
    std::function<T(const T &, const T &)> op;
    T e;
    std::vector<T> top, bottom, topfold, bottomfold;

    SwagBase(std::function<T(const T &, const T &)> operation, T identity)
        : op(std::move(operation)), e(identity), topfold{identity}, bottomfold{identity} {
    }

    void pushbottom(const T &x) {
        bottom.push_back(x);
        bottomfold.push_back(op(bottomfold.back(), x));
    }

    void pushtop(const T &x) {
        top.push_back(x);
        topfold.push_back(op(x, topfold.back()));
    }

    T popbottom() {
        bottomfold.pop_back();
        T x = std::move(bottom.back());
        bottom.pop_back();
        return x;
    }

    T poptop() {
        topfold.pop_back();
        T x = std::move(top.back());
        top.pop_back();
        return x;
    }

public:
    T fold() const {
        return op(topfold.back(), bottomfold.back());
    }

    Int len() const {
        return top.size() + bottom.size();
    }

    Int size() const {
        return len();
    }

    const T &operator[](Int i) const {
        if (i < 0 || i >= len())
            throw std::out_of_range("SWAG index");
        return i < Int(top.size()) ? top[top.size() - 1 - i] : bottom[i - top.size()];
    }

    CPLIB_BACKWARDS_INDEX_OVERLOADS
};

template <class T> std::string swag_string(const std::vector<T> &a) {
    std::ostringstream out;
    out << "@[";
    for (std::size_t i = 0; i < a.size(); ++i) {
        if (i)
            out << ", ";
        out << a[i];
    }
    out << ']';
    return out.str();
}
}
