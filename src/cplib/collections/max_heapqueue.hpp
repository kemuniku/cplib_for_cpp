#pragma once
#include <cplib/common.hpp>
#include <sstream>

namespace cplib {
template <class T> class MaxHeapQueue {
    std::vector<T> data;

    void up(Int i) {
        T x = std::move(data[i]);
        while (i && data[(i - 1) / 2] < x) {
            data[i] = std::move(data[(i - 1) / 2]);
            i = (i - 1) / 2;
        }
        data[i] = std::move(x);
    }

    void down(Int i) {
        T x = std::move(data[i]);
        for (Int c = i * 2 + 1; c < len(); c = i * 2 + 1) {
            if (c + 1 < len() && !(data[c + 1] < data[c]))
                ++c;
            if (!(x < data[c]))
                break;
            data[i] = std::move(data[c]);
            i = c;
        }
        data[i] = std::move(x);
    }

public:
    // 空のキューを作る。変数の宣言だけでも空に初期化される。O(1)。
    MaxHeapQueue() = default;

    // 配列の要素からキューを作る。O(N)。
    explicit MaxHeapQueue(std::span<const T> v) : data(v.begin(), v.end()) {
        for (Int i = len() / 2; i-- > 0;)
            down(i);
    }

    // 要素数を返す。O(1)。
    Int len() const {
        return data.size();
    }

    // 内部配列の i 番目を参照する。最大値は添字 0 で、全体は未整列。O(1)。
    const T &operator[](Int i) const {
        return data.at(i);
    }

    // 内部配列の順に全要素を列挙する。列挙中の要素数変更は不可。O(N)。
    auto begin() const {
        return data.begin();
    }

    auto end() const {
        return data.end();
    }

    // 要素を追加する。償却 O(log N)。
    void push(T x) {
        data.push_back(std::move(x));
        up(len() - 1);
    }

    // 最大値を取り除いて返す。空のキューには使用不可。O(log N)。
    T pop() {
        assert(len());
        T x = std::move(data[0]);
        T last = std::move(data.back());
        data.pop_back();
        if (len()) {
            data[0] = std::move(last);
            down(0);
        }
        return x;
    }

    // x と等しい最初の要素の添字を返す。存在しなければ -1。O(N)。
    Int find(const T &x) const {
        auto it = std::find(data.begin(), data.end(), x);
        return it == data.end() ? -1 : it - data.begin();
    }

    // x と等しい要素が存在するか返す。O(N)。
    bool contains(const T &x) const {
        return find(x) >= 0;
    }

    // 指定した添字の要素を削除する。O(log N)。
    void del(Int i) {
        assert(0 <= i && i < len());
        std::swap(data[i], data.back());
        data.pop_back();
        if (i < len()) {
            if (i && data[(i - 1) / 2] < data[i])
                up(i);
            else
                down(i);
        }
    }

    // 最大値を取り除いて返し、xを追加する。空には使用不可。O(log N)。
    T replace(T x) {
        assert(len());
        std::swap(x, data[0]);
        down(0);
        return x;
    }

    // xの追加後に最大値を取り除いて返す。空ならxを返す。O(log N)。
    T pushpop(T x) {
        if (len() && x < data[0]) {
            std::swap(x, data[0]);
            down(0);
        }
        return x;
    }

    // 全要素を削除する。要素の破棄を含めて O(N)。
    void clear() {
        data.clear();
    }

    // 内部配列の順に文字列化する。O(N + 出力文字列長)。
    std::string str() const {
        std::ostringstream s;
        s << '[';
        bool first = true;
        for (auto &x : data) {
            if (!first)
                s << ", ";
            first = false;
            s << x;
        }
        return s.str() + ']';
    }
};

template <class T> auto initMaxHeapQueue() {
    return MaxHeapQueue<T>();
}

template <class R> auto toMaxHeapQueue(const R &r) {
    using T = typename R::value_type;
    return MaxHeapQueue<T>(r);
}
}
