#pragma once
#include <cplib/common.hpp>

namespace cplib {
template <class T> class DoubleEndedPriorityQueue {
    std::vector<T> data_;
    static bool minLevel(Int index) {
        return (std::bit_width(UInt(index + 1)) & 1) != 0;
    }
    static bool less(const T &a, const T &b, bool minimum) {
        return minimum ? a < b : b < a;
    }
    void bubbleUp(Int pos, bool minimum) {
        while (pos >= 3) {
            Int grandparent = (pos - 3) / 4;
            if (!less(data_[pos], data_[grandparent], minimum))
                break;
            std::swap(data_[pos], data_[grandparent]);
            pos = grandparent;
        }
    }
    void trickleDown(Int pos, bool minimum) {
        while (pos * 2 + 1 < size()) {
            Int child = pos * 2 + 1, grandchild = pos * 4 + 3, best = child;
            if (child + 1 < size() && less(data_[child + 1], data_[best], minimum))
                best = child + 1;
            for (Int next = grandchild; next < std::min(grandchild + 4, size()); ++next)
                if (less(data_[next], data_[best], minimum))
                    best = next;
            if (!less(data_[best], data_[pos], minimum))
                break;
            std::swap(data_[pos], data_[best]);
            if (best < grandchild)
                break;
            Int parent = (best - 1) / 2;
            if (less(data_[parent], data_[best], minimum))
                std::swap(data_[parent], data_[best]);
            pos = best;
        }
    }
    Int maxIndex() const {
        if (size() == 1)
            return 0;
        return size() == 2 || data_[2] < data_[1] ? 1 : 2;
    }
    T popAt(Int index, bool minimum) {
        T result = std::move(data_[index]);
        if (index + 1 == size()) {
            data_.pop_back();
        } else {
            data_[index] = std::move(data_.back());
            data_.pop_back();
            trickleDown(index, minimum);
        }
        return result;
    }
public:
    // 空の両端優先度キューを構築する。計算量 O(1)。
    DoubleEndedPriorityQueue() = default;
    // 列からmin-max heapを構築する。時間・領域 O(N)。
    explicit DoubleEndedPriorityQueue(std::vector<T> values) : data_(std::move(values)) {
        for (Int i = size() / 2; i-- > 0;)
            trickleDown(i, minLevel(i));
    }
    // 要素数を返す。計算量 O(1)。
    Int size() const {
        return data_.size();
    }
    // 空かを返す。計算量 O(1)。
    bool isEmpty() const {
        return data_.empty();
    }
    // 要素を追加する。償却 O(log(N+1))。
    void push(T item) {
        data_.push_back(std::move(item));
        Int pos = size() - 1;
        if (pos == 0)
            return;
        Int parent = (pos - 1) / 2;
        bool minimum = minLevel(pos);
        if (less(data_[parent], data_[pos], minimum)) {
            std::swap(data_[parent], data_[pos]);
            bubbleUp(parent, !minimum);
        } else {
            bubbleUp(pos, minimum);
        }
    }
    // 空でないキューの最小値を参照する。計算量 O(1)。
    const T &min() const {
        assert(!isEmpty());
        return data_[0];
    }
    // 空でないキューの最大値を参照する。計算量 O(1)。
    const T &max() const {
        assert(!isEmpty());
        return data_[maxIndex()];
    }
    // 最小値を取り除いて返す。計算量 O(log N)。
    T popMin() {
        assert(!isEmpty());
        return popAt(0, true);
    }
    // 最大値を取り除いて返す。計算量 O(log N)。
    T popMax() {
        assert(!isEmpty());
        return popAt(maxIndex(), false);
    }
    // 全要素を破棄する。計算量 O(N)。
    void clear() {
        data_.clear();
    }
};
// 空の両端優先度キューを構築する。計算量 O(1)。
template <class T> auto initDoubleEndedPriorityQueue() {
    return DoubleEndedPriorityQueue<T>();
}
// 列から両端優先度キューを構築する。時間・領域 O(N)。
template <class T> auto toDoubleEndedPriorityQueue(std::vector<T> values) {
    return DoubleEndedPriorityQueue<T>(std::move(values));
}
}
