#pragma once

namespace cplib {
// 連結ノードの所有権は呼出し側が保持する。挿入はポインタの付替えのみでO(1)。
template <class T> struct DoublyLinkedNode {
    T value{};
    DoublyLinkedNode *prev = nullptr;
    DoublyLinkedNode *next = nullptr;
};

template <class T> struct DoublyLinkedList {
    DoublyLinkedNode<T> *head = nullptr;
    DoublyLinkedNode<T> *tail = nullptr;
};

template <class T>
void insertPrev(DoublyLinkedList<T> &l, DoublyLinkedNode<T> *a, DoublyLinkedNode<T> *b) {
    if (l.head == a) {
        l.head = b;
        a->prev = b;
        b->next = a;
    } else {
        a->prev->next = b;
        b->prev = a->prev;
        a->prev = b;
        b->next = a;
    }
}

template <class T>
void insert(DoublyLinkedList<T> &l, DoublyLinkedNode<T> *a, DoublyLinkedNode<T> *b) {
    if (l.tail == a) {
        l.tail = b;
        a->next = b;
        b->prev = a;
    } else {
        a->next->prev = b;
        b->next = a->next;
        a->next = b;
        b->prev = a;
    }
}
}
