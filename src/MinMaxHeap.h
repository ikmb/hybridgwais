/*
 *    Copyright (C) 2018-2023 by Lars Wienbrandt and Jan Christian Kässens,
 *    Institute of Clinical Molecular Biology, Kiel University
 *    
 *    This file is part of HybridGWAIS.
 *
 *    HybridGWAIS is free software: you can redistribute it and/or modify
 *    it under the terms of the GNU General Public License as published by
 *    the Free Software Foundation, either version 3 of the License, or
 *    (at your option) any later version.
 *
 *    HybridGWAIS is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *    GNU General Public License for more details.
 *
 *    You should have received a copy of the GNU General Public License
 *    along with HybridGWAIS. If not, see <https://www.gnu.org/licenses/>.
 */

#ifndef MINMAXHEAP_H
#define MINMAXHEAP_H

#include <functional>
#include <vector>

using namespace std;

//template<class T, class Compare = less<T>, class Eq = equal_to<T>>
template<class T, class Compare = less<T>>
class MinMaxHeap {
    vector<T> heap;
    size_t capacity;
    Compare cmp;
//    Eq eq;

    size_t parent(size_t index) const { return (index-1) / 2; }
    size_t left(size_t index) const { return 2*index + 1; }
    size_t right(size_t index) const { return 2*index + 2; }

    uint64_t log2 (uint64_t x) const {
        uint64_t result;
        if(x == 0) return 0;

        asm ("\tbsr %1, %0\n"
        : "=r"(result)
        : "r" (x));
        return result;
    }

    bool isMinLevel (size_t index) const { return (log2(index + 1) % 2) == 0; }
    bool isMaxLevel (size_t index) const { return !isMinLevel(index); }

    template<bool IsMaxLevel>
    void bubbleUpHelper(size_t index) {
        if(index == 0) return;

        size_t theparent = parent(index);
        if(theparent == 0) return;
        size_t grandparent = parent(theparent);

        if(cmp(heap[index], heap[grandparent]) ^ IsMaxLevel) {
            swap(heap[grandparent], heap[index]);
            bubbleUpHelper<IsMaxLevel>(grandparent);
        }
//        else if(eq(heap[index], heap[grandparent])) { // remove element if it is equal according to Eq function
//            deleteElement(index);
//        }
    }

    void bubbleUp(size_t index) {
        if(index == 0) return;

        size_t theparent = parent(index);
        if(isMinLevel(index))
            if(cmp(heap[theparent], heap[index])) {
                swap(heap[theparent], heap[index]);
                bubbleUpHelper<true>(theparent);
            } else {
                bubbleUpHelper<false>(index);
            }
        else
            if(cmp(heap[index], heap[theparent])) {
                swap(heap[theparent], heap[index]);
                bubbleUpHelper<false>(theparent);
            } else {
                bubbleUpHelper<true>(index);
            }
    }

    template<bool IsMaxLevel>
    void trickleDownHelper(size_t index) {
        // smallest on min levels, largest on max levels
        size_t extreme = index;
        size_t leftChild = left(index);


        if(leftChild < heap.size() && (cmp(heap[leftChild], heap[extreme])^IsMaxLevel))
            extreme = leftChild;
        if(leftChild+1 < heap.size() && (cmp(heap[leftChild+1], heap[extreme])^IsMaxLevel))
            extreme = leftChild+1;

        size_t leftGrandchild = left(leftChild);
        for(size_t i = 0; i < 4 && leftGrandchild + i < heap.size(); ++i)
            if(cmp(heap[leftGrandchild + i], heap[extreme]) ^ IsMaxLevel)
                extreme = leftGrandchild + i;

        if(index == extreme)
            return;

        swap(heap[index], heap[extreme]);

        if(extreme - leftChild > 1) {
            if(cmp(heap[parent(extreme)], heap[extreme])^IsMaxLevel)
                swap(heap[parent(extreme)], heap[extreme]);
            trickleDownHelper<IsMaxLevel>(extreme);
        }
    }

    void trickleDown(size_t index) {
        if(isMinLevel(index))
            trickleDownHelper<false>(index);
        else
            trickleDownHelper<true>(index);
    }

    size_t getMaxIndex() const {
        switch(heap.size()) {
        case 0: throw out_of_range("Heap is empty");
        case 1: return 0;
        case 2: return 1;
        default: return cmp(heap[1], heap[2]) ? 2 : 1;
        }
    }

    void deleteElement(size_t index) {

        if(index == heap.size() - 1) {
            heap.pop_back();
            return;
        }

        swap(heap[index], heap[heap.size() - 1]);
        heap.pop_back();
        trickleDown(index);
    }

public:
    using value_type = T;
    using reference = T&;
    using value_compare = Compare;

    explicit MinMaxHeap(size_t capacity_ = 0)
        : capacity(capacity_)
    {
        if(capacity > 0)
            heap.reserve(capacity);
    }

    T popMax() {
        if(heap.empty())
            throw out_of_range("Heap is empty");

        size_t index = getMaxIndex();
        T elem = move(heap[index]);
        deleteElement(index);
        return elem;
    }

    T popMin() {
        if(heap.empty())
            throw out_of_range("Heap is empty");

        T elem = move(heap[0]);
        deleteElement(0);
        return elem;
    }

    const T& getMax() const {
        if(heap.empty())
            throw out_of_range("Heap is empty");
        return heap[getMaxIndex()];
    }

    const T& getMin() const {
        if(heap.empty())
            throw out_of_range("Heap is empty");
        return heap[0];
    }

    bool isEmpty() const {
        return heap.empty();
    }

    size_t getSize() const {
        return heap.size();
    }

    void push(const T& elem) {
        heap.push_back(elem);
        bubbleUp(heap.size() - 1);
    }

    template<typename... Ts>
    void emplace(Ts&&... args) {
        heap.emplace_back(forward<Ts>(args)...);
        bubbleUp(heap.size() - 1);
    }

    const T * data() const {
        return heap.data();
    }
};

#endif /* MIMMAXHEAP_H */
