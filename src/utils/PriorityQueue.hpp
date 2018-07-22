#ifndef METRONOME_PRIORITY_QUEUE_HPP
#define METRONOME_PRIORITY_QUEUE_HPP

#include <boost/assert.hpp>
#include <functional>
#include <string>
#include <vector>

/**
 * @author Bence Cserna 
 */

namespace metronome {

template <typename T>
class PriorityQueue {
public:
    PriorityQueue(unsigned int capacity, std::function<int(const T&, const T&)> comparator)
            : capacity{capacity}, comparator{comparator}, queue{capacity}, size{0} {}

    void push(T& item) {
        if (size == capacity) {
            throw std::overflow_error("Priority queue reached its maximum capacity: " + std::to_string(capacity));
        }

        if (size == 0) {
            queue[0] = &item;
            item.index = 0;
        } else {
            siftUp(size, item);
        }

        ++size;
    }

    T* pop() {
        if (size == 0) {
            return nullptr;
        }

        --size;
        auto result = queue[0];
        auto& x = *queue[size];
        queue[size] = nullptr;

        if (size != 0) {
            siftDown(0, x);
        }

        BOOST_ASSERT_MSG(result->index == 0, "Internal index of top item was not null");
        result->index = std::numeric_limits<unsigned int>::max();
        return const_cast<T*>(result);
    }

    T* top() const {
        if (size == 0) {
            return nullptr;
        }

        return queue[0];
    }

    void clear() { size = 0; }

    void insertOrUpdate(T& item) {
        if (item.index == std::numeric_limits<unsigned int>::max()) {
            // Item is not in the queue yet
            push(item);
        } else {
            // Already in the queue
            update(item);
        }
    }

    void update(T& item) {
        auto index = item.index;

        assert(index != std::numeric_limits<unsigned int>::max() && "Cannot update a node that is not in the queue!");

        siftUp(item.index, item);

        if (item.index == index) {
            siftDown(item.index, item);
        }
    }

    void reorder(std::function<int(const T&, const T&)> comparator) {
        this->comparator = comparator;
        for (unsigned int i = 0; i < size; i++) {
            T& item = *queue[i];
            siftDown(i, item);
        }
    }

    void forEach(std::function<void(T*)> action) {
        for (unsigned int i = 0; i < size; i++) {
            action(queue[i]);
        }
    }

    unsigned int getCapacity() const { return capacity; }

    unsigned int getSize() const { return size; }

    bool isEmpty() const { return size == 0; }

    bool isNotEmpty() const { return size != 0; }

private:
    void siftUp(const unsigned int index, T& item) {
        unsigned int currentIndex = index;
        while (currentIndex > 0) {
            const auto parentIndex = (currentIndex - 1) >> 1;
            const auto parentItem = queue[parentIndex];

            if (comparator(item, *parentItem) >= 0) {
                break;
            }

            // Move parent down and update its index
            queue[currentIndex] = parentItem;
            parentItem->index = currentIndex;
            currentIndex = parentIndex;
        }

        queue[currentIndex] = &item;
        item.index = currentIndex;
    }

    void siftDown(const unsigned int index, T& item) {
        auto currentIndex = index;
        const unsigned int half = size >> 1;

        while (currentIndex < half) {
            auto childIndex = (currentIndex << 1) + 1;
            auto childItem = queue[childIndex];
            auto rightIndex = childIndex + 1;

            if (rightIndex < size && comparator(*childItem, *queue[rightIndex]) > 0) {
                childIndex = rightIndex;
                childItem = queue[rightIndex];
            }

            if (comparator(item, *childItem) <= 0) {
                break;
            }

            queue[currentIndex] = childItem;
            childItem->index = currentIndex;
            currentIndex = childIndex;
        }

        queue[currentIndex] = &item;
        item.index = currentIndex;
    }

    const unsigned int capacity;
    std::function<int(const T&, const T&)> comparator;
    std::vector<T*> queue;
    unsigned int size;
};
}

#endif // METRONOME_PRIORITY_QUEUE_HPP
