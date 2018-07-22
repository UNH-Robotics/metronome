#ifndef METRONOME_LINEARMEMORYPOOL_HPP
#define METRONOME_LINEARMEMORYPOOL_HPP

#include <cstdio>
#include <memory>
namespace metronome {

template <typename T, std::size_t N>
class StaticVector {
    using Storage = typename std::aligned_storage<sizeof(T), alignof(T)>::type;

    Storage* storage = new Storage[N];
    std::size_t size = 0;

public:
    template <typename... Args>
    T* construct(Args&&... args) {
        if (size >= N) // possible error handling
            throw std::overflow_error("ObjectPool reached its maximum capacity: " + std::to_string(N));

        auto* t = new (storage + size) T(std::forward<Args>(args)...);
        ++size;

        return t;
    }

    // Access an object in aligned storage
    const T& operator[](const std::size_t pos) const { return 
            *reinterpret_cast<const T*>(storage + pos); }

    // Delete objects from aligned storage
    ~StaticVector() {
        for (std::size_t pos = 0; pos < size; ++pos) {
            reinterpret_cast<const T*>(storage + pos)->~T();
        }
        
        delete[] storage;
    }
};
}

#endif // METRONOME_LINEARMEMORYPOOL_HPP
