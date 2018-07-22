#ifndef METRONOME_TIMEMEASUREMENT_HPP
#define METRONOME_TIMEMEASUREMENT_HPP

#include <chrono>
#include <functional>

namespace metronome {

long long int measureNanoTime(std::function<void()> action) {
    const auto begin = std::chrono::high_resolution_clock::now();

    action();

    const auto end = std::chrono::high_resolution_clock::now();

    return std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin).count();
}

long long int currentNanoTime() {
    const auto now = std::chrono::high_resolution_clock::now()
        .time_since_epoch();

    return std::chrono::duration_cast<std::chrono::nanoseconds>(now).count();
}

void logTime(const std::string message = "") {
    static long long int startTime = currentNanoTime();
    LOG(INFO) << (currentNanoTime() - startTime) / 1000000 << message;
}

class StopWatch {
    long long int reset() {
        auto end = std::chrono::high_resolution_clock::now();

        const auto elapsedTime =
            std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin).count();
        
        begin = end;
        return elapsedTime;
    }
    
    long long int elapsed() {
        const auto end = std::chrono::high_resolution_clock::now();
        return std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin).count();
    }

    std::chrono::high_resolution_clock::time_point begin = 
        std::chrono::high_resolution_clock::now();
};

}

#endif // METRONOME_TIMEMEASUREMENT_HPP
