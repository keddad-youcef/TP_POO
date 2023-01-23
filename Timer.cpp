//TD5//
#include <iostream>
#include <chrono>
#include "Timer.h"

// Start counting...
T Timer::start() {
    T start = std::chrono::system_clock::now();
    return start;
}


// Stop counting..
T Timer::stop() {
    T end = std::chrono::system_clock::now();
    return end;
}


//Print the difference
void Timer::print(T end, T start) {
    std::chrono::duration<double> elapsed_seconds = end - start;
    std::cout << "elapsed time: " << elapsed_seconds.count() << "s\n";
}