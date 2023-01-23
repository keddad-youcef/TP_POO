//TD5//

#include <chrono>

#ifndef TIMER_H_
#define TIMER_H_

using T =  std::chrono::time_point<std::chrono::system_clock>;

// Timer Class

struct Timer {
    T start();
    T stop();
    void print(T end, T start);
};

#endif