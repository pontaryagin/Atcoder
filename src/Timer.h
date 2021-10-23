#include "MyHeader.h"
#include <chrono>

class Timer {
    chrono::time_point<chrono::system_clock> start;
public:
    Timer() : start(chrono::system_clock::now()){}
    double seconds() { return chrono::duration<double>(chrono::system_clock::now() - start).count(); }
    void reset() {
        start = chrono::system_clock::now();
    }
};
