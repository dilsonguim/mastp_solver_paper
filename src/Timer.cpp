#include "Timer.h"

void Timer::Start() { start_time_ = chrono::steady_clock::now(); }

void Timer::Pause() {
  auto now = chrono::steady_clock::now();
  time_ += chrono::duration<double, milli>(now - start_time_).count();
}

void Timer::Reset() { time_ = 0.0; }

double Timer::Read() { return time_; }
