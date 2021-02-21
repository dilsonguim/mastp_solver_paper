#ifndef H_TIMER_
#define H_TIMER_

#include <chrono>

using namespace std;

class Timer {
public:
   Timer() : time_(0.0) {}

   void Start();

   void Pause();

   void Reset();

   double Read();

private:
  
   chrono::steady_clock::time_point start_time_;
   double time_;
};

#endif
