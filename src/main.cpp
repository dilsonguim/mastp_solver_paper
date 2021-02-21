#include "Arrangement.h"
#include "Instance.h"
#include "Solver.h"
#include "Timer.h"
#include <algorithm>
#include <cassert>
#include <iostream>
#include <map>
#include <set>
#include <vector>

using namespace std;

int main(int argc, char *argv[]) {
  Timer timer;
  timer.Start();

  Instance instance;
  cin >> instance.n_;
  for (int i = 0; i < instance.n_; i++) {
    int x, y;
    cin >> x >> y;
    instance.points_.push_back({x, y});
  }

  for (int i = 0; i < instance.n_; i++) {
    for (int j = i + 1; j < instance.n_; j++) {
      instance.edges_.push_back({i, j});
    }
  }

  BuildArrangementData(instance);

  timer.Pause();

  Solver solver(instance);
  solver.stats_.build_arrangement_time_ = timer.Read() / 1000.;
  solver.stats_.num_faces_ = instance.num_faces_;
  solver.BuildModel();
  solver.Run();
}
