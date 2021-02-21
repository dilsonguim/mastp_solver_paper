#include "Benders.h"

#include <cassert>
#include <iostream>
#include <set>

using namespace std;

void DFS(const Instance &instance, const vector<double> &x,
         set<pair<double, int>> &active, vector<bool> &visited, int face,
         vector<double> &w) {

  visited[face] = true;

  if (not active.empty()) {
    int edge = active.begin()->second;
    w[edge] += instance.face_area_[face];
  } else {
    assert(face == 0);
  }

  for (const auto &adj : instance.face_graph_[face]) {
    if (visited[adj.to_])
      continue;

    for (int edge : instance.circle_to_edges_[adj.circle_]) {
      pair<double, int> item(-x[edge], edge);
      if (adj.grows_) {
        active.insert(item);
      } else {
        active.erase(item);
      }
    }

    DFS(instance, x, active, visited, adj.to_, w);

    for (int edge : instance.circle_to_edges_[adj.circle_]) {
      pair<double, int> item(-x[edge], edge);
      if (adj.grows_) {
        active.erase(item);
      } else {
        active.insert(item);
      }
    }
  }
}

vector<double> SeparateBendersCut(const Instance &instance,
                                  const vector<double> &x, const double ub) {

  set<pair<double, int>> active;
  vector<double> w(instance.edges_.size());
  vector<bool> visited(instance.num_faces_);

  DFS(instance, x, active, visited, 0, w);

  double sum = 0.0;

  double ub_of_lb = 0.;
  for (int e = 0; e < instance.edges_.size(); e++) {
    ub_of_lb += w[e] * x[e];
  }
  if (ub_of_lb < ub) {
    return vector<double>();
  }
  double violation = -x[instance.edges_.size()];
  double obj = -violation;
  for (int i = 0; i < instance.edges_.size(); i++) {
    violation += x[i] * w[i];
    sum += w[i];
    ;
  }

  if (obj > 1.e-2)
    violation /= obj;

  if (violation < 1.e-4)
    return vector<double>();

  return w;
}
