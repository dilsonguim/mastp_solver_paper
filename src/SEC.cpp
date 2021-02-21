#include "SEC.h"
#include <cmath>
#include <iostream>
#include <limits>
#include <queue>
#include <vector>

using namespace std;

constexpr int SCALE = 10000000;

struct Edge {
  int v, rev;
  int cap;
  Edge(int v_, int cap_, int rev_) : v(v_), rev(rev_), cap(cap_) {}
};

struct MaxFlow {
  vector<vector<Edge>> g;
  vector<int> level;
  queue<int> q;
  int flow, n;

  MaxFlow(int n_) : g(n_), level(n_), n(n_) {}

  void addEdge(int u, int v, int cap) {
    if (u == v or cap == 0)
      return;
    Edge e(v, cap, int(g[v].size()));
    Edge r(u, 0, int(g[u].size()));
    g[u].push_back(e);
    g[v].push_back(r);
  }

  bool buildLevelGraph(int src, int sink) {
    fill(level.begin(), level.end(), -1);
    while (not q.empty())
      q.pop();
    level[src] = 0;
    q.push(src);
    while (not q.empty()) {
      int u = q.front();
      q.pop();
      for (auto e = g[u].begin(); e != g[u].end(); ++e) {
        if (not e->cap or level[e->v] != -1)
          continue;
        level[e->v] = level[u] + 1;
        if (e->v == sink)
          return true;
        q.push(e->v);
      }
    }
    return false;
  }

  int blockingFlow(int u, int sink, int f) {
    if (u == sink or not f)
      return f;
    int fu = f;
    for (auto e = g[u].begin(); e != g[u].end(); ++e) {
      if (not e->cap or level[e->v] != level[u] + 1)
        continue;
      int mincap = blockingFlow(e->v, sink, min(fu, e->cap));
      if (mincap) {
        g[e->v][e->rev].cap += mincap;
        e->cap -= mincap;
        fu -= mincap;
      }
    }
    if (f == fu)
      level[u] = -1;
    return f - fu;
  }

  int maxFlow(int src, int sink) {
    flow = 0;
    while (buildLevelGraph(src, sink))
      flow += blockingFlow(src, sink, numeric_limits<int>::max());
    return flow;
  }
};

int RoundCapacity(const double cap) { return cap * SCALE; }

vector<pair<double, vector<int>>> SeparateSEC(const Instance &instance,
                                              const vector<double> &x) {
  vector<pair<double, vector<int>>> cuts;

  vector<double> b(instance.n_);
  for (int i = 0; i < instance.edges_.size(); i++) {
    b[instance.edges_[i].first] += x[i];
    b[instance.edges_[i].second] += x[i];
  }

  const double INF = 2 * instance.n_;

  for (int k = 0; k < instance.n_; k++) {
    MaxFlow g(instance.n_ + 2);
    const int src = instance.n_;
    const int sink = instance.n_ + 1;

    for (int i = 0; i < instance.n_; i++) {
      const double src_cap = i == k ? INF : max(0.5 * b[i] - 1., 0.);
      g.addEdge(src, i, RoundCapacity(src_cap));

      const double sink_cap = i < k ? INF : max(1. - 0.5 * b[i], 0.);
      g.addEdge(i, sink, RoundCapacity(sink_cap));
    }

    for (int i = 0; i < instance.edges_.size(); i++) {
      const int a = instance.edges_[i].first;
      const int b = instance.edges_[i].second;
      const int cap = RoundCapacity(0.5 * x[i]);
      g.addEdge(a, b, cap);
      g.addEdge(b, a, cap);
    }

    const int flow = g.maxFlow(src, sink);

    vector<int> cut;
    for (int i = 0; i < instance.n_; i++) {
      if (g.level[i] != -1) {
        cut.push_back(i);
      }
    }

    double violation = 1 - int(cut.size());
    for (int i = 0; i < instance.edges_.size(); i++) {
      const int a = instance.edges_[i].first;
      const int b = instance.edges_[i].second;
      if (not binary_search(cut.begin(), cut.end(), a))
        continue;
      if (not binary_search(cut.begin(), cut.end(), b))
        continue;
      violation += x[i];
    }

    if (violation > 1.e-4) {
      cuts.push_back({violation, cut});
    }
  }

  return cuts;
}
