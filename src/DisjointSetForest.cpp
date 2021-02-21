#include "DisjointSetForest.h"

void DisjointSetForest::Reset(const int n) {
  parent.resize(n);
  for (int i = 0; i < n; i++)
    parent[i] = i;
}

int DisjointSetForest::Find(const int x) {
  if (parent[x] != x)
    parent[x] = Find(parent[x]);
  return parent[x];
}

void DisjointSetForest::Join(const int x, const int y) {
  parent[Find(x)] = Find(y);
}
