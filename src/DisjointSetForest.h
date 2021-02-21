#ifndef H_DISJOINT_SET_FOREST_
#define H_DISJOINT_SET_FOREST_

#include <vector>

using namespace std;

struct DisjointSetForest {
   vector<int> parent;

   void Reset(const int n);

   int Find(const int x);

   void Join(const int x, const int y);
};



#endif
