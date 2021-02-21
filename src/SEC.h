#ifndef H_SEC_
#define H_SEC_

#include <vector>
#include <algorithm>
#include "Instance.h"

using namespace std;

vector<pair<double, vector<int>>> SeparateSEC(const Instance& instance,
const vector<double>& x);

#endif
