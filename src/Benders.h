#ifndef H_BENDERS_
#define H_BENDERS_

#include <vector>
#include <algorithm>
#include "Instance.h"

using namespace std;

vector<double> SeparateBendersCut(const Instance& instance,
const vector<double>& x, const double ub = 0.);

#endif
