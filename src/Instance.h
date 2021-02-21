#ifndef H_INSTANCE_
#define H_INSTANCE_

#include <utility>
#include <vector>

using namespace std;

struct Instance {
  int n_;
  std::vector<pair<int, int>> points_;
  vector<pair<int, int>> edges_;
  vector<int> edge_to_circle_;
  vector<vector<int>> circle_to_edges_;

  struct FaceAdjacency {
    int to_;
    int circle_;
    bool grows_;
  };

  int num_faces_;
  vector<vector<FaceAdjacency>> face_graph_;
  vector<int> inclusion_set_size_;
  vector<double> face_area_;

  vector<vector<int>> hitting_set_;

};

#endif
