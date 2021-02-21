#include "Arrangement.h"
#include <cmath>
#include <iostream>
#include <map>

#include <CGAL/Arr_circle_segment_traits_2.h>
#include <CGAL/Arr_consolidated_curve_data_traits_2.h>
#include <CGAL/Arr_extended_dcel.h>
#include <CGAL/Arr_observer.h>
#include <CGAL/Arrangement_2.h>
#include <CGAL/Cartesian.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel_with_sqrt.h>
#include <CGAL/Exact_rational.h>
#include <CGAL/basic.h>

using NumberType = CGAL::Exact_rational;
using Kernel = CGAL::Cartesian<NumberType>;
using Point = Kernel::Point_2;
using Circle = Kernel::Circle_2;

using ArrTraits = CGAL::Arr_circle_segment_traits_2<Kernel>;
using DataTraits = CGAL::Arr_consolidated_curve_data_traits_2<ArrTraits, int>;
using Curve = DataTraits::Curve_2;
using DCEL = CGAL::Arr_face_extended_dcel<DataTraits, int>;
using Arrangement = CGAL::Arrangement_2<DataTraits, DCEL>;

class FaceIndexObserver : public CGAL::Arr_observer<Arrangement> {
private:
  int number_of_faces_;

public:
  FaceIndexObserver(Arrangement &arr)
      : CGAL::Arr_observer<Arrangement>(arr), number_of_faces_(0) {
    CGAL_precondition(arr.is_empty());
    arr.unbounded_face()->set_data(0);
    number_of_faces_++;
  }
  virtual void after_split_face(Face_handle, Face_handle new_face, bool) {
    new_face->set_data(number_of_faces_++);
  }
};

using namespace std;

void EdgeSenseDFS(Instance &instance, int face, set<int> &inclusion,
                  vector<bool> &vis) {
  instance.inclusion_set_size_[face] = inclusion.size();
  vis[face] = true;

  for (int circle : inclusion) {
    for (int edge : instance.circle_to_edges_[circle]) {
      instance.hitting_set_[edge].push_back(face);
    }
  }
  for (auto &adj : instance.face_graph_[face]) {
    adj.grows_ = not inclusion.count(adj.circle_);
    if (vis[adj.to_])
      continue;

    if (adj.grows_) {
      inclusion.insert(adj.circle_);
      EdgeSenseDFS(instance, adj.to_, inclusion, vis);
      inclusion.erase(adj.circle_);
    } else {
      inclusion.erase(adj.circle_);
      EdgeSenseDFS(instance, adj.to_, inclusion, vis);
      inclusion.insert(adj.circle_);
    }
  }
}

pair<double, double> PointDifference(const pair<double, double> &a,
                                     const pair<double, double> &b) {
  return {a.first - b.first, a.second - b.second};
}

double Dot(const pair<double, double> &a, const pair<double, double> &b) {
  return a.first * b.first + a.second * b.second;
}

double FaceArea(Instance &instance,
                Arrangement::Ccb_halfedge_const_circulator ccb,
                const vector<Circle> &circles) {
  double area = 0.;
  auto curr = ccb;
  Arrangement::Halfedge_const_handle he;
  do {
    he = curr;

    const auto source = he->source()->point();
    const auto target = he->target()->point();
    const int circle = *he->curve().data().begin();
    const auto squared_radius = circles[circle].squared_radius();
    const auto center = circles[circle].center();

    pair<double, double> approximate_source(CGAL::to_double(source.x()),
                                            CGAL::to_double(source.y()));
    pair<double, double> approximate_target(CGAL::to_double(target.x()),
                                            CGAL::to_double(target.y()));
    pair<double, double> approximate_center(CGAL::to_double(center.x()),
                                            CGAL::to_double(center.y()));
    double r2 = CGAL::to_double(squared_radius);

    auto p = PointDifference(approximate_source, approximate_center);
    auto q = PointDifference(approximate_target, approximate_center);
    double costheta = Dot(p, q) / r2;
    double theta = acos(costheta);

    int face = he->face()->data();
    int ngb_face = he->twin()->face()->data();
    if (instance.inclusion_set_size_[face] <
        instance.inclusion_set_size_[ngb_face]) {
      theta *= -1.;
    }

    auto qp = PointDifference(q, p);
    double delta = approximate_center.first * qp.second -
                   approximate_center.second * qp.first;
    double line_integral = (theta * r2 + delta) / 2.;
    area += line_integral;

    curr++;
  } while (curr != ccb);
  return area;
}

void AddFaceEdges(Instance &instance, int from,
                  Arrangement::Ccb_halfedge_const_circulator ccb) {
  set<int> neighbours;
  auto curr = ccb;
  Arrangement::Halfedge_const_handle he;
  do {
    he = curr;
    // int circ = he->curve().data().begin();

    int to = he->twin()->face()->data();
    if (not neighbours.count(to)) {
      neighbours.insert(to);
      int circle = *he->curve().data().begin();
      Instance::FaceAdjacency adj;
      adj.to_ = to;
      adj.circle_ = circle;
      adj.grows_ = false;
      instance.face_graph_[from].push_back(adj);
    }
    curr++;
  } while (curr != ccb);
}

void BuildArrangementData(Instance &instance) {

  vector<Curve> curves;
  vector<Circle> circles;
  for (int i = 0; i < instance.edges_.size(); i++) {
    const auto &edge = instance.edges_[i];
    Point a(instance.points_[edge.first].first,
            instance.points_[edge.first].second);
    Point b(instance.points_[edge.second].first,
            instance.points_[edge.second].second);
    Circle c(a, b);

    int circle_index = 0;
    while (circle_index < circles.size() and circles[circle_index] != c) {
      circle_index++;
    }

    instance.edge_to_circle_.push_back(circle_index);
    if (circle_index < circles.size()) {
      instance.circle_to_edges_[circle_index].push_back(i);
    } else {
      circles.push_back(c);
      curves.emplace_back(Curve(c, circle_index));
      instance.circle_to_edges_.push_back({i});
    }
  }

  Arrangement arrangement;
  FaceIndexObserver face_indexer(arrangement);
  insert(arrangement, curves.begin(), curves.end());

  instance.face_graph_.resize(arrangement.number_of_faces());
  for (auto fit = arrangement.faces_begin(); fit != arrangement.faces_end();
       ++fit) {
    if (fit->is_unbounded()) {
      for (auto hole = fit->holes_begin(); hole != fit->holes_end(); ++hole) {
        AddFaceEdges(instance, fit->data(), *hole);
      }
    } else {
      AddFaceEdges(instance, fit->data(), fit->outer_ccb());
    }
  }
  set<int> inclusion;
  vector<bool> vis(arrangement.number_of_faces());
  instance.inclusion_set_size_.resize(arrangement.number_of_faces());
  instance.hitting_set_.resize(instance.edges_.size());
  EdgeSenseDFS(instance, 0, inclusion, vis);

  instance.face_area_.resize(arrangement.number_of_faces());
  for (auto fit = arrangement.faces_begin(); fit != arrangement.faces_end();
       ++fit) {
    if (not fit->is_unbounded()) {
      instance.face_area_[fit->data()] =
          max(FaceArea(instance, fit->outer_ccb(), circles), 0.);
    }
  }

  instance.num_faces_ = arrangement.number_of_faces();
}
