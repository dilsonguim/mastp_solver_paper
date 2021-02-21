#include "Solver.h"
#include "DisjointSetForest.h"

#include <endian.h>
#include <iomanip>
#include <iostream>

using namespace std;

void Solver::BuildModel() {
  int open_cplex_status = 0;
  env_ = CPXopenCPLEX(&open_cplex_status);
  if (open_cplex_status) {
    cerr << "CPLEX: Can't create environment." << endl;
    exit(1);
  }

  if (CPXsetlogfilename(env_, "logfile.txt", "w")) {
    cerr << "CPLEX: Can't setup log file!" << endl;
    exit(1);
  }

  int create_prob_status = 0;
  model_ = CPXcreateprob(env_, &create_prob_status, "mastp");
  if (model_ == nullptr) {
    cerr << "CPLEX: Can't create model." << endl;
    exit(1);
  }

  // Variables: [edges, faces]

#ifdef BENDERS
  const int num_vars = instance_.edges_.size() + 1;
  vector<double> var_obj(num_vars);
  vector<char> var_type(num_vars);
  vector<int> binary_variables;
  for (int i = 0; i < instance_.edges_.size(); i++) {
    var_type[i] = CPX_BINARY;
    binary_variables.push_back(i);
  }
  const int z_var_index = instance_.edges_.size();
  var_type[z_var_index] = CPX_CONTINUOUS;
  var_obj[z_var_index] = 1.;

  if (CPXaddcols(env_, model_, num_vars, 0, var_obj.data(), nullptr, nullptr,
                 nullptr, nullptr, nullptr, nullptr)) {
    cerr << "CPLEX: Can't add variables to the model.\n";
    exit(1);
  }

  if (CPXchgctype(env_, model_, binary_variables.size(),
                  binary_variables.data(), var_type.data())) {
    cerr << "CPLEX: Can't change variables to binary!" << endl;
    exit(1);
  }

  // Warm start benders
  vector<double> x(num_vars);
  for (int k = 0; k < num_vars - 1 and false; k++) {
    x[k] = 1.;

    auto cut = SeparateBendersCut(instance_, x);

    vector<int> bc_ind(num_vars);
    vector<double> bc_val(num_vars);
    for (int i = 0; i < instance_.edges_.size(); i++) {
      bc_ind[i] = i;
      bc_val[i] = cut[i];
    }
    bc_ind.back() = instance_.edges_.size();
    bc_val.back() = -1.;

    int beg[1] = {0};
    char sense[1] = {'L'};
    double rhs[1] = {0.};

    if (CPXaddrows(env_, model_, 0, 1, bc_ind.size(), rhs, sense, beg,
                   bc_ind.data(), bc_val.data(), nullptr, nullptr)) {
      cerr << "CPLEX: Can't add bc warm start constraint." << endl;
      exit(1);
    }

    x[k] = 0.;
  }

#else
  const int num_vars = instance_.edges_.size() + instance_.num_faces_;
  vector<double> var_obj(num_vars);
  vector<char> var_type(num_vars);
  vector<int> binary_variables;
  for (int i = 0; i < instance_.edges_.size(); i++) {
    var_type[i] = CPX_BINARY;
    binary_variables.push_back(i);
  }
  for (int i = 0; i < instance_.num_faces_; i++) {
    const int var_index = instance_.edges_.size() + i;
    var_type[var_index] = CPX_CONTINUOUS;
    var_obj[var_index] = instance_.face_area_[i];
  }
  if (CPXaddcols(env_, model_, num_vars, 0, var_obj.data(), nullptr, nullptr,
                 nullptr, nullptr, nullptr, nullptr)) {
    cerr << "CPLEX: Can't add variables to the model.\n";
    exit(1);
  }

  if (CPXchgctype(env_, model_, binary_variables.size(),
                  binary_variables.data(), var_type.data())) {
    cerr << "CPLEX: Can't change variables to binary!" << endl;
    exit(1);
  }

  // Covering constraints
  vector<bool> visited(instance_.num_faces_);
  set<int> inclusion_set;

#ifdef PREPROCESS
  AddReducedCoveringConstraintsDFS(visited, 0, inclusion_set);
#else
  AddCoveringConstraintsDFS(visited, 0, inclusion_set);
#endif
#endif

  // Spanning tree constraints

  AddCardinalityConstraint();
}

int Solver::ArcIndex(const int root, const int edge,
                     const bool reversed) const {
  return instance_.edges_.size() * (2 * root + 1) + 2 * edge + reversed;
}

void Solver::AddImplicationConstraint(const int a, const int b) {
  const int a_var = instance_.edges_.size() + a;
  const int b_var = instance_.edges_.size() + b;
  int ind[2] = {b_var, a_var};
  double val[2] = {1.0, -1.0};
  int beg[1] = {0};
  char sense[1] = {'G'};

  if (CPXaddrows(env_, model_, 0, 1, 2, nullptr, sense, beg, ind, val, nullptr,
                 nullptr)) {
    cerr << "CPLEX: Can't add constraint." << endl;
    exit(1);
  }
}

void Solver::AddCoveringConstraint(const int face, const int edge) {
  const int face_var = instance_.edges_.size() + face;
  int ind[2] = {face_var, edge};
  double val[2] = {1.0, -1.0};
  int beg[1] = {0};
  char sense[1] = {'G'};

  if (CPXaddrows(env_, model_, 0, 1, 2, nullptr, sense, beg, ind, val, nullptr,
                 nullptr)) {
    cerr << "CPLEX: Can't add constraint." << endl;
    exit(1);
  }
}

void Solver::AddArcEdgeCouplingConstraints() {
  const int num_edges = instance_.edges_.size();
  for (int r = 0; r < instance_.n_; r++) {
    for (int e = 0; e < num_edges; e++) {
      int ind[3] = {e, ArcIndex(r, e, false), ArcIndex(r, e, true)};
      double val[3] = {-1., 1., 1.};
      int beg[1] = {0};
      char sense[1] = {'E'};
      double rhs[1] = {0.};
      if (CPXaddrows(env_, model_, 0, 1, 3, rhs, sense, beg, ind, val, nullptr,
                     nullptr)) {
        cerr << "CPLEX: Can't add constraint." << endl;
        exit(1);
      }
    }
  }
}

void Solver::AddCardinalityConstraint() {
  const int num_edges = instance_.edges_.size();
  vector<int> card_ind(num_edges);
  vector<double> card_val(num_edges, 1.0);
  for (int i = 0; i < num_edges; i++) {
    card_ind[i] = i;
  }

  int beg[1] = {0};
  char sense[1] = {'E'};
  double rhs[1] = {double(instance_.n_ - 1)};

  if (CPXaddrows(env_, model_, 0, 1, card_ind.size(), rhs, sense, beg,
                 card_ind.data(), card_val.data(), nullptr, nullptr)) {
    cerr << "CPLEX: Can't add constraint." << endl;
    exit(1);
  }
}

void Solver::AddCoveringConstraintsDFS(vector<bool> &visited, const int face,
                                       set<int> &circles) {
  visited[face] = true;

  for (const int circle : circles) {
    for (const int edge : instance_.circle_to_edges_[circle]) {
      AddCoveringConstraint(face, edge);
    }
  }

  for (const auto &e : instance_.face_graph_[face]) {
    if (visited[e.to_])
      continue;

    if (e.grows_) {
      circles.insert(e.circle_);
    } else {
      circles.erase(e.circle_);
    }

    AddCoveringConstraintsDFS(visited, e.to_, circles);

    if (e.grows_) {
      circles.erase(e.circle_);
    } else {
      circles.insert(e.circle_);
    }
  }
}

void Solver::AddReducedCoveringConstraintsDFS(vector<bool> &visited,
                                              const int face,
                                              set<int> &circles) {
  visited[face] = true;

  int parent = -1;
  int extra_circle = 0;

  for (const auto &e : instance_.face_graph_[face]) {
    if (visited[e.to_])
      continue;

    if (e.grows_) {
      circles.insert(e.circle_);
    } else {
      circles.erase(e.circle_);

      const int new_w = instance_.circle_to_edges_[e.circle_].size();
      const int old_w = instance_.circle_to_edges_[extra_circle].size();

      if (parent == -1 or new_w < old_w) {
        parent = e.to_;
        extra_circle = e.circle_;
      }
    }

    AddReducedCoveringConstraintsDFS(visited, e.to_, circles);

    if (e.grows_) {
      circles.erase(e.circle_);
    } else {
      circles.insert(e.circle_);
    }
  }

  if (parent == -1) {
    for (const int circle : circles) {
      for (const int edge : instance_.circle_to_edges_[circle]) {
        AddCoveringConstraint(face, edge);
      }
    }
  } else {
    AddImplicationConstraint(parent, face);
    for (const int edge : instance_.circle_to_edges_[extra_circle]) {
      AddCoveringConstraint(face, edge);
    }
  }
}

void AddBendersCut(Solver &solver, void *cbdata, int wherefrom,
                   const vector<double> &opt_cut) {
  const auto &instance = solver.instance_;

  vector<int> xs = {3, 5, 8, 18, 23, 28, 31, 47, 50, 55, 56};
  double norm = 0.;
  for (const double &e : opt_cut)
    norm += e * e;
  norm += 1;
  norm = sqrt(norm);

  vector<int> ind;
  vector<double> val;
  double xseval = 0.;
  for (int i = 0; i < instance.edges_.size(); i++) {
    if (binary_search(xs.begin(), xs.end(), i))
      xseval += opt_cut[i];
    if (fabs(opt_cut[i]) < 1.e-9)
      continue;

    ind.push_back(i);
    val.push_back(opt_cut[i]);
  }
  ind.push_back(instance.edges_.size());
  val.push_back(-1.);

  if (CPXcutcallbackadd(solver.env_, cbdata, wherefrom, ind.size(), 0., 'L',
                        ind.data(), val.data(), CPX_USECUT_PURGE)) {

    cerr << "CPLX: Can't add Benders cut" << endl;
    exit(1);
  }
}

void AddCutset(Solver &solver, void *cbdata, int wherefrom, const int root,
               const vector<int> &arcs) {

  const Instance &instance = solver.instance_;

  vector<int> ind = arcs;
  for (auto &i : ind)
    i += (root * 2 + 1) * instance.edges_.size();
  vector<double> val(arcs.size(), 1.);

  if (CPXcutcallbackadd(solver.env_, cbdata, wherefrom, ind.size(), 1., 'G',
                        ind.data(), val.data(), CPX_USECUT_PURGE)) {

    cerr << "CPLX: Can't add Cutset." << endl;
    exit(1);
  }
}

void AddSEC(Solver &solver, void *cbdata, int wherefrom, const vector<int> &S) {
  if (S.empty())
    return;

  const Instance &instance = solver.instance_;

  vector<int> ind;
  vector<double> val;
  vector<int> xs = {3, 5, 8, 18, 23, 28, 31, 47, 50, 55, 56};
  double xs_sec_vio = 0.;
  for (int i = 0; i < instance.edges_.size(); i++) {
    const int a = instance.edges_[i].first;
    const int b = instance.edges_[i].second;
    if (not binary_search(S.begin(), S.end(), a))
      continue;
    if (not binary_search(S.begin(), S.end(), b))
      continue;

    ind.push_back(i);
    val.push_back(1.0);
  }
  const double rhs = S.size() - 1;
  xs_sec_vio -= rhs;
  for (int i = 0; i < ind.size(); i++) {
    if (binary_search(xs.begin(), xs.end(), ind[i])) {
      xs_sec_vio += val[i];
    }
  }

  int err = -1;

  if (CPXcutcallbackadd(solver.env_, cbdata, wherefrom, ind.size(), rhs, 'L',
                        ind.data(), val.data(), CPX_USECUT_PURGE)) {

    cerr << "CPLX: Can't add SEC." << err << endl;
    exit(1);
  }
}

double Solver::DfsComputeObj(vector<bool> &visited, int face, int count,
                             double *x) {
  double obj = 0.;
  visited[face] = true;

  if (count) {
    obj += instance_.face_area_[face];
  }

  for (const auto &e : instance_.face_graph_[face]) {
    if (visited[e.to_])
      continue;

    int delta = 0;
    for (int i : instance_.circle_to_edges_[e.circle_]) {
      if (x[i] > 0.) {
        delta += e.grows_ ? 1 : -1;
        break;
      }
    }
    obj += DfsComputeObj(visited, e.to_, count + delta, x);
  }

  return obj;
}

void Solver::DfsComputeY(vector<bool> &visited, int face, int count,
                         double *x) {
  visited[face] = true;

  if (count) {
    x[instance_.edges_.size() + face] = 1.0;
  }

  for (const auto &e : instance_.face_graph_[face]) {
    if (visited[e.to_])
      continue;

    int delta = 0;
    for (int i : instance_.circle_to_edges_[e.circle_]) {
      if (x[i] > 0.) {
        delta += e.grows_ ? 1 : -1;
        break;
      }
    }
    DfsComputeY(visited, e.to_, count + delta, x);
  }
}

vector<int> Solver::ClosestTree(double *x) {
  vector<int> tree;
  vector<int> ord(instance_.edges_.size());
  for (int i = 0; i < instance_.edges_.size(); i++) {
    ord[i] = i;
  }
  sort(ord.begin(), ord.end(), [x](int i, int j) { return x[i] > x[j]; });

  DisjointSetForest components;
  components.Reset(instance_.n_);
  for (int i : ord) {
    int a = instance_.edges_[i].first;
    int b = instance_.edges_[i].second;

    if (components.Find(a) == components.Find(b))
      continue;

    components.Join(a, b);

    tree.push_back(i);
  }

  return tree;
}

int FindBestFaceInTwoCircles(const Instance &instance, const vector<double> &x,
                             const int c1, const int c2) {
  if (instance.hitting_set_[c1].size() > instance.hitting_set_[c2].size())
    return FindBestFaceInTwoCircles(instance, x, c2, c1);

  double best_value = numeric_limits<double>::infinity();
  int best_face = 0;

  for (const int f : instance.hitting_set_[c1]) {
    if (binary_search(instance.hitting_set_[c2].begin(),
                      instance.hitting_set_[c2].end(), f)) {

      const double value = x[instance.edges_.size() + f];
      if (value < best_value) {
        best_value = value;
        best_face = f;
      }
    }
  }

  return best_face;
}

int SeparateTriangleInequalities(Solver &solver, int wherefrom, void *cbdata,
                                 const vector<double> &x) {
  // The idea here is to fix one side of the triangle and compute
  // the best two faces for each vertex completing the triangle.

  const Instance &instance = solver.instance_;

  vector<pair<int, int>> best_face(instance.n_);
  vector<pair<int, int>> tip_edge(instance.n_);

  int num_vio_triang = 0;
  for (int k = 0; k < instance.edges_.size(); k++) {
    const int c_k = instance.edge_to_circle_[k];

    fill(best_face.begin(), best_face.end(), pair<int, int>(0, 0));

    const int r = instance.edges_[k].first;
    const int s = instance.edges_[k].second;

    for (int i = 0; i < instance.edges_.size(); i++) {
      if (i == k)
        continue;

      const int c_i = instance.edge_to_circle_[i];

      if (instance.edges_[i].first == r or instance.edges_[i].second == r) {
        const int t = instance.edges_[i].first == r ? instance.edges_[i].second
                                                    : instance.edges_[i].first;

        best_face[t].first = FindBestFaceInTwoCircles(instance, x, c_k, c_i);
        tip_edge[t].first = i;
      }

      if (instance.edges_[i].first == s or instance.edges_[i].second == s) {
        const int t = instance.edges_[i].first == s ? instance.edges_[i].second
                                                    : instance.edges_[i].first;

        best_face[t].second = FindBestFaceInTwoCircles(instance, x, c_k, c_i);
        tip_edge[t].second = i;
      }
    }

    for (int t = 0; t < instance.n_; t++) {
      if (best_face[t].first == 0 or best_face[t].second == 0)
        continue;

      const int i = tip_edge[t].first;
      const int j = tip_edge[t].second;
      const int f_a = best_face[t].first;
      const int f_b = best_face[t].second;
      const int var_f_a = instance.edges_.size() + f_a;
      const int var_f_b = instance.edges_.size() + f_b;

      const double violation = x[i] + x[j] + x[k] - x[var_f_a] - x[var_f_b];
      if (violation > 1.e-5) {
        num_vio_triang++;

        vector<int> ind = {i, j, k, var_f_a, var_f_b};
        vector<double> coef = {-1., -1., -1., 1., 1.};

        if (var_f_a == var_f_b) {
          ind.pop_back();
          coef.pop_back();
          coef.back() += 1.0;
          continue;
        }

        int err = 0;

        if (CPXcutcallbackadd(solver.env_, cbdata, wherefrom, ind.size(), 0.,
                              'G', ind.data(), coef.data(), CPX_USECUT_PURGE)) {

          cerr << "CPLX: Can't add Triangle Inequality." << err << endl;
          exit(1);
        }
      }
    }
  }

  return num_vio_triang;
}

void Solver::ImproveTree(vector<int> &tree) {
  vector<int> count(instance_.num_faces_);
  for (int edge : tree) {
    for (int face : instance_.hitting_set_[edge]) {
      count[face]++;
    }
  }

  DisjointSetForest dsf;

  while (true) {
    int best_remove = -1;
    int best_add = -1;
    double best_delta = 0.;
    for (int remove : tree) {
      double area_decrease = 0.;
      for (int face : instance_.hitting_set_[remove]) {
        count[face]--;
        if (count[face] == 0) {
          area_decrease += instance_.face_area_[face];
        }
      }

      dsf.Reset(instance_.n_);
      for (int edge : tree) {
        if (edge != remove) {
          dsf.Join(instance_.edges_[edge].first, instance_.edges_[edge].second);
        }
      }

      for (int add = 0; add < instance_.edges_.size(); add++) {
        if (add == remove)
          continue;
        if (dsf.Find(instance_.edges_[add].first) ==
            dsf.Find(instance_.edges_[add].second))
          continue;

        double area_increase = 0;
        for (int face : instance_.hitting_set_[add]) {
          if (count[face] == 0) {
            area_increase += instance_.face_area_[face];
          }
        }

        double delta = area_increase - area_decrease;
        if (delta < best_delta) {
          best_delta = delta;
          best_remove = remove;
          best_add = add;
        }
      }

      for (int face : instance_.hitting_set_[remove]) {
        count[face]++;
      }
    }

    if (best_delta > -1.e-2)
      break;

    for (int &edge : tree) {
      if (edge == best_remove) {
        edge = best_add;
        break;
      }
    }
    for (int face : instance_.hitting_set_[best_remove]) {
      count[face]--;
    }
    for (int face : instance_.hitting_set_[best_add]) {
      count[face]++;
    }
  }
}

int CPXPUBLIC HeuristicCallback(CPXCENVptr env, void *cbdata, int wherefrom,
                                void *cbhandle, double *objval_p, double *x,
                                int *checkfeas_p, int *useraction_p) {

  auto &solver = *(Solver *)cbhandle;
  const Instance &instance = solver.instance_;

  solver.stats_.num_runs_heuristic_++;

  vector<int> tree = solver.ClosestTree(x);
  solver.ImproveTree(tree);

  *checkfeas_p = 0;

  if (tree.size() != instance.n_ - 1) {
    *useraction_p = CPX_CALLBACK_DEFAULT;
    return 0;
  }

#ifdef BENDERS
  const int num_vars = instance.edges_.size() + 1;
#else
  const int num_vars = instance.edges_.size() + instance.num_faces_;
#endif
  fill(x, x + num_vars, 0);
  for (int i : tree) {
    x[i] = 1.0;
  }

  vector<bool> visited(instance.num_faces_);
#ifdef BENDERS
  *objval_p = solver.DfsComputeObj(visited, 0, 0, x);
  x[instance.edges_.size()] = *objval_p;
#else
  solver.DfsComputeY(visited, 0, 0, x);

  *objval_p = 0.;
  for (int i = 0; i < instance.num_faces_; i++) {
    if (x[i + instance.edges_.size()] == 1.) {
      *objval_p += instance.face_area_[i];
    }
  }
#endif

  *useraction_p = CPX_CALLBACK_SET;
  int node_depth = -1;
  if (CPXgetcallbacknodeinfo(env, cbdata, wherefrom, 0,
                             CPX_CALLBACK_INFO_NODE_DEPTH, &node_depth)) {

    cerr << "CPLEX - CutCallback - Failed to get node depth" << endl;
  }

  if (node_depth == 0) {
    solver.stats_.objective_upper_bound_root_ =
        max(solver.stats_.objective_upper_bound_root_, *objval_p);
  }

  solver.benders_early_stop_ub = min(solver.benders_early_stop_ub, *objval_p);

  return 0;
}

static int CPXPUBLIC CutCallback(CPXCENVptr env, void *cbdata, int wherefrom,
                                 void *cbhandle, int *useraction_p) {

  /*if (wherefrom == 106) {
    *useraction_p = CPX_CALLBACK_DEFAULT;
    return 0;
  }
  */
  int node_depth = -1;

  if (CPXgetcallbacknodeinfo(env, cbdata, wherefrom, 0,
                             CPX_CALLBACK_INFO_NODE_DEPTH, &node_depth)) {

    cerr << "CPLEX - CutCallback - Failed to get node depth" << endl;
  }

  auto &solver = *(Solver *)cbhandle;
  const Instance &instance = solver.instance_;
  auto &stats = solver.stats_;
#ifdef BENDERS
  const int num_vars = instance.edges_.size() + 1;
#else
  const int num_vars = instance.edges_.size() + instance.num_faces_;
#endif
  vector<double> x(num_vars);
  if (CPXgetcallbacknodex(env, cbdata, wherefrom, x.data(), 0, num_vars - 1)) {
    cerr << "CPLEX - CutCallback - Failed to get x!" << endl;
    exit(1);
  }

  if (wherefrom == CPX_CALLBACK_MIP_CUT_FEAS) {
    int sol_src = -1;
    if (CPXgetcallbacknodeinfo(env, cbdata, wherefrom, 0,
                               CPX_CALLBACK_INFO_LAZY_SOURCE, &sol_src)) {

      cerr << "CPLEX - CutCallback - Failed to get candidate solution source."
           << endl;
    }
  }

  int num_cuts_added = 0;

  auto secs = SeparateSEC(instance, x);
  stats.num_sec_separations_++;
  if (node_depth == 0) {
    stats.num_sec_separations_root_++;
  }

  sort(secs.begin(), secs.end());
  reverse(secs.begin(), secs.end());
  if (not secs.empty())
    secs.resize(1);
  for (const auto &sec : secs) {
    AddSEC(solver, cbdata, wherefrom, sec.second);
    num_cuts_added++;

    stats.num_sec_++;
    if (node_depth == 0) {
      stats.num_sec_root_++;
    }
  }

#ifdef SEPARATE_TRIANGLE_CUTS
  if (num_cuts_added == 0) {
    int num_triineq =
        SeparateTriangleInequalities(solver, wherefrom, cbdata, x);
    num_cuts_added++;

    stats.num_triineq_separations_++;
    stats.num_triineq_ += num_triineq;
    if (node_depth == 0) {
      stats.num_triineq_separations_root_++;
      stats.num_triineq_root_ += num_triineq;
    }
  }
#endif

#ifdef BENDERS

  bool integral = true;
  for (int e = 0; e < instance.edges_.size(); e++) {
    if (min(fabs(x[e] - 0), fabs(1 - x[e])) > 1.e-4) {
      integral = false;
      break;
    }
  }
  if (num_cuts_added == 0 or integral) {
    double early_stop_th = integral or solver.benders_early_stop_ub ==
                                           numeric_limits<double>::infinity()
                               ? 0.
                               : solver.benders_early_stop_ub;

    auto benders_cut = SeparateBendersCut(instance, x, early_stop_th);
    solver.stats_.num_benders_separations_++;
    if (node_depth == 0) {
      solver.stats_.num_benders_separations_root_++;
    }

    if (not benders_cut.empty()) {
      AddBendersCut(solver, cbdata, wherefrom, benders_cut);
      num_cuts_added++;
      double bvio = -x.back();
      for (int i = 0; i < instance.edges_.size(); i++) {
        bvio += x[i] * benders_cut[i];
      }

      solver.stats_.num_benders_cuts_++;
      if (node_depth == 0) {
        solver.stats_.num_benders_cuts_root_++;
      }
    }
  }
#endif

  if (num_cuts_added) {
    *useraction_p = CPX_CALLBACK_SET;
  } else {
    *useraction_p = CPX_CALLBACK_DEFAULT;
  }

  if (node_depth == 0 and false) {
    double dual_bound = 0.;
    if (CPXgetcallbacknodeinfo(env, cbdata, wherefrom, 0,
                               CPX_CALLBACK_INFO_NODE_OBJVAL, &dual_bound)) {

      cerr << "CPLEX - CutCallback - Failed to get dual bound" << endl;
    }

    solver.stats_.objective_lower_bound_root_ = dual_bound;
  }

  return 0;
}

void Solver::SetCPLEXParameters() {

  if (CPXsetdblparam(env_, CPXPARAM_TimeLimit, 3 * 60 * 60)) {
    cerr << "CPLEX - Cant't set time limit!" << endl;
    exit(1);
  }

  if (CPXsetintparam(env_, CPXPARAM_Threads, 1)) {
    cerr << "CPLEX - Cant't set the number of threads!" << endl;
    exit(1);
  }
}

void Solver::Run() {
  SetCPLEXParameters();

  stats_.num_sec_ = 0;
  stats_.num_sec_separations_ = 0;
  stats_.num_sec_root_ = 0;
  stats_.num_sec_separations_root_ = 0;
  stats_.num_triineq_ = 0;
  stats_.num_triineq_separations_ = 0;
  stats_.num_triineq_root_ = 0;
  stats_.num_triineq_separations_root_ = 0;
  stats_.num_runs_heuristic_ = 0;
  stats_.num_vars_ = CPXgetnumcols(env_, model_);
  stats_.num_constraints_ = CPXgetnumrows(env_, model_);
  stats_.num_benders_cuts_ = 0;
  stats_.num_benders_cuts_root_ = 0;
  stats_.num_benders_separations_ = 0;
  stats_.num_benders_separations_root_ = 0;

  benders_early_stop_ub = numeric_limits<double>::infinity();

  if (CPXwriteprob(env_, model_, "/tmp/cpx.lp", nullptr)) {
    cerr << "CPLEX: Can't write problem!" << endl;
    exit(1);
  }

  CPXsetintparam(env_, CPXPARAM_ScreenOutput, CPX_ON);

  if (CPXsetlazyconstraintcallbackfunc(env_, CutCallback, this)) {
    cerr << "CPLEX: Failed to setup lazy constraints callback!" << endl;
    exit(1);
  }

  if (CPXsetusercutcallbackfunc(env_, CutCallback, this)) {
    cerr << "CPLEX: Failed to setup user cuts callback!" << endl;
    exit(1);
  }

  if (CPXsetheuristiccallbackfunc(env_, HeuristicCallback, this)) {
    cerr << "CPLEX: Failed to setup  heuristic callback!" << endl;
    exit(1);
  }

  // Get start time
  double start_timestamp = 0.;
  if (CPXgettime(env_, &start_timestamp)) {
    cerr << "CPLEX: Can't get timestamp!" << endl;
    exit(1);
  }

  // Solve root

  if (CPXsetlongparam(env_, CPXPARAM_MIP_Limits_Nodes, 0)) {
    cerr << "CPLEX: Cant't limit the number of nodes!" << endl;
    exit(1);
  }
  if (CPXmipopt(env_, model_)) {
    cerr << "CPLEX: Error during optimization!" << endl;
    exit(1);
  }

  // Get root stats
  double root_end_timestamp = 0.;
  if (CPXgettime(env_, &root_end_timestamp)) {
    cerr << "CPLEX: Can't get timestamp!" << endl;
    exit(1);
  }
  stats_.runtime_root_ = root_end_timestamp - start_timestamp;

  if (CPXgetobjval(env_, model_, &stats_.objective_upper_bound_root_)) {
    cerr << "CPLEX: Can't get upper bound!" << endl;
    exit(1);
  }
  if (CPXgetbestobjval(env_, model_, &stats_.objective_lower_bound_root_)) {
    cerr << "CPLEX: Can't get lower bound!" << endl;
    exit(1);
  }

  vector<double> root_x(instance_.edges_.size() + 1);
  if (CPXgetx(env_, model_, root_x.data(), 0, root_x.size() - 1)) {
    cerr << "CPLEX: Can't get root x!" << endl;
    exit(1);
  }

  cout << endl
       << endl
       << "ROOT ENDED #########################################################"
       << endl
       << endl;

  // Solve the rest of the nodes
  if (CPXsetlongparam(env_, CPXPARAM_MIP_Limits_Nodes, 9223372036800000000ll)) {
    cerr << "CPLEX: Cant't restore the node limit!" << endl;
    exit(1);
  }
  if (CPXmipopt(env_, model_)) {
    cerr << "CPLEX: Error during optimization!" << endl;
    exit(1);
  }

  // Get final stats
  double end_timestamp = 0.;
  if (CPXgettime(env_, &end_timestamp)) {
    cerr << "CPLEX: Can't get timestamp!" << endl;
    exit(1);
  }
  stats_.runtime_ = end_timestamp - start_timestamp;
  if (CPXgetobjval(env_, model_, &stats_.objective_upper_bound_)) {
    cerr << "CPLEX: Can't get upper bound!" << endl;
    exit(1);
  }
  if (CPXgetbestobjval(env_, model_, &stats_.objective_lower_bound_)) {
    cerr << "CPLEX: Can't get lower bound!" << endl;
    exit(1);
  }
  stats_.num_nodes_ = CPXgetnodecnt(env_, model_);
  stats_.num_nodes_left_ = CPXgetnodeleftcnt(env_, model_);

  cout << fixed << setprecision(6);
  cout << "STATS:" << endl;
  cout << "{" << endl;
  cout << "\"objective_lower_bound_root\": "
       << stats_.objective_lower_bound_root_ << "," << endl;
  cout << "\"objective_upper_bound_root\": "
       << stats_.objective_upper_bound_root_ << "," << endl;
  cout << "\"runtime_root\": " << stats_.runtime_root_ << "," << endl;
  cout << "\"runtime\": " << stats_.runtime_ << "," << endl;
  cout << "\"num_nodes\": " << stats_.num_nodes_ << "," << endl;
  cout << "\"num_nodes_left\": " << stats_.num_nodes_left_ << "," << endl;
  cout << "\"num_sec\": " << stats_.num_sec_ << "," << endl;
  cout << "\"num_sec_separations\": " << stats_.num_sec_separations_ << ","
       << endl;
  cout << "\"num_sec_root\": " << stats_.num_sec_root_ << "," << endl;
  cout << "\"num_sec_separations_root\": " << stats_.num_sec_separations_root_
       << "," << endl;
  cout << "\"num_vars\": " << stats_.num_vars_ << "," << endl;
  cout << "\"num_constraints\": " << stats_.num_constraints_ << "," << endl;
  cout << "\"build_arrangement_time\": " << stats_.build_arrangement_time_
       << "," << endl;
  cout << "\"objective_lower_bound\": " << stats_.objective_lower_bound_ << ","
       << endl;
  cout << "\"objective_upper_bound\": " << stats_.objective_upper_bound_ << ","
       << endl;
  cout << "\"num_faces\": " << stats_.num_faces_ << "," << endl;
  cout << "\"num_triineq\": " << stats_.num_triineq_ << "," << endl;
  cout << "\"num_triineq_separations\": " << stats_.num_triineq_separations_
       << "," << endl;
  cout << "\"num_triineq_root\": " << stats_.num_triineq_root_ << "," << endl;
  cout << "\"num_triineq_separations_root\": "
       << stats_.num_triineq_separations_root_ << "," << endl;
  cout << "\"num_benders_separations\": " << stats_.num_benders_separations_
       << "," << endl;
  cout << "\"num_benders_separations_root\": "
       << stats_.num_benders_separations_root_ << "," << endl;
  cout << "\"num_benders_cuts\": " << stats_.num_benders_cuts_ << "," << endl;
  cout << "\"num_benders_cuts_root\": " << stats_.num_benders_cuts_root_ << ","
       << endl;
  cout << "\"num_runs_heuristic\": " << stats_.num_runs_heuristic_ << endl;

  cout << "}" << endl;
}
