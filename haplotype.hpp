#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <iostream>
#include <vector>

#include "dummy_gPBWT.hpp"

using namespace std;

// Dummy vector in place of h_iv from xg
vector<int64_t> h_iv;

//  RRMemo functions
//  Created by Jordan Eizenga on 6/21/16.
struct RRMemo {
private:

  std::vector<double> S_multipliers;
  double T_multiplier;

  std::vector< std::vector<double> > S;
  std::vector<double> T;

  double rho;
  double exp_rho;

  double S_value(int height, int width);
  double T_value(int width);

public:
  RRMemo(double recombination_penalty);
  ~RRMemo(void);

  double recombination_penalty();

  double rr_diff(int height, int width);
  double rr_same(int height, int width);
  double rr_adj(int width);
  double rr_all(int height, int width);
};

class rectangle {
private:
  // Search interval in B[] corresponding to strip
  // index of this strip's first node within A. We don't use this yet, but we
  // will when start talking about edits to haplo_d's
public:
  ThreadSearchState state;
  int a_index;
  // This lets us find I^a_b-1; R^a_b-1 without having to keep indices
  // consistent between cross-sections
  rectangle* prev = NULL;
  int J;
  int I;
  double R = 0;
  inline int get_next_J(int64_t next_id);
  inline void extend(int64_t next_id);
};

// A cross-section is a column of rectangles S^*_b, * <= b. Each "rectangle" in
// the sense of recomb-rectangle functions is a whole cross_section
struct cross_section {
private:
  int64_t id;
  int b_index;
public:
  vector<rectangle> S;
  int height;
  int width = 1;

  cross_section(int64_t new_height,int i,int64_t new_id);
  inline int64_t get_id();
};

// A haplo_d indexes |A| + 1 columns of rectangles S^*_b according in A-order
class haplo_d {
public:
  vector<cross_section> cs;
  haplo_d(vector<int64_t> h);
  // Needs to be called before the cross_sections have I values in their rectangles
  void calculate_Is(vector<int64_t> h);
  double probability(double recombination_penalty);
};
