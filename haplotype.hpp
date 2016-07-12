#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <iostream>
#include <vector>

using namespace std;

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
  ThreadSearchState state;
  int a_index;
  // This lets us find I^a_b-1; R^a_b-1 without having to keep indices
  // consistent between cross-sections
public:
  rectangle* prev = nullptr;
  int J = 0;
  int I = 0;
  double R = 0;
  int get_next_J(int64_t next_id);
  void extend(int64_t next_id);
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
  cross_section(int64_t new_height,int b,int64_t new_node_id);
  inline int64_t get_id();
};

// A haplo_d indexes |A| + 1 columns of rectangles S^*_b according in A-order
class haplo_d {
public:
  rectangle empty_rect;
  vector<cross_section> cs;
  haplo_d(const thread_t& t);
  // Needs to be called before the cross_sections have I values in their rectangles
  void calculate_Is();
  double probability(double recombination_penalty);
};
