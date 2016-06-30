#include "haplotype.hpp"

using namespace std;

projected_thread pt;

RRMemo::RRMemo(double recombination_penalty)  {
  rho = recombination_penalty;
  exp_rho = exp(-rho);
  S.push_back(std::vector<double>(1, 1.0));
  S_multipliers.push_back(1.0);
  T.push_back(1.0);
  T_multiplier = 1.0 - exp_rho;
}

RRMemo::~RRMemo(void) {

}

double RRMemo::recombination_penalty() {
  return rho;
}

double RRMemo::S_value(int height, int width) {

  while (S.size() < height) {
    S_multipliers.push_back(S_multipliers[S_multipliers.size() - 1] + exp_rho);
    S.push_back(std::vector<double>(1, 1.0));
  }
  std::vector<double>& S_row = S[height - 1];
  double S_multiplier = S_multipliers[height - 1];

  while (S_row.size() < width) {
    S_row.push_back(S_row[S_row.size() - 1] * S_multiplier);
  }

  return S_row[width - 1];
}

double RRMemo::T_value(int width) {

  while (T.size() < width) {
    T.push_back(T[T.size() - 1] * T_multiplier);
  }

  return T[width - 1];
}

double RRMemo::rr_diff(int height, int width) {

  if (height < 1 || width < 1) {
    cerr << "error:[RRMemo] height and width of recombination rectangle must be >= 1" << endl;
  }

  return (S_value(height, width) - T_value(width)) / height;
}

double RRMemo::rr_same(int height, int width) {

  if (height < 1 || width < 1) {
    cerr << "error:[RRMemo] height and width of recombination rectangle must be >= 1" << endl;
  }

  double T_val = T_value(width);
  return (S_value(height, width) - T_val) / height + T_val;
}

double RRMemo::rr_adj(int width) {

  if (width < 1) {
    cerr << "error:[RRMemo] height and width of recombination rectangle must be >= 1" << endl;
  }

  return T_value(width);
}

double RRMemo::rr_all(int height, int width) {

  if (height < 1 || width < 1) {
    cerr << "error:[RRMemo] height and width of recombination rectangle must be >= 1" << endl;
  }

  return exp_rho * S_value(height, width);
}

// unmemoized implementations

double rr_diff(int height, int width, double recombination_penalty) {
  double exp_rho = exp(-recombination_penalty);
  return (pow(1.0 + (height - 1.0) * exp_rho, width - 1.0) - pow(1.0 - exp_rho, width - 1.0)) / height;
}

double rr_same(int height, int width, double recombination_penalty) {
  double exp_rho = exp(-recombination_penalty);
  double T_val = pow(1.0 - exp_rho, width - 1.0);
  return (pow(1.0 + (height - 1.0) * exp_rho, width - 1.0) - T_val) / height + T_val;
}

double rr_adj(int width, double recombination_penalty) {
  return pow(1.0 - exp(-recombination_penalty), width - 1.0);
}

double rr_all(int height, int width, double recombination_penalty) {
  double exp_rho = exp(-recombination_penalty);
  return exp_rho * pow(1.0 + (height - 1.0) * exp_rho, width - 1.0);
}

cross_section::cross_section(rectangle R,int64_t new_height,int i,int64_t new_id) {
  S.push_back(R);
  b_index = i;
  height = new_height;
  id = new_id;
}

inline int64_t cross_section::get_id() {
  return id;
}

inline int rectangle::get_next_J(int64_t next_id) {
  // TO DO: make this work locally
  if(!state.is_empty()) { // When we start a new thread in current_threads_to_extend we make a new
    // ThreadSearchState. Set its interval to the entirety of the starting node
    state.range_start = 0;
    // int64_t next_side = id_to_rank(next_id)*2 + next_is_reverse;
    // state.range_end = h_iv[next_side];
    state.range_end = pt.h_iv(next_id);
  } else {
    // Not brand new so extend it
    state.range_start = pt.where_to(state.range_start, next_id);
    state.range_end = pt.where_to(state.range_end, next_id);
  }
  return state.count();
}

inline void rectangle::extend(int64_t next_id) {
  state.range_start = pt.where_to(state.range_start, next_id);
  state.range_end = pt.where_to(state.range_end, next_id);
}

haplo_d::haplo_d(vector<int64_t> h) {
  int prev_height = pt.h_iv(h[0]);
  int new_height;
  int curr_width = 0;
  for(int i = 1; i < h.size(); i++) {
    new_height = pt.h_iv(h[0]);
    curr_width++;
    rectangle test_R = cs.back().S[0];
    test_R.J = test_R.get_next_J(h[i]);
    if(prev_height > test_R.J || test_R.J < new_height) {
      cs.back().width = curr_width;
      curr_width = 0;
      cs.push_back(cross_section(test_R,new_height,i,h[i]));
    }
  }
}

// Because of how J-value queries must be done, we get a gPBWT search-state
// extension for free on all except for top-level queries
// This function actually also finishes building the cross_sections
void haplo_d::calculate_Is(vector<int64_t> h) {
  for(int b = 1; b < cs.size(); b++) {
    int64_t next_id = cs[b].get_id();
    int i = 0;
    int new_J = cs[b].S[0].get_next_J(next_id);
    int last_J = cs[b-1].S[0].J;
    bool more_changes = new_J < last_J;
    bool keep_going;
    // TO DO: right now this searches top-level down only--likely the most
    // efficient thing is to search bottom-up as well. Consider alternating.
    // The advantage is that top-level *and* bottom level strips are most
    // likely to change or leave
    while(more_changes) {
      last_J = cs[b-1].S[i].J;
      cs[b].S.push_back(cs[b-1].S[i]);
      cs[b].S.back().prev = &cs[b-1].S[i];
      new_J = cs[b].S.back().get_next_J(next_id);
      if(new_J == 0) {
        keep_going = 0;
        break;
      } else {
        if(new_J == last_J) {
          more_changes = 0;
          break;
        } else {
          (cs[b].S.end()--)->I = cs[b].S.back().J - new_J;
          if(i == cs[b-1].S.size() - 1) {
            more_changes = 0;
            break;
          } else {
            i++;
          }
        }
      }
    }
    if(more_changes) {
      for(i; i < cs[b-1].S.size(); i++) {
        cs[b].S.push_back(cs[b-1].S[i]);
        cs[b].S.back().prev = &cs[b-1].S[i];
        cs[b].S.back().extend(h[i]); //is this necessary?
      }
    }
  }
}

double haplo_d::probability(double recombination_penalty) {
  RRMemo memo = RRMemo(recombination_penalty);
  // defined same as in writeup
  double S1 = 0;
  double S1S2 = 0;
  // compute R for the first interval (which has no predecessor)
  // we are always working at the left edge of a cross_section
  cs[0].S[0].R = memo.rr_all(cs[0].height,cs[0].width);
  for (int b = 1; b < cs.size(); b++) {
    S1 = 0;
    S1S2 = 0;
    for(int a = 0; a < cs[b].S.size(); a++) {
      // N.B. that R's are r^a_b's rather that R^a_b's. Thus the I factor
      S1 += (cs[b].S[a].prev->R) * (cs[b].S[a].I);
      S1S2 += (cs[b].S[a].prev->R) * (cs[b].S[a].prev->I);
    }
    // calculate contributions from all continuing strips
    for(int a = 0; a < (cs[b].S.size() -1); a++) {
      cs[b].S[a].R =
      ((1 - memo.recombination_penalty()) * (S1 * memo.rr_diff(cs[b].height, cs[b].width)) +
      ((cs[b].S[a].prev->R) * memo.rr_adj(cs[b].width)) +
      (memo.recombination_penalty() * S1S2 * memo.rr_all(cs[b].height,cs[b].width)));
    }
    // calculate contribution of the new strip
    cs[b].S.back().R = memo.recombination_penalty() * S1S2 * memo.rr_all(cs[b].height,cs[b].width);
  }
  double total_probability_haplotype = 0;
  for(int a = 0; a < cs.back().S.size(); a++) {
    total_probability_haplotype += cs.back().S[a].R;
  }
  return total_probability_haplotype;
}

// tests
// TO DO: tests for haplo stuff

int main(void) {
  // RRMemo tests
  double epsilon = 0.0000001;

  double memo_val;
  double direct_val;

  for (double rho = 1.0; rho < 5.0; rho += 1.0) {

    RRMemo memo(rho);

    for (int c = 1; c < 10; c++) {
      for (int n = 1; n < 10; n++) {

        memo_val = memo.rr_diff(n, c);
        direct_val = rr_diff(n, c, rho);

        if (fabs(memo_val - direct_val) > epsilon) {
          cerr << "FAIL: rr_diff, n = " << n << ", c = " << c << ", rho = " << rho
          << ", memo value = " << memo_val << ", direct value = " << direct_val << endl;
          exit(1);
        }

        memo_val = memo.rr_same(n, c);
        direct_val = rr_same(n, c, rho);

        if (fabs(memo_val - direct_val) > epsilon) {
          cerr << "FAIL: rr_same, n = " << n << ", c = " << c << ", rho = " << rho
          << ", memo value = " << memo_val << ", direct value = " << direct_val << endl;
          exit(1);
        }

        memo_val = memo.rr_all(n, c);
        direct_val = rr_all(n, c, rho);

        if (fabs(memo_val - direct_val) > epsilon) {
          cerr << "FAIL: rr_all, n = " << n << ", c = " << c << ", rho = " << rho
          << ", memo value = " << memo_val << ", direct value = " << direct_val << endl;
          exit(1);
        }
      }

      memo_val = memo.rr_adj(c);
      direct_val = rr_adj(c, rho);

      if (fabs(memo_val - direct_val) > epsilon) {
        cerr << "FAIL: rr_adj, c = " << c << ", rho = " << rho
        << ", memo value = " << memo_val << ", direct value = " << direct_val << endl;
        exit(1);
      }
    }
  }

  cerr << "RR tests passed!" << endl;
}
