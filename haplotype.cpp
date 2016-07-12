#include "haplotype.hpp"

using namespace std;

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

// Rectangular decomposition building:

cross_section::cross_section(int64_t new_height,int b,int64_t new_node_id) {
  b_index = b;
  height = new_height;
  id = new_id;
}

void rectangle::extend(int64_t next_id) {
  if(state.current_side == 0) {
    state.range_start = 0;
    state.range_end = h_iv[(node_rank_as_entity(next_id) - 1) * 2 + next_is_reverse];
  } else {
    // Else, look at where the path goes to and apply the where_to function to shrink the range down.
    state.range_start = where_to(state.current_side, state.range_start, next_side);
    state.range_end = where_to(state.current_side, state.range_end, next_side);
  }
  state.current_side = next_side;
}

int rectangle::get_next_J(int64_t next_id) {
  extend(next_id);
  return state.count();
}

inline int64_t cross_section::get_id() {
  return id;
}

haplo_d::haplo_d(const thread_t& t) {
  rectangle rect;
  rect.J = rect.get_next_J(t[0].node_id);
  rect.I = rect.J;
  int last_height = rect.J;
  // TODO FIX INITIALIZATION OF CS
  cs.push_back(cross_section(rect.J,0,t[0].node_id));
  cs.back().S.push_back(rect);
  cs.back().S.back().prev = &empty_rect;
  int width = 0;
  int new_height;
  bool add_rectangle;
  bool add_A;
  for(int i = 1; i < t.size(); i++) {
    width += node_length(t[i-1].node_id);
    new_height = h_iv[(node_rank_as_entity(t[i].node_id) - 1) * 2 + t[i].is_reverse];
    rect = cs.back().S[0];
    rect.J = rect.get_next_J(t[i].node_id); // step this strip forward
    if(last_height > rect.J) {
      add_A = 1;
    }
    if(rect.J < new_height) {
      add_rectangle = 1;
      add_A = 1;
    }
    if(add_A) {
      cs.back().width = width;
      width = 0;
      //TODO FIX
      cs.push_back(cross_section(new_height,i,t[i].node_id));
    }
    if(add_rectangle) {
      rectangle new_rect;
      new_rect.extend(t[i].node_id);
      new_rect.J = new_height;
      cs.back().height = new_rect.J;
      cs.back().S.push_back(new_rect);
      cs.back().S.back().I = new_rect.J - rect.J;
      cs.back().S.back().prev = &empty_rect;
    }
    if(add_A) {
      cs.back().S.push_back(rect);
      cs.back().S.back().prev = &(cs.end()[-2].S[0]);
    }
    last_height = new_height;
    add_A = 0;
    add_rectangle = 0;
  }
}

void haplo_d::calculate_Is() {
  // node 1 is done
  for(int b = 1; b < cs.size(); b++) {
    int64_t next_id = cs[b].get_id();
    bool nonempty_J = (cs[b].S.back().J > 0);
    if(nonempty_J) {
      bool change_in_J = 1;
      int new_J;
      int old_J;
      for(int a = 1; a < cs[b-1].S.size(); a++) {
        if(change_in_J) {
          cs[b].S.push_back(cs[b-1].S[a]);
          cs[b].S.back().prev = &cs[b-1].S[a];
          old_J = cs[b].S.back().J;
          new_J = cs[b].S.back().get_next_J(next_id);
          cs[b].S.end()[-2].I = cs[b].S.end()[-2].J - new_J;
          if(old_J == new_J) {
            change_in_J = 0;
          } else if(new_J == 0) {
            change_in_J = 0;
            nonempty_J = 0;
            cs[b].S.pop_back();
          }
        } else if(nonempty_J) {
          cs[b].S.push_back(cs[b-1].S[a]);
          cs[b].S.back().prev = &cs[b-1].S[a];
        }
      }
    } else {
      cs[b].S.pop_back();
    }
    cs[b].S.back().I = cs[b].S.back().J;
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
    }
    for(int a = 0; a < cs[b-1].S.size(); a++) {
      S1S2 += (cs[b].S[a].prev->R) * (cs[b].S[a].prev->I);
    }
    // calculate contributions from all continuing strips
    for(int a = 0; a < cs[b].S.size(); a++) {
      cs[b].S[a].R =
      ((1 - memo.recombination_penalty()) * (S1 * memo.rr_diff(cs[b].height, cs[b].width)) +
      ((cs[b].S[a].prev->R) * memo.rr_adj(cs[b].width)) +
      (memo.recombination_penalty() * S1S2 * memo.rr_all(cs[b].height,cs[b].width)));
    }
  }
  double total_probability_haplotype = 0;
  for(int a = 0; a < cs.back().S.size(); a++) {
    total_probability_haplotype += cs.back().S[a].R;
  }
  return total_probability_haplotype;
}

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
