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

cross_section::cross_section(int64_t new_height,int i,int64_t new_id) {
  b_index = i;
  height = new_height;
  id = new_id;
}

inline int64_t cross_section::get_id() {
  return id;
}

inline int rectangle::get_next_J(int64_t next_id) {
  // TO DO: make this work locally
  if(state.is_empty()) { // When we start a new thread in current_threads_to_extend we make a new
    // ThreadSearchState. Set its interval to the entirety of the starting node
    //cerr << "Starting a new strip at " << next_id - 1 << endl;
    state.range_start = 1;
    // int64_t next_side = id_to_rank(next_id)*2 + next_is_reverse;
    // state.range_end = h_iv[next_side];
    state.range_end = pt.h_iv(next_id-1);
  } else {
    // Not brand new so extend it
    cerr << "extending " << state.range_start << " and ";
    state.range_start = pt.where_to(state.range_start, next_id);
    cerr << state.range_end << endl;
    state.range_end = pt.where_to(state.range_end, next_id);
  }
  if(state.range_start == 0) {
    return 0;
  } else {
    return state.count() + 1;
  }
}

inline void rectangle::extend(int64_t next_id) {
  state.range_start = pt.where_to(state.range_start, next_id);
  state.range_end = pt.where_to(state.range_end, next_id);
}

haplo_d::haplo_d(vector<int64_t> h) {
  rectangle rect;
  rect.J = pt.h_iv(0);
  rect.I = rect.J;
  int last_height = rect.J;
  rect.state.range_start = 1;
  rect.state.range_end = pt.h_iv(h[0]);
  cs.push_back(cross_section(rect.J,0,h[0]));
  cs.back().S.push_back(rect);
  cs.back().S.back().prev = &empty_rect;
  int width = 0;
  int new_height;
  bool add_rectangle;
  bool add_A;
  for(int i = 1; i < h.size(); i++) {
    //cerr << "node " << i << ":" << endl;
    width++;
    new_height = pt.h_iv(h[i]);
    rect = cs.back().S[0];
    rect.J = rect.get_next_J(h[i]); // step this strip forward
    if(last_height > rect.J) {
      //cerr << "last_height > rect.J" << endl;
      add_A = 1;
    }
    if(rect.J < new_height) {
      //cerr << "rect.J < new_height" << endl;
      add_rectangle = 1;
      add_A = 1;
    }
    if(add_A) {
      cerr << "adding A" << endl;
      cs.back().width = width;
      width = 0;
      cs.push_back(cross_section(new_height,i,h[i]));
    }
    if(add_rectangle) {
      rectangle new_rect;
      new_rect.state.range_start = 1;
      new_rect.state.range_end = new_height;
      new_rect.J = new_height;
      cs.back().height = new_rect.J;
      cs.back().S.push_back(new_rect);
      cs.back().S.back().I = new_rect.J - rect.J;
      cs.back().S.back().prev = &empty_rect;
      cerr << "at " << i << " added new rectangle with J =" << cs.back().S.back().J << " , I = " << cs.back().S.back().J << endl;
    }
    if(add_A) {
      cs.back().S.push_back(rect);
      cs.back().S.back().prev = &(cs.end()[-2].S.back());
      cerr << "at " << i << " added new rectangle with J =" << cs.back().S.back().J << endl;
    }
    last_height = new_height;
    add_A = 0;
    add_rectangle = 0;
  }
}

void haplo_d::calculate_Is(vector<int64_t> h) {
  // node 1 is done
  for(int b = 1; b < cs.size(); b++) {
    int64_t next_id = cs[b].get_id();
    bool nonempty_J = 1;
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
    cs[b].S.back().I = cs[b].S.back().J;
  }
  for(int i = 0; i < cs.size(); i++) {
    cerr << "at node " << i;
    for(int j = 0; j < cs[i].S.size(); j++) {
       cerr <<", J = " << cs[i].S[j].J;
       cerr <<", I = " << cs[i].S[j].I;
    }
    cerr << " & ht = " << cs[i].height << " & width = " << cs[i].width << endl;
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
  cerr << "R so far " << cs[0].S[0].R << endl;
  for (int b = 1; b < cs.size(); b++) {
    S1 = 0;
    S1S2 = 0;
    for(int a = 0; a < cs[b].S.size(); a++) {
      // N.B. that R's are r^a_b's rather that R^a_b's. Thus the I factor
      cerr << "previous R^" << a << "_" << b - 1 << " is " << cs[b].S[a].prev->R << endl;
      S1 += (cs[b].S[a].prev->R) * (cs[b].S[a].I);
      cerr << "S1 is " << S1 << endl;
    }
    for(int a = 0; a < cs[b-1].S.size(); a++) {
      S1S2 += (cs[b].S[a].prev->R) * (cs[b].S[a].prev->I);
      cerr << " and S1S2 is " << S1S2 << endl;
    }
    // calculate contributions from all continuing strips
    for(int a = 0; a < cs[b].S.size(); a++) {
      cs[b].S[a].R =
      ((1 - memo.recombination_penalty()) * (S1 * memo.rr_diff(cs[b].height, cs[b].width)) +
      ((cs[b].S[a].prev->R) * memo.rr_adj(cs[b].width)) +
      (memo.recombination_penalty() * S1S2 * memo.rr_all(cs[b].height,cs[b].width)));
    }
    // calculate contribution of the new strip
    //cs[b].S.back().R = memo.recombination_penalty() * S1S2 * memo.rr_all(cs[b].height,cs[b].width);
    cerr << "R so far " << cs[b].S.back().R << endl;
  }
  double total_probability_haplotype = 0;
  for(int a = 0; a < cs.back().S.size(); a++) {
    total_probability_haplotype += cs.back().S[a].R;
  }
  return total_probability_haplotype;
}

// tests
// TO DO: tests for haplo stuff

projected_thread array_to_projected_thread(vector<vector<int> > test_array) {
  projected_thread thread;
  thread.h_nodes = vector<projected_node>(test_array[0].size());
  for(int i = 0; i < test_array.size(); i++) {
    thread.h_nodes[i].node = i;
    for(int j = 0; j < test_array[i].size() - 1; j++) {
      int a = test_array[i][j];
      int b = test_array[i][j+1];
      int left = a - b;
      int tot = thread.h_nodes[i+j].next_ranks.back();
      for(int k = 1; k <= b; k++) {
        thread.h_nodes[i+j].next_ranks.push_back(tot + k);
      }
      for(int l = 0; l < left; l++) {
        thread.h_nodes[i+j].next_ranks.push_back(thread.h_nodes[i+j].next_ranks.back());
      }
    }
    for(int j = 0; j < test_array[i].back(); j++) {
      thread.h_nodes.back().next_ranks.push_back(0);
    }
  }
  return thread;
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

  vector<vector<int> > test_down
    { {4,3,3,3,3,2,1}
    };
  // Should have I = 4, 3, 2, 1, 0
  vector<int64_t> h_down {0,1,2,3,4,5,6};
  vector<int> h_iv {4,3,3,3,3,2,1};
  pt = array_to_projected_thread(test_down);
  cerr << "made the mock gPBWT" <<endl;
  for(int i = 0; i < pt.h_nodes.size(); i++) {
    cerr << "next_ranks @ node" << i << " = { ";
    for(int j = 0; j < pt.h_nodes[i].next_ranks.size(); j++) {
      cerr << pt.h_nodes[i].next_ranks[j] << " ";
    }
    cerr << "}" << endl;
  }
  haplo_d hd_down = haplo_d(h_down);
  hd_down.calculate_Is(h_down);
  double prob_down = hd_down.probability(0.2);
  cerr << "Calculated " << prob_down << " as P(h|G,H) for test case 'down'" << endl;

  vector<vector<int> > test_up
    { {3,3,3,3,3},
        {4,4,4,4},
          {5,5,5},
            {6,6},
              {7}
    };

  vector<int64_t> h_up {0,1,2,3,4};
  h_iv = {3,7,12,18,25};
  pt = array_to_projected_thread(test_up);

  for(int i = 0; i < pt.h_nodes.size(); i++) {
    cerr << "next_ranks @ node" << i << " = { ";
    for(int j = 0; j < pt.h_nodes[i].next_ranks.size(); j++) {
      cerr << pt.h_nodes[i].next_ranks[j] << " ";
    }
    cerr << "}" << endl;
  }

  haplo_d hd_up = haplo_d(h_up);
  hd_up.calculate_Is(h_up);
  double prob_up = hd_up.probability(0.2);
  cerr << "Calculated " << prob_up << " as P(h|G,H) for test case 'up'" << endl;


  vector<vector<int> > test_tri
    { {5,5,5},
        {2,0}
    };
  vector<int64_t> h_tri {0,1,2};
  h_iv = {5,7,5};
  pt = array_to_projected_thread(test_tri);
  for(int i = 0; i < pt.h_nodes.size(); i++) {
    cerr << "next_ranks @ node" << i << " = { ";
    for(int j = 0; j < pt.h_nodes[i].next_ranks.size(); j++) {
      cerr << pt.h_nodes[i].next_ranks[j] << " ";
    }
    cerr << "}" << endl;
  }

  haplo_d hd_tri = haplo_d(h_tri);
  hd_tri.calculate_Is(h_tri);
  double prob_tri = hd_tri.probability(0.2);
  cerr << "Calculated " << prob_tri << " as P(h|G,H) for test case 'up'" << endl;



  vector<vector<int> > test_switch
    { {8,0},
        {4}
    };
  vector<int64_t> h_switch {0,1};
  h_iv = {8,4};
  pt = array_to_projected_thread(test_switch);
  for(int i = 0; i < pt.h_nodes.size(); i++) {
    cerr << "next_ranks @ node" << i << " = { ";
    for(int j = 0; j < pt.h_nodes[i].next_ranks.size(); j++) {
      cerr << pt.h_nodes[i].next_ranks[j] << " ";
    }
    cerr << "}" << endl;
  }

  haplo_d hd_switch = haplo_d(h_switch);
  hd_switch.calculate_Is(h_switch);
  double prob_switch = hd_switch.probability(0.2);
  cerr << "Calculated " << prob_switch << " as P(h|G,H) for test case 'up'" << endl;

  vector<vector<int> > leave_together
    { {5,4,3},
        {6,5},
          {7}
    };

  vector<vector<int> > test_snake
    { {4,4,4,0,0},
        {3,3,3,0},
          {2,2,2}
    };
}
