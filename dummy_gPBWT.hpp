#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <vector>

using namespace std;

struct projected_node {
  int64_t node;
  vector<int> next_ranks;
  projected_node();
};

projected_node::projected_node() {
  next_ranks.push_back(0);
}

struct projected_thread{
  vector<projected_node> h_nodes;
  inline int where_to(int current_rank, int64_t node) {
    // reserve 0 for null
    return h_nodes[node-1].next_ranks[current_rank];
  }
  inline int h_iv(int64_t node) {
    return h_nodes[node].next_ranks.size();
  }
};

struct ThreadSearchState {
    // What is the first visit at that side that is selected?
    int range_start = 0;
    // And what is the past-the-last visit that is selected?
    int range_end = 0;

    // How many visits are selected?
    inline int64_t count() {
        return range_end - range_start;
    }

    // Return true if the range has nothing selected.
    inline bool is_empty() {
        return range_end == range_start;
    }
};

projected_thread convert_test_structure(vector<vector<int> > test_array) {
  projected_thread thread;
  thread.h_nodes = vector<projected_node>(test_array.size(),projected_node());
  for(int i = 0; i < test_array.size(); i++) {
    thread.h_nodes[i].node = i;
    for(int j = 0; j < test_array[i].size()-1; j++) {
      int a = test_array[i][j];
      int b = test_array[i][j+1];
      int left = b - a;
      int tot = thread.h_nodes[i+j].next_ranks.back();
      for(int l = 0; l < left; l++) {
        thread.h_nodes[i+j].next_ranks.push_back(0);
      }
      for(int k = 1; k <= b; b++) {
        thread.h_nodes[i+j].next_ranks.push_back(tot + k);
      }
    }
  }
  return thread;
}
