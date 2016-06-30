#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <vector>

using namespace std;

struct projected_node {
  vector<int> next_ranks;
};

struct projected_thread{
  vector<projected_node> h_nodes;
  inline int where_to(int current_rank, int64_t node) {
    // reserve 0 for null
    return h_nodes[node-1].next_ranks[current_rank-1];
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
