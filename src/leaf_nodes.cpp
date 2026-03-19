#include <Rcpp.h>
#include <unordered_map>
using namespace Rcpp;

// [[Rcpp::export]]
IntegerMatrix get_leaf_nodes_cpp(List forest, NumericMatrix newdata) {
  int n = newdata.nrow();
  int ntrees = as<int>(forest["_num_trees"]);
  IntegerMatrix result(n, ntrees);

  List root_nodes   = forest["_root_nodes"];
  List child_nodes  = forest["_child_nodes"];
  List split_vars_l = forest["_split_vars"];
  List split_vals_l = forest["_split_values"];
  List sml_list     = forest["_send_missing_left"];

  for (int tr = 0; tr < ntrees; tr++) {
    int root = as<int>(root_nodes[tr]);
    List children = child_nodes[tr];
    IntegerVector left  = children[0];
    IntegerVector right = children[1];
    IntegerVector svars = split_vars_l[tr];
    NumericVector svals = split_vals_l[tr];
    LogicalVector sml   = sml_list[tr];

    for (int i = 0; i < n; i++) {
      int node = root;

      // Traverse from root to leaf following grf's Tree.cpp::find_leaf_node
      while (left[node] != 0 || right[node] != 0) {
        int sv = svars[node];
        double sval = svals[node];
        double x_val = newdata(i, sv);
        bool x_na = NumericVector::is_na(x_val);
        bool s_na = NumericVector::is_na(sval);
        bool go_left;

        if (!x_na && !s_na) {
          go_left = (x_val <= sval);
        } else if (x_na && sml[node]) {
          go_left = true;
        } else if (x_na && s_na) {
          go_left = true;
        } else {
          go_left = false;
        }

        node = go_left ? left[node] : right[node];
      }

      // Store 1-based node index to match R convention
      result(i, tr) = node + 1;
    }
  }

  return result;
}


// [[Rcpp::export]]
IntegerMatrix remap_leaves_cpp(IntegerMatrix leaf_matrix) {
  int n = leaf_matrix.nrow();
  int B = leaf_matrix.ncol();
  IntegerMatrix out(n, B);
  int offset = 0;

  for (int b = 0; b < B; b++) {
    std::unordered_map<int, int> leaf_map;
    int next_id = 1;

    for (int i = 0; i < n; i++) {
      int leaf = leaf_matrix(i, b);
      auto it = leaf_map.find(leaf);
      if (it == leaf_map.end()) {
        leaf_map[leaf] = next_id;
        next_id++;
      }
      out(i, b) = leaf_map[leaf] + offset;
    }

    offset += next_id - 1;
  }

  return out;
}
