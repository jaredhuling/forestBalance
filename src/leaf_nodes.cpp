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


// Build sparse CSC matrix Z from leaf assignments, such that K = Z Z^T / B.
// Combines leaf remapping and CSC construction in a single pass,
// avoiding the overhead of R's sparseMatrix() triplet sort.
// [[Rcpp::export]]
List build_Z_cpp(IntegerMatrix leaf_matrix) {
  int n = leaf_matrix.nrow();
  int B = leaf_matrix.ncol();

  // Pass 1: remap leaf IDs per tree and count total columns
  std::vector<std::vector<int>> tree_remap(B);
  std::vector<int> tree_ncols(B);
  int total_cols = 0;

  for (int b = 0; b < B; b++) {
    std::unordered_map<int, int> leaf_map;
    int next_id = 0;
    tree_remap[b].resize(n);
    for (int i = 0; i < n; i++) {
      int leaf = leaf_matrix(i, b);
      auto it = leaf_map.find(leaf);
      if (it == leaf_map.end()) {
        leaf_map[leaf] = next_id;
        tree_remap[b][i] = next_id;
        next_id++;
      } else {
        tree_remap[b][i] = it->second;
      }
    }
    tree_ncols[b] = next_id;
    total_cols += next_id;
  }

  int nnz = n * B;

  // Pass 2: count entries per column (for column pointers)
  IntegerVector p_vec(total_cols + 1, 0);
  {
    int offset = 0;
    for (int b = 0; b < B; b++) {
      for (int i = 0; i < n; i++) {
        p_vec[offset + tree_remap[b][i] + 1]++;
      }
      offset += tree_ncols[b];
    }
  }
  for (int j = 1; j <= total_cols; j++) {
    p_vec[j] += p_vec[j - 1];
  }

  // Pass 3: fill row indices (sorted within each column because we
  // iterate i = 0..n-1 in order, and each column belongs to one tree)
  IntegerVector i_vec(nnz);
  std::vector<int> pos(total_cols);
  for (int j = 0; j < total_cols; j++) pos[j] = p_vec[j];

  {
    int offset = 0;
    for (int b = 0; b < B; b++) {
      for (int i = 0; i < n; i++) {
        int col = offset + tree_remap[b][i];
        i_vec[pos[col]++] = i;  // 0-based row index
      }
      offset += tree_ncols[b];
    }
  }

  return List::create(
    Named("i") = i_vec,
    Named("p") = p_vec,
    Named("x") = NumericVector(nnz, 1.0),
    Named("nrow") = n,
    Named("ncol") = total_cols
  );
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
