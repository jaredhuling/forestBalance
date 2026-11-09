## Update to 'forestBalance' -- 0.1.1

* Compatibility with Eigen 5: `src/cg_solve.cpp` uses `Eigen::Map<Eigen::SparseMatrix>`
  in place of the removed `Eigen::MappedSparseMatrix` (contributed by Dirk Eddelbuettel
  ahead of the RcppEigen update).
* New estimands (ATT, ATC), a ridge penalty (default 0.1), and solver diagnostics and
  checks; see NEWS.md.

## Test environments

* local macOS Sequoia (R 4.5.1)

## R CMD check results

0 errors | 0 warnings | 2 notes (local environment only: the current time could not be
verified, and the local HTML Tidy is too old for manual validation)
