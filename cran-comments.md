## Update to 'forestBalance' -- 0.1.1

* Compatibility with Eigen 5: `src/cg_solve.cpp` uses `Eigen::Map<Eigen::SparseMatrix>`
  in place of the removed `Eigen::MappedSparseMatrix` (contributed by Dirk Eddelbuettel
  ahead of the RcppEigen update).
* New estimands (ATT, ATC), a ridge penalty (default 0.1), and solver diagnostics and
  checks; see NEWS.md.

## Test environments

* local macOS Sequoia (R 4.5.1), also against the development RcppEigen
  (0.4.9.9-5, Eigen 5.0.1)
* R-hub: linux (R-devel), windows (R-devel), macos-arm64 (R-devel), clang-asan
* win-builder: R-devel (2026-11-07 r90850) and R-release (4.6.2)

## R CMD check results

R-hub and win-builder: 0 errors | 0 warnings | 0 notes.

Local check: 0 errors | 0 warnings | 2 notes, both from the local environment (the current
time could not be verified, and the local HTML Tidy is too old for manual validation).
