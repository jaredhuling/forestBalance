test_that("every solve is checked and truncated CG fails visibly", {
  set.seed(712); A <- rep(0:1, each=30); lm <- matrix(sample(1:6,60*20,TRUE),60,20)
  Z <- leaf_node_kernel_Z(lm); K <- Matrix::tcrossprod(Z)/20
  for (target in c("ATE","ATT","ATC")) {
    d <- kernel_balance(A,kern=K,solver="direct",lambda=.1,estimand=target,tol=1e-10)
    for (sol in c("cg","bj")) {
      f <- kernel_balance(A,Z=Z,leaf_matrix=lm,num.trees=20,solver=sol,lambda=.1,estimand=target,tol=1e-10)
      expect_equal(f$weights,d$weights,tolerance=1e-6)
      expect_length(f$diagnostics$solves,if(target=="ATE")4 else 2)
      expect_true(all(vapply(f$diagnostics$solves,`[[`,logical(1),"converged")))
    }
  }
  expect_error(kernel_balance(A,Z=Z,num.trees=20,solver="cg",lambda=.1,maxiter=1),"did not converge")
  expect_error(kernel_balance(A,Z=Z,num.trees=20,lambda=Inf),"nonnegative")
})
test_that("both forests and complete diagnostics retain outcome exclusion", {
  set.seed(716); n <- 160; X <- matrix(rnorm(n*8),n,8); A <- rbinom(n,1,.5); Y <- X[,1]+rnorm(n)
  run <- function(y,kr) {
    set.seed(717); forest_balance(X,A,y,num.trees=100,min.node.size=5,solver="cg",lambda=.1,
      tol=1e-10,kernel.response=kr,crossfit.balance="full",num.threads=1,seed=718)
  }
  for (kr in c("joint","treatment","outcome")) {
    f <- run(Y,kr); yy<-Y; yy[f$fold_ids==1]<-yy[f$fold_ids==1]+100
    g <- run(yy,kr)
    expect_length(f$forests,2);expect_length(f$fold_diagnostics,2)
    expect_equal(f$fold_diagnostics[[1]]$balance_weights,g$fold_diagnostics[[1]]$balance_weights,tolerance=1e-12)
    expect_equal(f$weights[f$fold_ids==1],g$weights[g$fold_ids==1],tolerance=1e-12)
    expect_equal(f$ate,mean(f$fold_ates),tolerance=1e-12)
    expect_true(all(vapply(f$fold_diagnostics,function(z)length(z$diagnostics$solves)==4,logical(1))))
  }
  expect_error(forest_balance(X,A,Y,lambda=.1,num.folds=1.5),"num.folds")
  ids<-rep(1:2,each=n/2); aa<-as.numeric(ids==2)
  expect_error(forestBalance:::.fit_one_fold(1,ids,X,aa,Y,100,5,FALSE,NULL,TRUE,"ATE","cg",.1,"joint",1e-10),"Empty")
  expect_error(forestBalance:::.compute_ate(Y,A,rep(0,n),FALSE),"arm mass")
})

test_that("sparse finiteness guards do not allocate a dense matrix", {
  z<-Matrix::sparseMatrix(i=1L,j=1L,x=1,dims=c(100000L,100000L))
  expect_true(forestBalance:::.matrix_all_finite(z))
  z@x[1]<-Inf;expect_false(forestBalance:::.matrix_all_finite(z))
})
