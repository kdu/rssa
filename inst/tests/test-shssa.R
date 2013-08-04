library(testthat);
library(Rssa);
source(system.file("extdata", "common.test.methods.R", package = "Rssa"));
context("2dSSA");


test_that("simple 2d-ssa test", {
  F0 <- matrix(c(0, 1, 0, 1, 1, 1), 2,3);
  F <- rbind(cbind(F0,F0), cbind(F0,F0));
  
  # Test simple matrix multiplication
  vtest <- as.vector(F0);
  ht_mul_v <- c(4,2,4,3,2,3,3,2,3,4,2,4);
  utest <- c(0,0,0,0,1,0,1,0,1,0,0,0);
  h_mul_u <- c(3,2,1,3,1,2);
  Frec <- outer(c(1,3,9,27), c(1, 2, 4, 8, 16, 32));
 
  # Test 2D-SSA
  s <-  ssa(F, kind="2d-ssa", force.decompose=FALSE);
  h <- .get.or.create.hbhmat(s);
  expect_equal(round(hbhmatmul(h, vtest, TRUE)), ht_mul_v, 
               label = "tmatmul, 2D-SSA");
  expect_equal(round(hbhmatmul(h, utest)), h_mul_u, label = "matmul, 2D-SSA");
  Frec_res <- .Call("hbhankelize_one_fft", as.vector(Frec[1:2,1:3]), 
                    as.vector(Frec[1:3,1:4]), h);
  expect_equal(as.vector(Frec), round(Frec_res), label = "hankelize, 2D-SSA");
  
  # Test R-SH-SSA
  s2 <- ssa(F, kind="shaped2d-ssa",umask=matrix(TRUE,2,3), force.decompose=FALSE);
  h2 <- .sh2hbhmat(s2);
  
  expect_equal(round(ematmul(h2, vtest, TRUE)), ht_mul_v, 
               label = "tmatmul, SH-SSA");
  expect_equal(round(ematmul(h2, utest)), h_mul_u, label = "matmul, SH-SSA");

  Frec_res <- .hankelize.one.shaped2d.ssa(s2, as.vector(Frec[1:2,1:3]), 
                    as.vector(Frec[1:3,1:4]));
  expect_equal(as.vector(Frec), round(as.vector(Frec_res)), label = "hankelize, SH-SSA");

  # Test 2D-SH-SSA (hbhmat with call to 2D-SSA C routine)
  s3 <-  ssa(F, kind="2d-ssa",umask=matrix(TRUE,2,3), force.decompose=FALSE);
  h3 <- .get.or.create.hbhmat(s3);
  expect_equal(round(hbhmatmul(h3, vtest, TRUE)), ht_mul_v, 
               label = "tmatmul, 2D-SSA");
  expect_equal(round(hbhmatmul(h3, utest)), h_mul_u, label = "matmul, 2D-SH-SSA");
  Frec_res <- .Call("hbhankelize_one_fft", as.vector(Frec[1:2,1:3]), 
                    as.vector(Frec[1:3,1:4]), h3);
  expect_equal(as.vector(Frec), round(Frec_res), label = "hankelize, 2D-SH-SSA");
});


test_that("Shaped 2D-SSA test", {
  F0 <- matrix(c(0, 1, 0, 1, 1, 1), 2,3);
  F <- rbind(cbind(F0,F0), cbind(F0,F0));
  F <- rbind(cbind(F,F), cbind(F,F));
  
  fmask <- matrix(TRUE, dim(F)[1], dim(F)[2]);
  fmask[dim(F)[1],1] <- FALSE;
  fmask[dim(F)[1],2] <- FALSE;
  
  # Test simple matrix multiplication
  umask <- matrix(c(FALSE,TRUE, TRUE, TRUE, 
                    TRUE, FALSE, TRUE, FALSE, 
                    FALSE, FALSE, FALSE, TRUE),3,4);
  
  utest <- 1:6;
  vtest <- 1:52;
 
  # Test R-SH-SSA
  s2 <- ssa(F, kind="shaped2d-ssa",mask=fmask,umask=umask, force.decompose=FALSE);
  h2 <- .sh2hbhmat(s2);
  s3 <- ssa(F, kind="2d-ssa",mask=fmask,umask=umask, force.decompose=FALSE);
  h3 <- .get.or.create.hbhmat(s3);
 
  expect_equal(ematmul(h2, utest,TRUE), hbhmatmul(h3, utest, TRUE), 
               label = "tmatmul, SH-SSA <-> 2D-SH-SSA", tolerance=1e-5);
  expect_equal(ematmul(h2, vtest), hbhmatmul(h3, vtest), 
               label = "matmul, SH-SSA <-> 2D-SH-SSA", tolerance=1e-5);
});





