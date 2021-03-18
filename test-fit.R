test_that("basics", {
  
  fit <- blblm(mpg ~ wt  +   drat , data = mtcars, m = 3, B = 150)
  expect_s3_class(fit, "blblm")
  expect_equal(length(coef(fit)), 3)
  
    
  fit <- blblm(mpg ~ wt  +  drat , data = mtcars, m = 3, B = 400)
  expect_s3_class(fit, "blblm")
  expect_equal(dim(confint(fit, c("wt","drat"))), c(2,2))
  
  
})
