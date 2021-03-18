#' Linear Regression with Little Bag of Bootstraps
#'
#'
#' @importFrom magrittr %>%
#' @import stats
#' @import purrr
#' @import furrr
#' @import future 
#' @importFrom utils capture.output
#' @importFrom parallel makeCluster
#' @details
#' Linear Regression with Little Bag of Bootstraps
#' @docType package
#' @name blblm
"_PACKAGE"


utils::globalVariables("fit")
utils::globalVariables(".")


#' linear regression with little bag of bootstraps
#'
#'
#' @param formula a formula of model.
#' @param data a data frame.
#' @param m data splits, use a default to be 10.
#' @param B bootstraps, use a  default to be 5000.
#' @param parallel logical value indicate use parallelization or not, default is that not use.
#' @return res A list including estimates and sigma.
#'
#' @export
blblm <- function(formula, data, m = 10, B = 5000, parallel = FALSE) {

   data_list <- split_data(data, m)
    
   if(!parallel){ #not use parallelization
      
    estimates <-  map(data_list,
        ~ lm_each_subsample(formula = formula, data = ., n = nrow(.), B = B))
    
    } else { #use parallelization
      
      cl <- parallel::makeCluster(4, outfile = "")
      plan(cluster, workers = cl)
      estimates <- map(data_list,
        ~ lm_each_subsample(formula = formula, data = .x, n = nrow(.x), B = B))
    }
  
  res <- list(estimates = estimates, formula = formula)
  class(res) <- "blblm"
  invisible(res)
}


#split data into m parts of approximated equal sizes
split_data <- function(data, m) {
  idx <- sample.int(m, nrow(data), replace = TRUE)
  data %>% split(idx)
}


#compute the estimates
lm_each_subsample <- function(formula, data, n, B) {
  replicate(B, lm_each_boot(formula, data, n), simplify = FALSE)
}


#compute the regression estimates for a blb dataset
lm_each_boot <- function(formula, data, n) {
  freqs <- rmultinom(1, n, rep(1, nrow(data)))
  lm1(formula, data, freqs)
}



#A improved Cpp version of lm1
lm1 <- function(formula, data, freqs) {

  environment(formula) <- environment()
  n <- nrow(data)
  X <-   model.matrix(formula, data[rep(1:n, freqs),])
  Y <-   model.response(model.frame(formula, data[rep(1:n, freqs),]))
  fit <- lmcpp(X, Y,
               freqs
               )
 res <- blbcoef(fit)
  coef <- res[,1]
  names(coef) <- colnames(X)
  result <- list(coef = coef, sigma = blbsigma(fit))
  return(result)
}


blbcoef <- function(fit) {
  return(fit$coefficients)
}

blbsigma <- function(fit) {
  p <- fit$rank
  e <- fit$residuals 
  w <- fit$weights
  return(sqrt(sum(w * (e^2)) / (sum(w) - p)))
}


#' Generalized linear regression with little bag of bootstraps
#'
#'
#' @param formula a formula of model.
#' @param family family for GLM  
#' @param data a data frame.
#' @param m data splits, default to be 10.
#' @param B bootstraps, default to be 5000.
#' @param parallel logical value indicate use parallelization or not, default is not to use.
#' @return res A list including estimates and sigma.
#'
#' @export
blbglm <- function(formula, family, data,  m = 10, B = 5000, parallel = FALSE) {
  
  data_list <- split_data(data, m)
  
  if(!parallel){ 
    
    estimates <-  map(data_list,
                      ~ glm_each_subsample(formula = formula,family, data = ., n = nrow(.), B = B))
    
  } else { 
    
    cl <- parallel::makeCluster(4, outfile = "")
    plan(cluster, workers = cl)
    estimates <- map(data_list,
                     ~ glm_each_subsample(formula = formula, family, data = .x, n = nrow(.x), B = B))
  }
  
  res <- list(estimates = estimates, formula = formula)
  class(res) <- c("blbglm")
  invisible(res)
}

glm_each_subsample <- function(formula, family, data, n, B) {
  replicate(B, glm_each_boot(formula, family, data, n), simplify = FALSE)
}

glm_each_boot <- function(formula, family, data, n) {
  environment(formula) <- environment()
  glm(formula,  family, data, rmultinom(1, nrow(data), rep(1, nrow(data))))
}


#R version
lm2 <- function(formula, data, freqs) {
  environment(formula) <- environment()
  fit <- lm(formula, data, weights = freqs)
  result <-   list(coef = blbcoef_R(fit), sigma = blbsigma_R(fit))
  return(result)
}


blbcoef_R <- function(fit) {
  return(coef(fit))
}

blbsigma_R <- function(fit) {
  p <- fit$rank
  y <- model.extract(fit$model, "response")
  e <- fitted(fit) - y
  w <- fit$weights
  result <- sqrt(sum(w * (e^2)) / (sum(w) - p))
  return( result )
}

lm_each_subsample_R <- function(formula, data, n, B) {
  replicate(B, lm_each_boot_R(formula, data, n), simplify = FALSE)
}
#compute the regression estimates for a blb dataset
lm_each_boot_R <- function(formula, data, n) {
  freqs <- rmultinom(1, n, rep(1, nrow(data)))
  lm2(formula, data, freqs)
}


#' @export
#' @method print blblm
print.blblm <- function(x, ...) {
  cat("blblm model:", capture.output(x$formula))
  cat("\n")
}


#' @export
#' @method sigma blblm
sigma.blblm <- function(object, confidence = FALSE, level = 0.95, ...) {
  est <- object$estimates
  sigma <- mean(map_dbl(est, ~ mean(map_dbl(., "sigma"))))
  if (confidence) {
    alpha <- 1 - 0.95
    limits <- est %>%
      map_mean(~ quantile(map_dbl(., "sigma"), c(alpha / 2, 1 - alpha / 2))) %>%
      set_names(NULL)
    return(c(sigma = sigma, lwr = limits[1], upr = limits[2]))
  } else {
    return(sigma)
  }
}

#' @export
#' @method coef blblm
coef.blblm <- function(object, ...) {
  est <- object$estimates
  #for bootstrap, there is possible only one observation repeated for small subsample,
  #so use na.rm = TRUE to drop these NA estimates
  map_mean(est, ~ map_cbind(., "coef") %>% rowMeans(na.rm = T))
}

#' @export
#' @method coef blbglm
coef.blbglm <- function(object, ...) {
  est <- object$estimates
  #for bootstrap, there is possible only one observation repeated for small subsample,
  #so use na.rm = TRUE to drop these NA estimates
  do.call("rbind",lapply(est[[1]],coef)) %>% colMeans(na.rm = T)
}

#' @export
#' @method confint blblm
confint.blblm <- function(object, parm = NULL, level = 0.95, ...) {
  if (is.null(parm)) {
    parm <- attr(terms(fit$formula), "term.labels")
  }
  alpha <- 1 - level
  est <- object$estimates
  out <- map_rbind(parm, function(p) {
    map_mean(est, ~ map_dbl(., list("coef", p)) %>% quantile(c(alpha / 2, 1 - alpha / 2), na.rm = T))
  })
  if (is.vector(out)) {
    out <- as.matrix(t(out))
  }
  dimnames(out)[[1]] <- parm
  out
}

#' @export
#' @method predict blblm
predict.blblm <- function(object, new_data, confidence = FALSE, level = 0.95, ...) {
  est <- object$estimates
  X <- model.matrix(reformulate(attr(terms(object$formula), "term.labels")), new_data)
  if (confidence) {
    map_mean(est, ~ map_cbind(., ~ X %*% .$coef) %>%
      apply(1, mean_lwr_upr, level = level) %>%
      t())
  } else {
    map_mean(est, ~ map_cbind(., ~ X %*% .$coef) %>% rowMeans(na.rm = T))
  }
}



mean_lwr_upr <- function(x, level = 0.95) {
  alpha <- 1 - level
  c(fit = mean(x, na.rm = TRUE), quantile(x, c(alpha / 2, 1 - alpha / 2), na.rm = T) %>% set_names(c("lwr", "upr")))
}

map_mean <- function(.x, .f, ...) {
  (map(.x, .f, ...) %>% reduce(`+`)) / length(.x)
}

map_cbind <- function(.x, .f, ...) {
  map(.x, .f, ...) %>% reduce(cbind)
}

map_rbind <- function(.x, .f, ...) {
  map(.x, .f, ...) %>% reduce(rbind)
}
