
df <- with(mtcars, data.frame(y=mpg, x1=disp, x2=hp, x3=wt))

nll_lm <- function() {
  y <- data$y
  X <- as.matrix(data[, -1]) 
  beta <- par[1:(ncol(X))]  
  sigma <- par[ncol(X) + 1]  
  
  if (sigma <= 0) {
    return(Inf) 
  }
  
  residuals <- y - X %*% beta
  
  if (any(is.na(residuals)) || any(is.infinite(residuals))) {
    return(Inf)
  }
  
  nll <- (length(y) / 2) * log(2 * pi * sigma^2) + sum(residuals^2) / (2 * sigma^2)
  
  if (!is.finite(nll)) {
    return(Inf)  # If nll is not finite, return a large value to indicate failure
  }
  
  return(nll)
}

initial_guess <- c(rep(mean(df$y), ncol(df) - 1), 1)  

optim_result <- optim(par = initial_guess, fn = nll_lm, data = df, hessian = TRUE, control = list(fnscale = -1), 
                      lower = c(rep(-Inf, ncol(df)-1), 0.0001),  
                      upper = c(rep(Inf, ncol(df)-1), 100))  

if (optim_result$convergence != 0) {
  stop("Optimization did not converge. Try adjusting the initial guesses or optimization settings.")
}

hessian_matrix <- optim_result$hessian

hessian_det <- det(hessian_matrix)
if (abs(hessian_det) < 1e-6) {
  stop("Hessian matrix is singular or near-singular. The standard errors cannot be computed.")
}

cov_matrix <- solve(hessian_matrix)

se <- sqrt(diag(cov_matrix))  

cat("Standard errors of the regression coefficients:", se, "\n")

lm_fit <- lm(y ~ x1 + x2 + x3, data = df)  

lm_coefficients <- coef(lm_fit)
lm_sigma <- summary(lm_fit)$sigma
cat("Residual standard deviation (sigma) from lm():", lm_sigma, "\n")
