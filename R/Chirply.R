rm(list = ls())

# Define the f function
chirp_est<- function(y_t) {
  n <- length(y_t)
  t <- 1:n

  # Generate grid for parameters
  s1 <- seq(0, (n - 1) * pi / n, length.out = n - 1)
  s2 <- seq(0, (n^2 - 1) * pi / n^2, length.out = n^2 - 1)

  # Define the function to minimize
  f1 <- function(theta) {
    a <- theta[1]
    b <- theta[2]
    return((-1 / n) * ((sum(y_t * cos(a * t + b * t^2)))^2 + (sum(y_t * sin(a * t + b * t^2)))^2))
  }

  # Grid search for initial estimates
  param_grid <- expand.grid(a = s1, b = s2)
  param_grid$f1_values <- apply(param_grid, 1, f1)
  max_param <- param_grid[which.min(param_grid$f1_values), ]
  alpha_int <- as.numeric(max_param[1])
  beta_int <- as.numeric(max_param[2])

  # Refine estimates using optimization
  res_1 <- optim(par = c(alpha_int, beta_int), fn = f1, method = "BFGS")
  alpha.upd <- res_1$par[1]
  beta.upd <- res_1$par[2]

  thres <- 0.0001
  max.iter <- 15
  iter <- 1
  alpha.upd_1 <- alpha_int
  beta.upd_1 <- beta_int

  # Iterative optimization until convergence
  while (abs(alpha.upd - alpha.upd_1) > thres && abs(beta.upd - beta.upd_1) > thres && iter <= max.iter) {
    result <- optim(par = c(alpha.upd, beta.upd), fn = f1, method = "BFGS")
    alpha.upd_1 <- alpha.upd
    beta.upd_1 <- beta.upd
    alpha.upd <- result$par[1]
    beta.upd <- result$par[2]
    iter <- iter + 1
  }

  # Final parameter estimates
  alpha_hat <- alpha.upd
  beta_hat <- beta.upd
  A_hat <- (2 / n) * sum(y_t * cos(alpha_hat * t + beta_hat * t^2))
  B_hat <- (2 / n) * sum(y_t * sin(alpha_hat * t + beta_hat * t^2))
  y_hat <- A_hat * cos(alpha_hat * t + beta_hat * t^2) + B_hat * sin(alpha_hat * t + beta_hat * t^2)
  sigma.sq_hat <- mean((y_t - y_hat)^2)
  y_t <- y_t - y_hat

  return(list(
    est = c(alpha_int, beta_int, alpha_hat, beta_hat, A_hat, B_hat, sigma.sq_hat),
    res = y_t
  ))
}



chirp_model_fit <- function(d, k) {
  # Initialize variables
  y <- d
  es <- c()

  # Perform iterative fitting for k components
  for (n in 1:k) {
    out <- chirp_est(y)  # Call the estimation function
    es <- cbind(es, out$est)  # Collect the estimates
    y <- out$res  # Update residuals for the next iteration
  }



  # Prepare data for fitting
  data <- d - mean(d)  # Center the data
  t <- 1:length(data)

  # Compute fitted values for each component
  c <- matrix(0, nrow = length(data), ncol = k)
  for (i in 1:k) {
    c[, i] <- es[5, i] * cos(es[3, i] * t + es[4, i] * t^2) +
      es[6, i] * sin(es[3, i] * t + es[4, i] * t^2)
  }

  # Sum fitted values for all components
  y_f_hat <- rowSums(c)

  # Residual analysis
  r <- data - y_f_hat

  # Calculate mean RSS
  rss_m <- mean(r^2)

  # Stationarity tests
  t1 <- kpss.test(r)  # KPSS test
  t2 <- adf.test(r)   # ADF test

  # Return results as a list
  return(list(
    rss_mean = rss_m,
    kpss_test = t1,
    adf_test = t2
  ))
}


