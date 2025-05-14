y <- 0:10
freq <- c(49, 8, 5, 8, 7, 6, 4, 5, 4, 3, 1)
n <- sum(freq)  
loglik <- function(params) {
  omega <- params[1]
  theta <- params[2]
  phi <- params[3]
  
if (omega <= 0 || omega >= 1 || theta <= 0 || phi <= 0) {
    return(Inf)
  }
  
  log_probs <- numeric(length(y))
  
  for (i in seq_along(y)) {
    yi <- y[i]
    if (yi == 0) {
      prob <- omega + (1 - omega) * exp(-theta / (1 + phi))
    } else {
      prob <- (1 - omega) * theta * (theta + phi * yi)^(yi - 1) / ((1 + phi)^yi * factorial(yi)) * 
        exp(-(theta + phi * yi) / (1 + phi))
    }
    
    if (prob <= 0) {
      return(Inf)
    }
    
    log_probs[i] <- log(prob)
  }
  
  ll <- sum(freq * log_probs)
  
  return(-ll)  
}

init_params <- c(0.5, 1, 0.5)

#Maximize log-likelihood (hessian = TRUE is added)
fit <- optim(
  par = init_params,
  fn = loglik,
  method = "L-BFGS-B",
  lower = c(0.0001, 0.0001, 0.0001),
  upper = c(0.9999, Inf, Inf),
  hessian = TRUE
)

# Step 5: Results
omega_hat <- fit$par[1]
theta_hat <- fit$par[2]
phi_hat <- fit$par[3]

cat("Estimated parameters:\n")
cat("omega =", omega_hat, "\n")
cat("theta =", theta_hat, "\n")
cat("phi =", phi_hat, "\n")
#######################################
# fitted probabilities
fitted_probs <- numeric(length(y))

for (i in seq_along(y)) {
  yi <- y[i]
  if (yi == 0) {
    prob <- omega_hat + (1 - omega_hat) * exp(-theta_hat / (1 + phi_hat))
  } else {
    prob <- (1 - omega_hat) * theta_hat * (theta_hat + phi_hat * yi)^(yi - 1) /
      ((1 + phi_hat)^yi * factorial(yi)) * exp(-(theta_hat + phi_hat * yi) / (1 + phi_hat))
  }
  fitted_probs[i] <- prob
}
# Observed relative frequencies
observed_probs <- freq / sum(freq)
##Plot observed vs fitted
barplot(
  rbind(observed_probs, fitted_probs),
  beside = TRUE,
  names.arg = y,
  col = c("skyblue", "tomato"),
  legend.text = c("Observed", "Fitted"),
  args.legend = list(x = "topright"),
  xlab = "y",
  ylab = "Probability",
  main = "Observed vs Fitted Probabilities"
)
###################################
# Step 8: Standard errors
if (!is.null(fit$hessian)) {
  varcov <- tryCatch(
    solve(fit$hessian),
    error = function(e) {
      message("Warning: Hessian not invertible.")
      return(NULL)
    }
  )
  
  if (!is.null(varcov)) {
    se <- sqrt(diag(varcov))
    cat("\nStandard errors:\n")
    cat("SE(omega) =", se[1], "\n")
    cat("SE(theta) =", se[2], "\n")
    cat("SE(phi) =", se[3], "\n")
  }
}

############################ 
Fisher_Matrix <- matrix(c(
  369.9195, -5.9568, 19.3975,
  -5.9568, 5.86535, 0.0826,
  19.3975, 0.0826, 36.1404
), nrow = 3, ncol = 3, byrow = TRUE)
rownames(Fisher_Matrix) <- c("omega", "theta", "phi")
colnames(Fisher_Matrix) <- c("omega", "theta", "phi")
Fisher_Matrix_inv <- solve(Fisher_Matrix)
print(Fisher_Matrix_inv)
standard_errors <- sqrt(diag(Fisher_Matrix_inv))
print(standard_errors)

