# Plot for Figure 2 (training and testing MSPE curves for Simulation 1)

library(crf)
library(ggplot2)

generate_data = function() {
  n = 4
  I = 10000
  N = I*n
  f <- function(x) tanh(x)
  x = matrix(NA,nrow=0,ncol=1)
  CovX = matrix(0,n,n) + diag(1,n)
  sqrtCovX = expm::sqrtm(CovX)
  for (i in 1:I) {
    x_ind = matrix(rnorm(n,0,1),n,1)
    x_ind = sqrtCovX %*% x_ind
    x = rbind(x, x_ind)
  }
  epsilon = rnorm(N)
  for (i in 1:I) {
    CovY = matrix(0,n,n)
    for (j1 in 1:4) {
      for (j2 in 1:4) {
        if (j1==j2) CovY[j1,j2] = (1/4+1/(1+exp(4*x[(i-1)*n+j1]))) * (1/4+1/(1+exp(4*x[(i-1)*n+j2])))
        if (j1!=j2) CovY[j1,j2] = 0.8 * (1/4+1/(1+exp(4*x[(i-1)*n+j1]))) * (1/4+1/(1+exp(4*x[(i-1)*n+j2])))
      }
    }
    sqrtCovY = expm::sqrtm(CovY)
    epsilon[((i-1)*n+1):(i*n)] = sqrtCovY %*% epsilon[((i-1)*n+1):(i*n)]
  }
  y = f(x) + epsilon
  rdf <- data.frame(y=y,x=x, id=rep(1:(N/n),each=n))
  return (rdf)
}
generate_train_data = function() {
  n = 4
  I = 10000
  N = I*n
  f <- function(x) tanh(x)
  x = matrix(NA,nrow=0,ncol=1)
  CovX = matrix(0,n,n) + diag(1,n)
  sqrtCovX = expm::sqrtm(CovX)
  for (i in 1:I) {
    x_ind = matrix(rnorm(n,0,1),n,1)
    x_ind = sqrtCovX %*% x_ind
    x = rbind(x, x_ind)
  }
  epsilon = rnorm(N)
  for (i in 1:I) {
    CovY = matrix(0,n,n)
    for (j1 in 1:4) {
      for (j2 in 1:4) {
        if (j1==j2) CovY[j1,j2] = (1/4+1/(1+exp(4*x[(i-1)*n+j1]))) * (1/4+1/(1+exp(4*x[(i-1)*n+j2])))
        if (j1!=j2) CovY[j1,j2] = 0.8 * (1/4+1/(1+exp(4*x[(i-1)*n+j1]))) * (1/4+1/(1+exp(4*x[(i-1)*n+j2])))
      }
    }
    sqrtCovY = expm::sqrtm(CovY)
    epsilon[((i-1)*n+1):(i*n)] = sqrtCovY %*% epsilon[((i-1)*n+1):(i*n)]
  }
  y = f(x) + epsilon
  rdf <- data.frame(y=y,x=x, id=rep(1:(N/n),each=n))
  return (rdf)
}
generate_test_data <- function() {
  n = 4
  I = 10000
  N = I*n
  f <- function(x) tanh(x)
  x = matrix(NA,nrow=0,ncol=1)
  for (i in 1:I) {
    x_ind = matrix(runif(n,1,2),n,1)
    x = rbind(x, x_ind)
  }
  epsilon = rnorm(N)
  CovY = toeplitz(ARMAacf(ar=c(0.3, 0.5),ma=c(-0.3),lag.max=n-1))
  sqrtCovY = expm::sqrtm(CovY)
  for (i in 1:I) epsilon[((i-1)*n+1):(i*n)] = sqrtCovY %*% epsilon[((i-1)*n+1):(i*n)]
  y = f(x) + epsilon
  rdf <- data.frame(y=y,x=x, id=rep(1:(N/n),each=n))
  return (rdf)
}

fixed_test_errors=numeric()
fixed_train_errors=numeric()
x = seq(-0.2,0.9,by=0.01)

set.seed(1)
rdf <- generate_data()
rdf_test <- generate_test_data()
rdf_train <- generate_train_data()
for (rho in x) {
  print(rho)
  ind_tree <- crf(y~x, rdf, x0=NULL, B=500, L=1, 0.9, weight_optimiser="Training MSE", correlation="equicorr", fixrho=rho)

  PREDS_ALL_00 = data.frame()
  for (bb in 1:500) PREDS_ALL_00 = rbind(PREDS_ALL_00, predict(ind_tree$forest[[1]][[bb]], rdf_test))
  fixed_test_errors = append(fixed_test_errors, mean((colSums(PREDS_ALL_00)/500-tanh(rdf_test$x))^2))

  PREDS_ALL_00 = data.frame()
  for (bb in 1:500) PREDS_ALL_00 = rbind(PREDS_ALL_00, predict(ind_tree$forest[[1]][[bb]], rdf_train))
  fixed_train_errors = append(fixed_train_errors, mean((colSums(PREDS_ALL_00)/500-tanh(rdf_train$x))^2))
}
y_test = fixed_test_errors
y_train = fixed_train_errors


# Test MSPE Plot
test_plot <- ggplot(data.frame(x=x[16:111], y=y_test[17:112]), aes(x = x, y = y)) +
  geom_line(color="red", size=0.8) +
  scale_y_continuous(
    labels = function(x) x * 10^4,
    limits = c(0,12*10^(-4))
  ) +
  labs(title = "Test MSPE",
       x = expression(rho),
       y = expression("Test MSPE " * ~ ("\u00D7"~10^-4))) +
  theme_minimal()

# Train MSPE Plot
train_plot <- ggplot(data.frame(x=x[16:111], y=y_train[17:112]), aes(x = x, y = y)) +
  geom_line(color="blue", size=0.8) +
  scale_y_continuous(
    labels = function(x) x * 10^3,
    limits = c(0,5*10^(-3))
  ) +
  labs(title = "Train MSPE",
       x = expression(rho),
       y = expression("Train MSPE " * ~ ("\u00D7"~10^-3))) +
  theme_minimal()

gridExtra::grid.arrange(test_plot, train_plot, ncol=2)

