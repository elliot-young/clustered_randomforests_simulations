# Reporoducibility code for simulation in introduction

library(rpart)
library(RColorBrewer)
coul <- brewer.pal(9, "Set1")
library(ggplot2)
library(crf)

generate_data = function() {
  n = 2
  I = 10000
  N = I*n
  x = matrix(rnorm(2*N),nrow=N,ncol=2)
  epsilon = rnorm(N)
  CovY = matrix(0.8,n,n) + diag(1-0.8,n)
  sqrtCovY = expm::sqrtm(CovY)
  for (i in 1:I) epsilon[((i-1)*n+1):(i*n)] = sqrtCovY %*% epsilon[((i-1)*n+1):(i*n)]
  y = tanh(x[,1]) + tanh(x[,2]) + epsilon
  rdf <- data.frame(y=y, x=x, id=rep(1:(N/n),each=n))
  return (rdf)
}

BETAS = seq(0.7,0.95,by=0.05)
ERRORS = list()
for (betas in BETAS) ERRORS[[betas*100]] = data.frame()
update_errors <- function(ERRORS) {
  set.seed(dim(ERRORS[[100*BETAS[1]]])[1]+1)
  data = generate_data()
  for (betas in BETAS) {
    FOREST.0 <- crf(y~x.1+x.2, data, L=1, B=500, weight_optimiser="Training MSE", cp=0, beta=betas, maxdepth=30, fixrho=0)
    FOREST.crf <- crf(y~x.1+x.2, data, L=1, B=500, weight_optimiser="Pointwise variance", x0=data.frame(x.1=1,x.2=1), cp=0, beta=betas, maxdepth=30)
    datapoint <- data.frame(x.1=1,x.2=1)
    ERRORS[[100*betas]] = rbind(ERRORS[[100*betas]],
                                data.frame(ones.0=predict(FOREST.0, datapoint)$fitted,
                                           ones.crf=predict(FOREST.crf, datapoint)$fitted
                                           )
                                )
  }
  return(ERRORS)
}
for (arp in 1:500) ERRORS = update_errors(ERRORS)

vars.0 = sqbias.0 = mses.0 = numeric()
vars.crf = sqbias.crf = mses.crf = numeric()
XVALS = seq(0.6,0.95,by=0.05)
for (tags in BETAS) {
  vars.0 = append(vars.0, var(ERRORS[[100*tags]]$ones.0))
  sqbias.0 = append(sqbias.0, mean(ERRORS[[100*tags]]$ones.0-2*tanh(1))^2)
  mses.0 = append(mses.0, mean((ERRORS[[100*tags]]$ones.0-2*tanh(1))^2))
  vars.crf = append(vars.crf, var(ERRORS[[100*tags]]$ones.crf))
  sqbias.crf = append(sqbias.crf, mean(ERRORS[[100*tags]]$ones.crf-2*tanh(1))^2)
  mses.crf = append(mses.crf, mean((ERRORS[[100*tags]]$ones.crf-2*tanh(1))^2))
}

XVALS <- BETAS
FINER_XVALS <- seq(0.7,0.99,by=0.001)
data <- data.frame(
  Method = rep(c("RF", "CRF"), each = 3*length(FINER_XVALS)),
  Betas = rep(FINER_XVALS, times = 3),
  Value = c(splinefun(XVALS, mses.0, method = "natural")(FINER_XVALS),
            splinefun(XVALS, vars.0, method = "natural")(FINER_XVALS),
            splinefun(XVALS, sqbias.0, method = "natural")(FINER_XVALS),
            splinefun(XVALS, mses.crf, method = "natural")(FINER_XVALS),
            splinefun(XVALS, vars.crf, method = "natural")(FINER_XVALS),
            splinefun(XVALS, sqbias.crf, method = "natural")(FINER_XVALS)
  ),
  Metric = rep(c("Mean Squared Error", "Variance", "Squared Bias"), each = length(FINER_XVALS))
)
ggplot(data, aes(x = Betas, y = Value, color = Method, linetype = Metric)) +
  geom_line(size = 0.8) +
  labs(title = "", x = expression(beta), y = "Performance Metric") +
  scale_color_manual(values = c("RF"=coul[3], "CRF"=coul[5])) +
  theme_minimal()
