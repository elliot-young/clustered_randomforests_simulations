# Reproducibility code for Simulation 2

library(crf)

# Data Generating Mechanism
generate_data <- function(n,I,d) {
  N = I*n
  X = matrix(rnorm(N*d), nrow=N, ncol=d)
  epsilon = 1*rnorm(N)
  y = rep(NA,N)
  CovY = toeplitz(ARMAacf(ar=c(0.6,0.3), lag.max=n-1))
  sqrtCovY = expm::sqrtm(CovY)
  for (i in 1:I) epsilon[((i-1)*n+1):(i*n)] = sqrtCovY %*% epsilon[((i-1)*n+1):(i*n)]
  for (i in 1:N) y[i] = 4*sin(X[i,1]) + epsilon[i]
  rdf <- data.frame(y=y, x=X, id=rep(1:(N/n),each=n))
  return (rdf)
}

d = 1 # covariate dimension
if (d==1) {formula=y~x} else {formula=as.formula(paste("y ~", paste(paste0("x.", 1:d), collapse = " + ")))}


estimators_crf_ptwse <- ses_crf_ptwse <- numeric()
for (reps in 1:1000) {
  set.seed(reps)
  rdf <- generate_data(5,1000,d)
  OUTPUT <- crf(formula, rdf, x0=rep(1,d), B=500, L=100, beta=0.9, correlation="ar1", weight_optimiser="Pointwise variance")
  OUTPUT_analys <- predict(OUTPUT, data.frame(x=matrix(rep(1,d),nrow=1)), sderr=TRUE)
  estimators_crf_ptwse = append(estimators_crf_ptwse, OUTPUT_analys$fitted)
  ses_crf_ptwse = append(ses_crf_ptwse, OUTPUT_analys$sderr)
}

estimators_crf_unw <- ses_crf_unw <- numeric()
for (reps in 1:1000) {
  set.seed(reps)
  rdf <- generate_data(5,1000,d)
  OUTPUT <- crf(formula, rdf, x0=rep(1,d), B=500, L=100, beta=0.9, correlation="ar1", fixrho=0)
  OUTPUT_analys <- predict(OUTPUT, data.frame(x=matrix(rep(1,d),nrow=1)), sderr=TRUE)
  estimators_crf_unw = append(estimators_crf_unw, OUTPUT_analys$fitted)
  ses_crf_unw = append(ses_crf_unw, OUTPUT_analys$sderr)
}

write.csv(data.frame(estimators_unw=estimators_unw, ses_unw=ses_unw, estimators_crf=ses_crf_ptwse, ses_crf=ses_crf_ptwse),paste0("~/Downloads/output_",d,".csv"))
