# Reproducible Code for Simulation 1

# Load relevant functions
source("sim1.functions.R")
# ^Note we load individual functions from crf for computational ease - we can run the common practices between all of RF, TRAIN, CRF, and REEM only once.

# Train data
generate_data <- function() {
  n = 4
  I = 10000
  N = I*n
  f <- function(x) tanh(x)
  x = matrix(NA,nrow=0,ncol=1)
  CovX = diag(1,n)
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
    CovY = CovY
    sqrtCovY = expm::sqrtm(CovY)
    epsilon[((i-1)*n+1):(i*n)] = sqrtCovY %*% epsilon[((i-1)*n+1):(i*n)]
  }

  y = f(x) + epsilon

  rdf <- data.frame(y=y,x=x, id=rep(1:(N/n),each=n))
  return (rdf)
}
# Test data
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

formula = y~x
beta = 0.9 # subsampling rate
L = 1
# Additional hyperparameters
B = 500
maxdepth = 30
minbucket = 10
cp = 0

all_train_errors = all_test_errors = matrix(NA, nrow=2000, ncol=4)

for (reps in 1:2000) {

    print(paste0("Running number ",reps," of ",2000))
    set.seed(reps)
    data = generate_data()
    test_data = generate_test_data()
    train_data = generate_data()

    XWX_XWY_calc <- XWX_XWY_equicorr_cpp

    forest_train <- forest_test <- forest_nlme <- forest_unw <- list()

    full_ids <- data$id
    unique_ids <- unique(full_ids)
    I <- length(unique_ids)
    s <- ceiling(I^beta)

    for (b in 1:B) {
      set.seed(b)
      ids_selected <- sample(unique_ids, s)
      rdf <- data[full_ids %in% ids_selected,]

      unique_ids_3partition <- split(ids_selected, cut(seq_along(ids_selected), breaks = 3, labels = FALSE))
      rdf_split <- rdf[rdf$id %in% unique_ids_3partition[[1]],]
      rdf_evalfix <- rdf[rdf$id %in% unique_ids_3partition[[2]],]
      rdf_evalrand <- rdf[rdf$id %in% unique_ids_3partition[[3]],]

      # Create initial (unweighted) tree using rpart
      initial_tree <- rpart::rpart(formula, rdf_split, cp=cp, maxdepth=maxdepth, minbucket=minbucket)
      # Create eval design matrix
      leaves <- which(initial_tree$frame$var=="<leaf>")
      num_leaves <- length(leaves)
      initial_tree_where <- initial_tree$where
      inital_tree_where_1ton_indexed <- match(initial_tree_where, leaves)

      leaves_1to <- seq_len(num_leaves)
      lookup_preds <- setNames(leaves_1to, initial_tree$frame$yval[leaves])
      lookup_preds <- as.list(lookup_preds)
      convert_pred_to_leaf <- function(x) lookup_preds[[as.character(x)]]
      ####################################
      preds_evalrand <- predict(initial_tree, rdf_evalrand)
      nodes_of_evalrand <- vapply(preds_evalrand, convert_pred_to_leaf, numeric(1))
      num_evalrand <- nrow(rdf_evalrand)
      nis_evalrand <- c()
      unique_ids_rand <- unique(rdf_evalrand$id)
      for (i in 1:length(unique_ids_rand)) nis_evalrand <- append(nis_evalrand, sum(rdf_evalrand$id==unique_ids_rand[i]))
      I.evalrand <- length(nis_evalrand)
      epsilon_evalrand <- rdf_evalrand$y - preds_evalrand

      ###############################################################

      rho_optim_train <- optimise_randeff_trainmse_equicorr(num_leaves, I.evalrand, nis_evalrand, nodes_of_evalrand, epsilon_evalrand)
      rho_optim_test <- optimise_randeff_testmse_equicorr(num_leaves, I.evalrand, nis_evalrand, nodes_of_evalrand, epsilon_evalrand, tree=initial_tree, convert_pred_to_leaf=convert_pred_to_leaf, test_data=test_data)
      rho_unw <- 0
      NLME <- nlme::gls(yminus~1,cbind(rdf_evalrand,yminus=rdf_evalrand$y-predict(initial_tree,rdf_evalrand)), correlation=nlme:::corCompSymm(form=~1|id))
      rho_optim_nlmee <- as.numeric(coef(NLME$modelStruct, unconstrained=FALSE))

      preds_evalfix <- predict(initial_tree, rdf_evalfix)
      nodes_of_evalfix <- vapply(preds_evalfix, convert_pred_to_leaf, numeric(1))
      num_evalfix <- nrow(rdf_evalfix)

      nis_evalfix <- c()
      unique_ids_fix <- unique(rdf_evalfix$id)
      for (i in 1:length(unique_ids_fix)) nis_evalfix <- append(nis_evalfix, sum(rdf_evalfix$id==unique_ids_fix[i]))
      I.evalfix <- length(nis_evalfix)

      # TRAIN
      XWX_XWY_mats <- XWX_XWY_calc(rho_optim_train, num_leaves, I.evalfix, nis_evalfix, nodes_of_evalfix, rdf_evalfix$y)
      final_tree <- initial_tree
      final_tree$frame$yval[final_tree$frame$var=="<leaf>"] <- solve(XWX_XWY_mats[[1]]+diag(1e-6,num_leaves,num_leaves), XWX_XWY_mats[[2]])
      forest_train[[b]] <- final_tree
      # TEST
      XWX_XWY_mats <- XWX_XWY_calc(rho_optim_test, num_leaves, I.evalfix, nis_evalfix, nodes_of_evalfix, rdf_evalfix$y)
      final_tree <- initial_tree
      final_tree$frame$yval[final_tree$frame$var=="<leaf>"] <- solve(XWX_XWY_mats[[1]]+diag(1e-6,num_leaves,num_leaves), XWX_XWY_mats[[2]])
      forest_test[[b]] <- final_tree
      # UNW
      XWX_XWY_mats <- XWX_XWY_calc(rho_unw, num_leaves, I.evalfix, nis_evalfix, nodes_of_evalfix, rdf_evalfix$y)
      final_tree <- initial_tree
      final_tree$frame$yval[final_tree$frame$var=="<leaf>"] <- solve(XWX_XWY_mats[[1]]+diag(1e-6,num_leaves,num_leaves), XWX_XWY_mats[[2]])
      forest_unw[[b]] <- final_tree
      # NLME
      XWX_XWY_mats <- XWX_XWY_calc(rho_optim_nlmee, num_leaves, I.evalfix, nis_evalfix, nodes_of_evalfix, rdf_evalfix$y)
      final_tree <- initial_tree
      final_tree$frame$yval[final_tree$frame$var=="<leaf>"] <- solve(XWX_XWY_mats[[1]]+diag(1e-6,num_leaves,num_leaves), XWX_XWY_mats[[2]])
      forest_nlme[[b]] <- final_tree
    }


    dim_testdata_1 = dim(test_data)[1]
    PREDS_ALL_train = PREDS_ALL_test = PREDS_ALL_nlme = PREDS_ALL_unw = matrix(0, B, dim_testdata_1)
    for (bb in 1:B) {
      PREDS_ALL_train[bb,] = predict(forest_train[[bb]], test_data)
      PREDS_ALL_test[bb,] = predict(forest_test[[bb]], test_data)
      PREDS_ALL_nlme[bb,] = predict(forest_unw[[bb]], test_data)
      PREDS_ALL_unw[bb,] = predict(forest_nlme[[bb]], test_data)
    }

    TRAIN = mean((colSums(PREDS_ALL_train)/B-tanh(test_data$x))^2)
    TEST = mean((colSums(PREDS_ALL_test)/B-tanh(test_data$x))^2)
    NLME = mean((colSums(PREDS_ALL_nlme)/B-tanh(test_data$x))^2)
    UNW = mean((colSums(PREDS_ALL_unw)/B-tanh(test_data$x))^2)

    all_test_errors[reps,] = c(TRAIN, TEST, NLME, UNW)


    dim_testdata_1 = dim(train_data)[1]
    PREDS_ALL_train = PREDS_ALL_test = PREDS_ALL_nlme = PREDS_ALL_unw = matrix(0, B, dim_testdata_1)
    for (bb in 1:B) {
      PREDS_ALL_train[bb,] = predict(forest_train[[bb]], train_data)
      PREDS_ALL_test[bb,] = predict(forest_test[[bb]], train_data)
      PREDS_ALL_nlme[bb,] = predict(forest_unw[[bb]], train_data)
      PREDS_ALL_unw[bb,] = predict(forest_nlme[[bb]], train_data)
    }

    TRAIN = mean((colSums(PREDS_ALL_train)/B-tanh(train_data$x))^2)
    TEST = mean((colSums(PREDS_ALL_test)/B-tanh(train_data$x))^2)
    NLME = mean((colSums(PREDS_ALL_nlme)/B-tanh(train_data$x))^2)
    UNW = mean((colSums(PREDS_ALL_unw)/B-tanh(train_data$x))^2)

    all_train_errors[reps,] = c(TRAIN, TEST, NLME, UNW)

}
