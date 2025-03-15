
library(Rcpp)

cppFunction('List XWX_XWY_equicorr_cpp(double rho, int num_leaves, int I, NumericVector n_i, NumericVector nodesis, NumericVector Y) {
  double theta = rho/(1-rho);
  NumericVector nodesis_shift = nodesis-1;
  NumericMatrix XWX(num_leaves, num_leaves);
  NumericVector XWY(num_leaves);

  int index_tracker1 = 0;
  int index_tracker2 = 0;
  int node_gg1;
  int node_gg2;

  for(int i = 0; i < I; i++) {
    int n = n_i(i);
    double phi = theta/(1+n*theta);
    NumericVector nodes_i(n);
    NumericVector y_i(n);

    index_tracker2 = index_tracker1+n-1;
    for (int gg = 0; gg < n; gg++) {
      nodes_i(gg) = nodesis_shift(index_tracker1 + gg);
      y_i(gg) = Y(index_tracker1 + gg);
    }
    index_tracker1 = index_tracker2+1;
    double phi_sum_y_i = phi * sum(y_i);

    // Comp. improvements over commented code (for inv equicorr structure)
    for (int gg1 = 0; gg1 < n; gg1++) {
      node_gg1 = nodes_i(gg1);
      XWX(node_gg1, node_gg1) += 1;
      XWY(node_gg1) += y_i(gg1) - phi_sum_y_i;
      for (int gg2 = 0; gg2 < n; gg2++) {
        node_gg2 = nodes_i(gg2);
        XWX(node_gg1, node_gg2) -= phi;
      }
    }

  }

  List list_XWX_XWY = List::create(XWX, XWY);

  return list_XWX_XWY;
}
')

cppFunction('List XWX_XWSWX_XX_equicorr_cpp(double rho, int num_leaves, int I, NumericVector n_i, NumericVector nodesis, NumericVector epsilon) {
  double theta = rho/(1-rho);
  NumericVector nodesis_shift = nodesis-1;
  NumericMatrix XWSWX(num_leaves, num_leaves);
  NumericMatrix XWX(num_leaves, num_leaves);
  NumericMatrix XX(num_leaves, num_leaves);

  int index_tracker1 = 0;
  int index_tracker2 = 0;
  int node_gg1;
  int node_gg2;

  for(int i = 0; i < I; i++) {
    int n = n_i(i);
    double phi = theta/(1+n*theta);
    NumericVector nodes_i(n);
    NumericVector tilde_blank_i(n);
    NumericVector epsilon_i(n);
    NumericVector epsilon_pert_i(n);

    index_tracker2 = index_tracker1+n-1;
    // Extract node numbers and epsilon vector of a single group (index i)
    // NTS: Could do the below in one single loop if updating epsilon_i at
    //      gg+1 extra along each line (this format an easier baseline for
    //      adapting to other weights though).
    for (int gg = 0; gg < n; gg++) {
      nodes_i(gg) = nodesis_shift(index_tracker1 + gg);
      epsilon_i(gg) = epsilon(index_tracker1 + gg);
    }
    double phi_sum_eps_i = phi * sum(epsilon_i);
    for (int gg = 0; gg < n; gg++) {
      epsilon_pert_i(gg) = epsilon_i(gg) - phi_sum_eps_i;
    }
    index_tracker1 = index_tracker2+1;


    for (int gg1 = 0; gg1 < n; gg1++) {
      node_gg1 = nodes_i(gg1);
      XWX(node_gg1, node_gg1) += 1;
      XX(node_gg1, node_gg1) += 1;
      for (int gg2 = 0; gg2 < n; gg2++) {
        node_gg2 = nodes_i(gg2);
        XWX(node_gg1, node_gg2) -= phi;
        XWSWX(node_gg1, node_gg2) += epsilon_pert_i(gg1) * epsilon_pert_i(gg2);
      }
    }

  }

  List list_XWX_XWSWX_XX = List::create(XWX, XWSWX, XX);

  return list_XWX_XWSWX_XX;
}
')

cppFunction('List XWX_XWSWX_XtestX_equicorr_cpp(double rho, int num_leaves, int I, NumericVector n_i, NumericVector nodesis, NumericVector epsilon, int N_test, NumericVector nodesis_test) {
  double theta = rho/(1-rho);
  NumericVector nodesis_shift = nodesis-1;
  NumericVector nodesis_shift_test = nodesis_test-1;
  NumericMatrix XWSWX(num_leaves, num_leaves);
  NumericMatrix XWX(num_leaves, num_leaves);
  NumericMatrix XX(num_leaves, num_leaves);

  int index_tracker1 = 0;
  int index_tracker2 = 0;
  int node_gg1;
  int node_gg2;

  for(int i = 0; i < I; i++) {
    int n = n_i(i);
    double phi = theta/(1+n*theta);
    NumericVector nodes_i(n);
    NumericVector tilde_blank_i(n);
    NumericVector epsilon_i(n);
    NumericVector epsilon_pert_i(n);

    index_tracker2 = index_tracker1+n-1;
    // Extract node numbers and epsilon vector of a single group (index i)
    // NTS: Could do the below in one single loop if updating epsilon_i at
    //      gg+1 extra along each line (this format an easier baseline for
    //      adapting to other weights though).
    for (int gg = 0; gg < n; gg++) {
      nodes_i(gg) = nodesis_shift(index_tracker1 + gg);
      epsilon_i(gg) = epsilon(index_tracker1 + gg);
    }
    double phi_sum_eps_i = phi * sum(epsilon_i);
    for (int gg = 0; gg < n; gg++) {
      epsilon_pert_i(gg) = epsilon_i(gg) - phi_sum_eps_i;
    }
    index_tracker1 = index_tracker2+1;


    for (int gg1 = 0; gg1 < n; gg1++) {
      node_gg1 = nodes_i(gg1);
      XWX(node_gg1, node_gg1) += 1;
      for (int gg2 = 0; gg2 < n; gg2++) {
        node_gg2 = nodes_i(gg2);
        XWX(node_gg1, node_gg2) -= phi;
        XWSWX(node_gg1, node_gg2) += epsilon_pert_i(gg1) * epsilon_pert_i(gg2);
      }
    }

  }

  for(int i = 0; i < N_test; i++) {
    int loc = nodesis_shift_test(i);
    XX(loc,loc) += 1;
  }

  List list_XWX_XWSWX_XX = List::create(XWX, XWSWX, XX);

  return list_XWX_XWSWX_XX;
}
')

optimise_randeff_trainmse_equicorr <- function(num_leaves, I, nis, nodesis, epsilon, ...) {
  LOSS <- function(rho) {
    XWX_XWSWX_mats <- XWX_XWSWX_XX_equicorr_cpp(rho, num_leaves, I, nis, nodesis, epsilon)
    solveXWX <- solve(XWX_XWSWX_mats[[1]]+diag(1e-6,num_leaves,num_leaves))
    loss <- sum(diag(  XWX_XWSWX_mats[[3]] %*% solveXWX %*% XWX_XWSWX_mats[[2]] %*% solveXWX  ))
    return(loss)
  }
  rho_optim <- optim(0, LOSS, method="Brent", lower=-0.9, upper=0.9)$par
}

optimise_randeff_testmse_equicorr <- function(num_leaves, I, nis, nodesis, epsilon, test_data, test_density, tree, convert_pred_to_leaf, ...) {
  if (is.null(test_data) && !is.null(test_density)) {
    error("test_density not currently supported. Create your own dataset following the chosen test_density and input it in test_data")
  }
  N.TEST <- nrow(test_data)
  preds_test <- predict(tree, test_data)
  nodesisTEST <- vapply(preds_test, convert_pred_to_leaf, numeric(1))
  LOSS <- function(rho) {
    XWX_XWSWX_mats <- XWX_XWSWX_XtestX_equicorr_cpp(rho, num_leaves, I, nis, nodesis, epsilon, N.TEST, nodesisTEST)
    solveXWX <- solve(XWX_XWSWX_mats[[1]]+diag(1e-6,num_leaves,num_leaves))
    loss <- sum(diag(  XWX_XWSWX_mats[[3]] %*% solveXWX %*% XWX_XWSWX_mats[[2]] %*% solveXWX  ))
    return(loss/N.TEST)
  }
  rho_optim <- optim(0, LOSS, method="Brent", lower=-0.9, upper=0.9)$par

  return(rho_optim)
}


