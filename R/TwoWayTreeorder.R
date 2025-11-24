#' @title Two-Way Likelihood Ratio Test under Tree Order Restriction
#' @description
#' Performs the likelihood ratio test (LRT) for testing
#' \eqn{H0: \alpha_0 = \alpha_1 = ... = \alpha_a} versus the
#' tree-ordered alternative \eqn{H1: \alpha_0 \le \alpha_i} for all i,
#' with at least one strict inequality.
#'
#' @param sample_data A list-of-lists: sample_data[[i]][[j]] = numeric vector for cell (i,j).
#' @param significance_level Significance level (default = 0.05).
#' @importFrom stats quantile rnorm var
#' @import Matrix
#' @import MASS
#' @export
#' @export
#' @export
TwoWayTreeLRT <- function(sample_data, significance_level){
  set.seed(123)
  a <- length(sample_data)
  b <- length(sample_data[[1]])

  R_MLE <- function(X, n) {
    X1 <- X[-1]
    n1 <- n[-1]
    sorted_indices <- order(X1)
    X1_sorted <- X1[sorted_indices]
    n1_sorted <- n1[sorted_indices]
    if (length(X)==2){
      if (all(X1 < X[1])){
        new_X <- c((X[1]*n[1]+X[2]*n[2])/sum(n),(X[1]*n[1]+X[2]*n[2])/sum(n))
      }
      else{
        new_X <- X
      }
      return(new_X)
    }
    else {
      A <- numeric(length(X1_sorted))
      for (j in 2:length(X)) {
        A[j-1] <- (n[1] * X[1] + sum(n1_sorted[1:(j - 1)] * X1_sorted[1:(j - 1)])) /
          (n[1] + sum(n1_sorted[1:(j - 1)]))
      }
      if (all(X1 >= X[1])) {
        new_X <- X
      } else if (A[length(X)-2] >= X1_sorted[length(X)-1]) {
        X <- rep(A[length(X)-1], length(X))
        new_X <- X
      } else {
        comparisons <- logical(length(X1_sorted) - 1)
        comparisons1 <- logical(length(X1_sorted) - 1)
        stoblack_values <- numeric(0)
        for (k in 1:(length(X1_sorted) - 1)) {
          comparisons[k] <- A[k] < X1_sorted[k + 1]
          if(comparisons1[k] <- A[k] < X1_sorted[k + 1]) {
            for (s in 1:k) {
              stoblack_values[s] <- X1_sorted[s]
            }
            break
          }
        }
        selected_A_values <- A[comparisons]
        X[1] <- selected_A_values[1]
        for (l in 2:length(X)) {
          if (X[l] %in% stoblack_values) {
            X[l] <- selected_A_values[1]
          }
        }
        new_X <- X
      }
      return(new_X)
    }
  }

  Func_MLE_0<-function(sample_data_list){
    cell_means <- sapply(1:a, function(i)
      sapply(1:b, function(j) mean(sample_data_list[[i]][[j]]))
    )
    cell_means<-matrix(cell_means, nrow = a, ncol = b, byrow = TRUE)
    cell_ns <- sapply(1:a, function(i)
      sapply(1:b, function(j) length(sample_data_list[[i]][[j]]))
    )
    cell_ns <- matrix(cell_ns, nrow = a, ncol = b, byrow = TRUE)
    cell_vars <- sapply(1:a, function(i)
      sapply(1:b, function(j) ((cell_ns[i,j] -1)/cell_ns[i,j])*var(sample_data_list[[i]][[j]]))
    )
    cell_vars <- matrix(cell_vars, nrow = a, ncol = b, byrow = TRUE)
    beta_0_0<-colSums(cell_means)
    var_0_0<-cell_vars
    repeat{
      beta_1_0 <- sapply(1:b, function(j) {
        sum_1 <- sum(cell_ns[,j] / var_0_0[,j])
        sum_2 <- sum((cell_ns[,j] / var_0_0[,j]) * cell_means[,j])
        sum_2 / sum_1
      })
      var_1_0 <- sapply(1:a, function(i)
        sapply(1:b, function(j) (sum((sample_data_list[[i]][[j]]-beta_1_0[j])^2))/cell_ns[i,j]
        ))
      # ensure var_1_0 is matrix a x b
      var_1_0 <- matrix(var_1_0, nrow = a, ncol = b, byrow = TRUE)

      if(max(abs(beta_1_0-beta_0_0))<0.00001)
      {
        break
      }
      beta_0_0<-beta_1_0
      var_0_0 <- var_1_0
    }
    return(var_1_0)
  }

  Func_MLE_1<-function(sample_data_list){
    cell_means <- sapply(1:a, function(i)
      sapply(1:b, function(j) mean(sample_data_list[[i]][[j]]))
    )
    cell_means<-matrix(cell_means, nrow = a, ncol = b, byrow = TRUE)
    cell_ns <- sapply(1:a, function(i)
      sapply(1:b, function(j) length(sample_data_list[[i]][[j]]))
    )
    cell_ns <- matrix(cell_ns, nrow = a, ncol = b, byrow = TRUE)
    cell_vars <- sapply(1:a, function(i)
      sapply(1:b, function(j) ((cell_ns[i,j] -1)/cell_ns[i,j])*var(sample_data_list[[i]][[j]]))
    )
    cell_vars <- matrix(cell_vars, nrow = a, ncol = b, byrow = TRUE)
    alpha_0_1 <- rowSums(cell_means)
    beta_0_1 <- colSums(cell_means)-mean(cell_means)
    var_0_1 <- cell_vars
    u_0_1 <- cell_ns/cell_vars

    repeat{
      w <- rowSums(u_0_1)
      x <- sapply(1:a, function(i) {
        sum_4 <- sum(u_0_1[i, ] * cell_means[i, ])
        sum_5 <- sum((u_0_1[i, b]-u_0_1[i, 1:(b-1)]) * beta_0_1[1:(b-1)])
        (1 / w[i]) * (sum_4 + sum_5)
      })
      alpha_1<- R_MLE(x,w)  #### Updated value of alpha

      # Compute t_b
      t_b <- sapply(1:(b-1), function(j) {
        sum(u_0_1[, j] * (cell_means[, j] - alpha_1) - u_0_1[, b] * (cell_means[, b] - alpha_1))
      })

      # Compute Q and solve for beta_1
      Q_1 <- a * colMeans(u_0_1[, 1:(b-1), drop = FALSE])
      Q <- diag(Q_1, nrow = b - 1) + a * mean(u_0_1[, b]) * matrix(1, nrow = b - 1, ncol = b - 1)
      beta_1_old <- solve(Q, t_b)

      # Complete beta vector
      beta_1 <- c(beta_1_old, -sum(beta_1_old))

      # Update variances
      var_1 <- sapply(1:a, function(i) {
        sapply(1:b, function(j) {
          (1 / cell_ns[i, j]) * sum((sample_data_list[[i]][[j]] - alpha_1[i] - beta_1[j])^2)
        })
      })
      # ensure var_1 is matrix a x b
      var_1 <- matrix(var_1, nrow = a, ncol = b, byrow = TRUE)

      if(max(abs(alpha_1-alpha_0_1))<0.00001 & max(abs(beta_1-beta_0_1))<0.00001)  #### stopping criteria
      {
        break
      }
      alpha_0_1<-alpha_1
      beta_0_1<-beta_1
      var_0_1 <- var_1

    }#### repeat
    return(var_1)
  }### Func_MLE_1

  Calculate_lambda <- function(sample_data_list) {
    cell_means <- sapply(1:a, function(i)
      sapply(1:b, function(j) mean(sample_data_list[[i]][[j]]))
    )
    cell_means<-matrix(cell_means,nrow = a,ncol = b,byrow = TRUE)
    cell_ns <- sapply(1:a, function(i)
      sapply(1:b, function(j) length(sample_data_list[[i]][[j]]))
    )
    cell_ns <- matrix(cell_ns, nrow = a, ncol = b, byrow = TRUE)
    cell_vars <- sapply(1:a, function(i)
      sapply(1:b, function(j) ((cell_ns[i,j] -1)/cell_ns[i,j])*var(sample_data_list[[i]][[j]]))
    )
    cell_vars <- matrix(cell_vars, nrow = a, ncol = b, byrow = TRUE)
    # Compute MLE variance estimates under H1 and H0
    mle1 <- Func_MLE_1(sample_data_list)
    mle0 <- Func_MLE_0(sample_data_list)

    # coerce to matrices with same dims and protect against zeros
    mle1 <- matrix(mle1, nrow = a, ncol = b, byrow = TRUE)
    mle0 <- matrix(mle0, nrow = a, ncol = b, byrow = TRUE)
    cell_ns <- matrix(cell_ns, nrow = a, ncol = b, byrow = TRUE)

    # Compute test statistic (likelihood ratio) elementwise then product
    elementwise <- (mle1 / mle0)^(cell_ns / 2)
    prod_val <- prod(elementwise, na.rm = TRUE)

    return(prod_val)
  }

  sample_data <- lapply(sample_data, function(row) {
    lapply(row, function(cell) cell[!is.na(cell)])
  })
  n_sample<- 10000
  n<- sapply(1:a, function(i)
    sapply(1:b, function(j) length(sample_data[[i]][[j]]))
  )
  n <- matrix(n, nrow = a, ncol = b, byrow = TRUE)
  lambda_values_star <- numeric(n_sample)
  for (i in 1:n_sample) {
    bootstrap_samples <- lapply(1:a, function(i){
      lapply(1:b, function(j){
        rnorm(n[i,j], mean = 0, sd = sqrt(var(sample_data[[i]][[j]])))
      })
    } )
    lambda_values_star[i] <- Calculate_lambda(bootstrap_samples)
  }
  sort_lambda_star <- sort(lambda_values_star)
  quantile_value <- quantile(sort_lambda_star, probs = significance_level)

  lambda <- Calculate_lambda(sample_data)
  if (lambda < quantile_value) {
    result <- "Reject null hypothesis"
  } else {
    result <- "Do not reject null hypothesis"
  }
  return(
    list(
      critical_value = quantile_value,
      test_statistic = lambda,
      decision = result
    )
  )
}



#' @title MAX Test (M1) for Two-Way Tree-Ordered Data
#' @description Performs the MAX (M1) test for two-way ANOVA with tree-ordered treatment effects.
#' @param sample_data List-of-lists of numeric vectors: sample_data[[i]][[j]] = data in cell (i,j)
#' @param significance_level Significance level (default = 0.05)
#' @return A list with critical_value, test_statistic, and decision
#' @examples
#' TreeMax(sample_data, significance_level)
#' @export

TreeMax <- function(sample_data, significance_level = 0.05) {

  ## -----------------------------------
  ## DIMENSIONS
  ## -----------------------------------
  a <- length(sample_data)
  b <- length(sample_data[[1]])

  ## Remove NA values in every cell
  sample_data <- lapply(sample_data, function(row) {
    lapply(row, function(cell) cell[!is.na(cell)])
  })

  ## -----------------------------------
  ## Max Test Function  (Recomputed inside)
  ## -----------------------------------
  Max_Func <- function(sample_data_list) {

    # sample sizes (recomputed inside)
    cell_ns <- sapply(1:a, function(i)
      sapply(1:b, function(j) length(sample_data_list[[i]][[j]]))
    )
    cell_ns<-matrix(cell_ns,nrow = a,ncol = b,byrow = TRUE)
    # cell means
    cell_mean_data <- sapply(1:a, function(i)
      sapply(1:b, function(j) mean(sample_data_list[[i]][[j]]))
    )
    cell_mean_data<-matrix(cell_mean_data,nrow = a,ncol = b,byrow = TRUE)
    # row means
    cell_mean <- sapply(1:a, function(i) mean(cell_mean_data[i, ]))

    # variance component
    var1 <- sapply(1:(a - 1), function(i) {
      sum(sapply(1:b, function(j) {
        var(sample_data_list[[i + 1]][[j]]) / cell_ns[i + 1, j] +
          var(sample_data_list[[1]][[j]])     / cell_ns[1, j]
      }))
    })

    # Test statistic
    D <- sapply(1:(a - 1), function(i) {
      (cell_mean[i + 1] - cell_mean[1]) / sqrt(var1[i] / b^2)
    })

    return(max(D))
  }

  ## OBSERVED TEST STATISTIC
  ## -----------------------------------
  obs_max <- Max_Func(sample_data)
  ## BOOTSTRAP SAMPLING
  boots_sample <- 10000

  # Precompute cell variances & sizes (based on original data)
  cell_vars <- sapply(1:a, function(i)
    sapply(1:b, function(j) var(sample_data[[i]][[j]]))
  )
  cell_vars<-matrix(cell_vars,nrow = a,ncol = b, byrow = TRUE)
  cell_ns <- sapply(1:a, function(i)
    sapply(1:b, function(j) length(sample_data[[i]][[j]]))
  )
  cell_ns<-matrix(cell_ns,nrow = a,ncol = b, byrow = TRUE)
  max_boot <- numeric(boots_sample)

  for (k in 1:boots_sample) {

    # generate bootstrap sample with same sizes & variances
    bootstrap_samples <- lapply(1:a, function(i) {
      lapply(1:b, function(j) {
        rnorm(
          n = cell_ns[i, j],
          mean = 0,
          sd = sqrt(cell_vars[i, j])
        )
      })
    })

    # compute max statistic in bootstrap
    max_boot[k] <- Max_Func(bootstrap_samples)
  }
  critical_val <- quantile(max_boot, probs = 1 - significance_level)
  decision <- ifelse(obs_max > critical_val,
                     "Reject null hypothesis",
                     "Do not reject null hypothesis")
  return(list(
    critical_value = critical_val,
    test_statistic = obs_max,
    decision = decision
  ))
}





#' @title Min Test (M2) under Tree Order Restriction
#' @description
#' Computes the minimum-type test statistic for tree-ordered treatment effects.
#' @param sample_data list-of-lists cells
#' @param significance_level numeric (default 0.05)
#' @return list(critical_value, test_statistic, decision)
#' @importFrom stats var rnorm quantile
#' @export


TreeMin <- function(sample_data, significance_level = 0.05) {

  ## -----------------------------------
  ## DIMENSIONS
  ## -----------------------------------
  a <- length(sample_data)
  b <- length(sample_data[[1]])

  ## Remove NA values in every cell
  sample_data <- lapply(sample_data, function(row) {
    lapply(row, function(cell) cell[!is.na(cell)])
  })

  ## -----------------------------------
  ## Max Test Function  (Recomputed inside)
  ## -----------------------------------
  Min_Func <- function(sample_data_list) {

    # sample sizes (recomputed inside)
    cell_ns <- sapply(1:a, function(i)
      sapply(1:b, function(j) length(sample_data_list[[i]][[j]]))
    )
    cell_ns<-matrix(cell_ns,nrow = a,ncol = b,byrow = TRUE)
    # cell means
    cell_mean_data <- sapply(1:a, function(i)
      sapply(1:b, function(j) mean(sample_data_list[[i]][[j]]))
    )
    cell_mean_data<-matrix(cell_mean_data,nrow = a,ncol = b,byrow = TRUE)
    # row means
    cell_mean <- sapply(1:a, function(i) mean(cell_mean_data[i, ]))

    # variance component
    var1 <- sapply(1:(a - 1), function(i) {
      sum(sapply(1:b, function(j) {
        var(sample_data_list[[i + 1]][[j]]) / cell_ns[i + 1, j] +
          var(sample_data_list[[1]][[j]])     / cell_ns[1, j]
      }))
    })

    # Test statistic
    D <- sapply(1:(a - 1), function(i) {
      (cell_mean[i + 1] - cell_mean[1]) / sqrt(var1[i] / b^2)
    })

    return(min(D))
  }

  ## OBSERVED TEST STATISTIC
  ## -----------------------------------
  obs_min <- Min_Func(sample_data)
  ## BOOTSTRAP SAMPLING
  boots_sample <- 10000

  # Precompute cell variances & sizes (based on original data)
  cell_vars <- sapply(1:a, function(i)
    sapply(1:b, function(j) var(sample_data[[i]][[j]]))
  )
  cell_vars<-matrix(cell_vars,nrow = a,ncol = b, byrow = TRUE)
  cell_ns <- sapply(1:a, function(i)
    sapply(1:b, function(j) length(sample_data[[i]][[j]]))
  )
  cell_ns<-matrix(cell_ns,nrow = a,ncol = b, byrow = TRUE)
  min_boot <- numeric(boots_sample)

  for (k in 1:boots_sample) {

    # generate bootstrap sample with same sizes & variances
    bootstrap_samples <- lapply(1:a, function(i) {
      lapply(1:b, function(j) {
        rnorm(
          n = cell_ns[i, j],
          mean = 0,
          sd = sqrt(cell_vars[i, j])
        )
      })
    })

    # compute max statistic in bootstrap
    min_boot[k] <- Min_Func(bootstrap_samples)
  }
  critical_val <- quantile(min_boot, probs = 1 - significance_level)
  decision <- ifelse(obs_min > critical_val,
                     "Reject null hypothesis",
                     "Do not reject null hypothesis")
  return(list(
    critical_value = critical_val,
    test_statistic = obs_min,
    decision = decision
  ))
}


#' @title Run all three tree-ordered tests
#' @description
#' Convenience wrapper that runs LRT, Min, and Max tests and returns all results.
#' @param sample_data list-of-lists cells
#' @param significance_level numeric (default 0.05)
#' @return list(LRT=..., Min=..., Max=...)
#' @export
TwoWayTreeAll <- function(sample_data, significance_level) {
  lrt <- TwoWayTreeLRT(sample_data, significance_level = significance_level)
  M1  <- TreeMax(sample_data, significance_level = significance_level)
  M2  <- TreeMin(sample_data, significance_level = significance_level)

  list(LRT = lrt, Max = M1, Min = M2)
}
