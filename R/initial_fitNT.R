#' @title get.mediation.initdataNT
#' @description helper function called within SDE_tmleNT to perform initial regressions for TMLE and
#' estimating equation (one step estimator). 
#' @export
get.mediation.initdataNT = function(data, covariates, sl) {
  
  df_YM1 = data
  df_YM1$M = 1

  
  df_YM0 = df_YM1
  df_YM0$M = 0
  
  # Used for fitting and predicting Y
  task_Y <- sl3_Task$new(
    data = data.table::copy(data),
    covariates = covariates$covariates_Y,
    outcome = "Y"
  )
  
  task_YM1 <- sl3_Task$new(
    data = data.table::copy(df_YM1),
    covariates = covariates$covariates_Y,
    outcome = "Y"
  )
  
  task_YM0 <- sl3_Task$new(
    data = data.table::copy(df_YM0),
    covariates = covariates$covariates_Y,
    outcome = "Y"
  )
  
  # used for fitting M
  task_M <- sl3_Task$new(
    data = data.table::copy(data),
    covariates = covariates$covariates_M,
    outcome = "M"
  )

  task_A <- sl3_Task$new(
    data = data.table::copy(data),
    covariates = covariates$covariates_A,
    outcome = "A"
  )
  
  # Y, M, A fits 
  Afit = sl$train(task_A)

  # fitting M and Y on subsets where S is 1.  This is all that is needed for M and Y
  Mfit = sl$train(task_M)
  Yfit = sl$train(task_Y)
  
  # propensity scores
  A_ps = Afit$predict()
  M_ps = Mfit$predict(task_M)
  
  # Predict Y for whole data, also with M = 1 and 0, S does not matter since
  # Y is not a function of that variable so won't be used for prediction
  Y_init = pmin(pmax(Yfit$predict(task_Y), 0.001), .999)
  Y_init_M1 = pmin(pmax(Yfit$predict(task_YM1), 0.001), .999)
  Y_init_M0 = pmin(pmax(Yfit$predict(task_YM0), 0.001), .999)
  return(list(initdata = list(M_ps = M_ps, A_ps = A_ps), 
              Y_preds = list(Y_init = Y_init, 
                             Y_init_M1 = Y_init_M1, 
                             Y_init_M0 = Y_init_M0)))
}

