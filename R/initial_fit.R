#' @title get.mediation.initdata
#' @description helper function called within SDE_tmle to perform initial regressions for TMLE and
#' estimating equation (one step estimator). 
#' @export
get.mediation.initdata = function(data, covariates, sl) {
  
  df_YM1S1 = data
  df_YM1S1$M = 1
  df_YM1S1$S = 1
  
  df_YM0S1 = df_YM1S1
  df_YM0S1$M = 0
  
  df_ZS0 = data
  df_ZS0$S = 0
  
  # to be used for fitting Y
  task_YsubS1 <- sl3_Task$new(
    data = data.table::copy(data[data$S==1,]),
    covariates = covariates$covariates_Y,
    outcome = "Y"
  )
  
  # Used for predicting Y
  task_Y <- sl3_Task$new(
    data = data.table::copy(data),
    covariates = covariates$covariates_Y,
    outcome = "Y"
  )
  
  task_YM1S1 <- sl3_Task$new(
    data = data.table::copy(df_YM1S1),
    covariates = covariates$covariates_Y,
    outcome = "Y"
  )
  
  task_YM0S1 <- sl3_Task$new(
    data = data.table::copy(df_YM0S1),
    covariates = covariates$covariates_Y,
    outcome = "Y"
  )
  
  # used for fitting M
  task_MsubS1 <- sl3_Task$new(
    data = data.table::copy(data[data$S == 1,]),
    covariates = covariates$covariates_M,
    outcome = "M"
  )
  
  # used for predicting M
  task_MS1 <- sl3_Task$new(
    data = data.table::copy(df_YM1S1),
    covariates = covariates$covariates_M,
    outcome = "M"
  )
  
  
  task_Z <- sl3_Task$new(
    data = data.table::copy(data),
    covariates = covariates$covariates_Z,
    outcome = "Z"
  )
  
  task_ZS0 <- sl3_Task$new(
    data = data.table::copy(df_ZS0),
    covariates = covariates$covariates_Z,
    outcome = "Z"
  )
  
  task_A <- sl3_Task$new(
    data = data.table::copy(data),
    covariates = covariates$covariates_A,
    outcome = "A"
  )
  
  task_ZS1 <- sl3_Task$new(
    data = data.table::copy(df_YM1S1),
    covariates = covariates$covariates_Z,
    outcome = "Z"
  )
  
  task_S <- sl3_Task$new(
    data = data.table::copy(data),
    covariates = covariates$covariates_S,
    outcome = "S"
  )
  
  # Y, M, Z, A, S fits
  Sfit = sl$train(task_S)
  Afit = sl$train(task_A)
  Zfit = sl$train(task_Z)
  
  # fitting M and Y on subsets where S is 1.  This is all that is needed for M and Y
  Mfit = sl$train(task_MsubS1)
  Yfit = sl$train(task_YsubS1)
  
  # propensity scores
  S_ps = Sfit$predict()
  PS0 = mean(data$S == 0)
  A_ps = Afit$predict()
  Z_ps = Zfit$predict()
  ZS0_ps = Zfit$predict(task_ZS0)
  # might as well predict for S = 1 on whole data because only S = 1 subset is relevant
  # as clev cov is 0 otherwise 
  M_ps = Mfit$predict(task_MS1)
  # 1st clever cov FOR SDE, ALSO NEED THE ONE FOR SIE
  
  # Predict Y for whole data, also with M = 1 and 0, S does not matter since
  # Y is not a function of that variable so won't be used for prediction
  Y_init = pmin(pmax(Yfit$predict(task_Y), 0.001), .999)
  Y_init_M1 = pmin(pmax(Yfit$predict(task_YM1S1), 0.001), .999)
  Y_init_M0 = pmin(pmax(Yfit$predict(task_YM0S1), 0.001), .999)
  return(list(initdata = list(M_ps = M_ps, ZS0_ps = ZS0_ps, Z_ps = Z_ps, A_ps = A_ps, 
                              S_ps = S_ps, PS0 = PS0), 
              Y_preds = list(Y_init = Y_init, 
                             Y_init_M1 = Y_init_M1, 
                             Y_init_M0 = Y_init_M0)))
}

