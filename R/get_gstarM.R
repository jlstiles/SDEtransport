#' @export
get_gstarM  = function(data, sl, V=10, covariates) 
{
  nn = nrow(data)
  # for fitting Mstar
  task_Mstar <- sl3_Task$new(
    data = data.table::copy(data[data$S == 1,]),
    covariates = covariates$covariates_M,
    outcome = "M"
  )
  
  # for fitting Zstar
  task_Zstar <- sl3_Task$new(
    data = data.table::copy(data),
    covariates = covariates$covariates_Z,
    outcome = "Z"
  )
  
  Mstarfit = sl$train(task_Mstar)
  Zstarfit = sl$train(task_Zstar)
  
  df_ZA1 = df_ZA0 = data
  df_ZA1$A = 1
  df_ZA0$A = 0
  
  task_ZA1 <- sl3_Task$new(
    data = data.table::copy(df_ZA1),
    covariates = covariates$covariates_Z,
    outcome = "Z"
  )
  
  task_ZA0 <- sl3_Task$new(
    data = data.table::copy(df_ZA0),
    covariates = covariates$covariates_Z,
    outcome = "Z"
  )
  
  df_MZ1 = df_MZ0 = data
  df_MZ1$Z = 1
  df_MZ0$Z = 0
  
  task_MZ1 =   sl3_Task$new(
    data = data.table::copy(df_MZ1),
    covariates = covariates$covariates_M,
    outcome = "M"
  )
  
  task_MZ0 =   sl3_Task$new(
    data = data.table::copy(df_MZ0),
    covariates = covariates$covariates_M,
    outcome = "M"
  )
  
  pred_MZ1 = Mstarfit$predict(task_MZ1)
  pred_MZ0 = Mstarfit$predict(task_MZ0)
  
  # predicting Z for S = 1, Astar = 0, since we assume we only have M for S = 1
  
  df_ZA0 = df_ZA1 = data
  df_ZA0$A = 0
  df_ZA1$A = 1
  
  task_ZA0 <- sl3_Task$new(
    data = data.table::copy(df_ZA0),
    covariates = covariates$covariates_Z,
    outcome = "Z"
  )
  
  task_ZA1 <- sl3_Task$new(
    data = data.table::copy(df_ZA1),
    covariates = covariates$covariates_Z,
    outcome = "Z"
  )
  
  pred_ZA0 = Zstarfit$predict(task_ZA0)
  pred_ZA1 = Zstarfit$predict(task_ZA1)
  
  gstarM_astar0 = pred_MZ1*pred_ZA0 + pred_MZ0*(1 - pred_ZA0)
  gstarM_astar1 = pred_MZ1*pred_ZA1 + pred_MZ1*(1 - pred_ZA1)
  
  return(list(gstarM_astar1 = gstarM_astar1, gstarM_astar0 = gstarM_astar0))
}
