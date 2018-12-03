#' @title mediation.step1
#' @description helper function called within SDE_tmle to perform TMLE update for the first regression
#' @export
mediation.step1 = function(initdata, Y_preds, data, gstarM_astar, a, transport) {
  if (!transport) {
    H = with(data, with(initdata, ((A == a)*((M == 1)*gstarM_astar + (M == 0)*(1 - gstarM_astar))/
                                     (((M == 1)*M_ps + (M == 0)*(1 - M_ps))*(A_ps*A + (1 - A)*(1 - A_ps))))))
  } else {
    H = with(data, with(initdata, ((S == 1)*(A == a)*
                                     ((M == 1)*gstarM_astar + (M == 0)*(1 - gstarM_astar))*
                                     ((Z == 1)*ZS0_ps + (Z == 0)*(1 - ZS0_ps))*(1 - S_ps))/
                          (((M == 1)*M_ps + (M == 0)*(1 - M_ps))*
                             ((Z == 1)*Z_ps + (Z == 0)*(1 - Z_ps))*
                             (A_ps*A + (1 - A)*(1 - A_ps))*S_ps*PS0)))
  }
  
  # updates
  Qfit = try(suppressWarnings(glm(data$Y ~ 1 + offset(qlogis(Y_preds$Y_init)), family = binomial,
                                  weights = H)), silent = TRUE)
  
  if (class(Qfit)[1]=="try-error") eps = 0 else eps = Qfit$coefficients
  
  est_iptw = mean(data$Y*H)
  IC_iptw = data$Y*H - est_iptw
  return(list(Qstar_M  = plogis(qlogis(Y_preds$Y_init) + eps),
              Qstar_M1 = plogis(qlogis(Y_preds$Y_init_M1) + eps),
              Qstar_M0 = plogis(qlogis(Y_preds$Y_init_M0) + eps),
              IC_iptw = IC_iptw, est_iptw = est_iptw,
              Hm = H,
              eps = eps))
}


#' @title mediation.step2
#' @description helper function called within SDE_tmle to perform TMLE update for the 2nd regression
#' and also performs estimating equation (one step estimator)
#' @export
mediation.step2 = function(data, sl, Qstar_M, Qstar_Mg, covariates_QZ, Hm, A_ps, a, tmle = TRUE,
                           EE = FALSE, transport) {
  
  if (transport) PS0 = mean(data$S==0)
  df_QZ = data
  df_QZ$Qstar_Mg = Qstar_Mg
  
  # for tmle 
  task_QZsubAa <- sl3_Task$new(
    data = data.table::copy(df_QZ[df_QZ$A == a, ]),
    covariates = covariates_QZ,
    outcome = "Qstar_Mg",
    outcome_type = "continuous"
  )
  
  task_data <- sl3_Task$new(
    data = data.table::copy(data),
    covariates = covariates_QZ,
    outcome = "Y"
  )
  
  QZfit = sl$train(task_QZsubAa)
  
  
  # compute the clever covariate and update if tmle
  if (tmle) {
    if (!transport) {
      Hz = with(data, (A == a)/(A*A_ps + (1 - A)*(1 - A_ps)))
    } else {
      Hz = with(data, (A == a)*(S == 0)/(A*A_ps + (1 - A)*(1 - A_ps))/PS0)  
    }
    QZ_preds_a = pmin(pmax(QZfit$predict(task_data), .001), .999)
    # update
    QZfit_tmle = try(suppressWarnings(glm(Qstar_Mg ~ 1 + offset(qlogis(QZ_preds_a)), family = binomial,
                                          weights = Hz)), silent = TRUE)
    if (class(QZfit_tmle)[1]=="try-error") eps2 = 0 else eps2 = QZfit_tmle$coefficients
    
    QZstar_a = plogis(qlogis(QZ_preds_a) + eps2)
    if (transport) {
      est = mean(QZstar_a[data$S==0])
    } else {
      est = mean(QZstar_a)
    }
    
    D_Y = with(data, Hm*(Y - Qstar_M))
    D_Z = Hz*(Qstar_Mg - QZstar_a)
    if (transport) {
    D_W = with(data, (QZstar_a - est)*(S ==0)/PS0)
    } else {
      D_W = with(data, (QZstar_a - est))
    }
    D = D_Y + D_Z + D_W
    
  } 
  
  # regress if EE or mle, EE updates the estimate, mle does not
  if (EE) {
    QZstar_a = pmin(pmax(QZfit$predict(task_data), .001), .999) 
    
    if (transport) {
      init_est = mean(QZstar_a[data$S==0])
    } else
    {
      init_est = mean(QZstar_a)
    }
    D_Y1s = with(data, Hm*(Y - Qstar_M))
    if (!transport) {
      Hz = with(data, (A == a)/(A*A_ps + (1 - A)*(1 - A_ps)))
    } else {
      Hz = with(data, (A == a)*(S == 0)/(A*A_ps + (1 - A)*(1 - A_ps))/PS0)  
    }
    D_Z1s = Hz*(Qstar_Mg - QZstar_a)
    if (transport) { 
      D_W1s = with(data, (QZstar_a - init_est)*(S ==0)/PS0)
    } else {
      D_W1s = with(data, (QZstar_a - init_est))
    }
    D = D_Y1s + D_Z1s + D_W1s
    # update the estimate
    est = init_est + mean(D)
  }
  
  if (tmle) return(list(est = est, IC = D, eps2 = eps2)) else {
    return(list(est = est, IC = D, init_est = init_est))
  }
  
}

