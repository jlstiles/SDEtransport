#' @title SDE_tmle
#' @description computes the sequential regression, targeted maximum likelihood estimate
#' for the stochastic direct effect and stochastic indirect effect when the outcome and 
#' mediator model are only available on site 1 (S = 1).  This is a data adaptive parameter
#' as the stochastic direct effect has a model for the mediator is determined on the data for 
#' site 1. We note, we are working on future generalizations to allow the mediator to be present
#' for both sites and also to allow specification of separate superlearners for all models that are 
#' fitted.  
#' @param data, data.frame of variables in time ordering from left to right. Make sure if using a 
#' continuous outcome that it is scaled to be between 0 and 1 via (Y - minY)/(maxY - minY).  The mediator
#' M, intermediate confounder Z, treatment A,and site S are all binaries and must be named as such.  
#' @param a, the treatment intervention of interest, either 0 or 1
#' @param a_star, the treatment intervention of the stochastic model for Mediator, M, either 0 or 1.
#' @param sl the sl3 superlearner defined, see sl3 documentation for defining a superlearner
#' and the example below
#' @param V number of folds for cross-validation (fixed to 10 for now)
#' @param covariates, dataframe of covariates for each necessary model, going backwards from the 
#' outcome. 
#' @param truncate, a list with elements lower and upper to truncate the various p-scores  
#' not functional at present
#' @param alpha = type I error rate
#' @param glm_only = fit using main terms glm
#' @return  a list with the following elements
#' @param CI confidence intervals for tmle, one step estimator, iptw estimators for stochastic
#' direct and indirect effects
#' @param IC dataframe of influence curve approximations for tmle, the one step estimator, 
#' and iptw estimator for both stochastic direct and indirect effects
#' @param gcomp_ests stochastic direct and indirect effects for gcomp, IC inference is 
#' generally not valid here so no confidence interval offered.
#' @param fluctuation_eps The tmle here uses the clever covariate as a weight in an intercept
#' logistic regression model. This returns the intercepts for the three parameters defined as 
#' means, respectively, under the stochastic interventions astar = a = 0, astar = a = 1 and
#' astar = 0 with a = 1.
#' @example /inst/example_SDE_sl3.R 
#' @export
SDE_tmle = function(data, sl, V=10, covariates, truncate = list(lower =.0001, upper = .9999), 
                    alpha = .05, glm_only = FALSE) 
{
  
  # n = 1e3
  # data = gendata.SDEtransport(n, f_W = f_W, f_S = f_S, f_A = f_A, f_Z = f_Z, f_M = f_M, f_Y = f_Y)
  # 
  if (glm_only) sl = make_learner(Lrnr_glm, family = binomial())
  
  L = truncate$lower
  U = truncate$upper
  
  # get the stochastic dist of M and true params if you want 
  gstar_info = get_gstarM(data = data, sl, V=10, covariates = covariates)
  gstarM_astar = list(gstarM_astar0 = gstar_info$gstarM_astar0, gstarM_astar1 = gstar_info$gstarM_astar1)
  
  # perform initial fits for the first regression
  init_info = get.mediation.initdata(data = data, covariates = covariates, sl = sl)
  
  Y_preds = init_info$Y_preds
  est_info = lapply(0:1, FUN = function(astar) {
    lapply(0:1, FUN = function(a) {
      # get tmle info
      # get iptw here while I'm at it
      update = mediation.step1(initdata = init_info$initdata, init_info$Y_preds, data = data, 
                               gstarM_astar[[astar+1]], a)
      iptw_info = list(update$IC_iptw, update$est_iptw)
      
      Y_Mg = get.stochasticM(gstarM_astar[[astar+1]], Y_preds[[2]], Y_preds[[3]]) 
      A_ps = init_info$initdata$A_ps
      EE_gcomp_info = mediation.step2(data = data, sl = sl, Qstar_M = Y_preds[[1]], 
                                      Qstar_Mg = Y_Mg, covariates$covariates_QZ, Hm = update$Hm, A_ps = A_ps, 
                                      a = a, tmle = FALSE,
                                      EE = TRUE)
      
      Qstar_Mg = get.stochasticM(gstarM_astar[[astar+1]], update$Qstar_M1, update$Qstar_M0) 
      # compute Qstar_Mg here
      
      tmle_info = mediation.step2(data = data, sl = sl, Qstar_M = update$Qstar_M, 
                                  Qstar_Mg = Qstar_Mg, covariates$covariates_QZ, Hm = update$Hm, A_ps = A_ps, 
                                  a = a, tmle = TRUE,
                                  EE = FALSE)
      
      # compile all estimates
      tmle_info$eps1 = update$eps
      return(list(EE_gcomp_info = EE_gcomp_info, iptw_info = iptw_info, tmle_info = tmle_info))
    }) 
  })
  
  # run the bootstrap 500 times
  # bootstrap from the data and thus the gstarM_astar1 and gstarM_astar0
  n = nrow(data)
  
  ests_astar0a1 =   c(tmle = est_info[[1]][[2]]$tmle_info$est,
                      EE = est_info[[1]][[2]]$EE_gcomp_info$est,
                      iptw = est_info[[1]][[2]]$iptw_info[[2]],
                      gcomp = est_info[[1]][[2]]$EE_gcomp_info$init_est)
  
  ests_astar0a0 =   c(tmle = est_info[[1]][[1]]$tmle_info$est,
                      EE = est_info[[1]][[1]]$EE_gcomp_info$est,
                      iptw = est_info[[1]][[1]]$iptw_info[[2]],
                      gcomp = est_info[[1]][[1]]$EE_gcomp_info$init_est)
  
  ests_astar1a1 =   c(tmle = est_info[[2]][[2]]$tmle_info$est,
                      EE = est_info[[2]][[2]]$EE_gcomp_info$est,
                      iptw = est_info[[2]][[2]]$iptw_info[[2]],
                      gcomp = est_info[[2]][[2]]$EE_gcomp_info$init_est)
  
  D_SDE = est_info[[1]][[2]]$tmle_info$IC - est_info[[1]][[1]]$tmle_info$IC
  D_SIE = est_info[[2]][[2]]$tmle_info$IC- est_info[[1]][[2]]$tmle_info$IC
  
  D_SDE_1s = est_info[[1]][[2]]$EE_gcomp_info$IC - est_info[[1]][[1]]$EE_gcomp_info$IC
  D_SIE_1s = est_info[[2]][[2]]$EE_gcomp_info$IC- est_info[[1]][[2]]$EE_gcomp_info$IC
  
  D_SDE_iptw = est_info[[1]][[2]]$iptw_info[[1]] - est_info[[1]][[1]]$iptw_info[[1]]
  D_SIE_iptw = est_info[[2]][[2]]$iptw_info[[1]]- est_info[[1]][[2]]$iptw_info[[1]]
  
  SE_SDE = sd(D_SDE)/sqrt(n)
  SE_SIE = sd(D_SIE)/sqrt(n)
  
  SE_SDE_1s = sd(D_SDE_1s)/sqrt(n)
  SE_SIE_1s = sd(D_SIE_1s)/sqrt(n)
  
  SE_SDE_iptw = sd(D_SDE_iptw)/sqrt(n)
  SE_SIE_iptw = sd(D_SIE_iptw)/sqrt(n)
  
  SDE_ests = ests_astar0a1 - ests_astar0a0
  SIE_ests = ests_astar1a1 - ests_astar0a1
  
  z_alpha = qnorm(1-alpha/2,0,2)
  CI_tmle = rbind(c(SDE_ests[1], SE_SDE, SDE_ests[1] - z_alpha*SE_SDE, SDE_ests[1] + z_alpha*SE_SDE),
                  c(SIE_ests[1], SE_SIE, SIE_ests[1] - z_alpha*SE_SIE, SIE_ests[1] + z_alpha*SE_SIE))
  
  CI_1s = rbind(c(SDE_ests[2], SE_SDE_1s, SDE_ests[2] - z_alpha*SE_SDE_1s, 
                      SDE_ests[2] + z_alpha*SE_SDE_1s),c(SIE_ests[2], SE_SIE_1s, 
                                                         SIE_ests[2] - z_alpha*SE_SIE_1s, 
                                                         SIE_ests[2] + z_alpha*SE_SIE_1s))
  
  CI_iptw = rbind(c(SDE_ests[3], SE_SDE_iptw, SDE_ests[3] - z_alpha*SE_SDE_iptw, 
                    SDE_ests[3] + z_alpha*SE_SDE_iptw), c(SIE_ests[3], SE_SIE_iptw, 
                                                          SIE_ests[3] - z_alpha*SE_SIE_iptw, 
                                                          SIE_ests[3] + z_alpha*SE_SIE_iptw))
  
  # naming CI slots
  colnames(CI_tmle) = colnames(CI_1s) = colnames(CI_iptw) = 
    c("est", "se", "left", "right")
  
  fluctuation_eps = rbind(c(eps1_astar0a1 = est_info[[1]][[2]]$tmle_info$eps1,
                      eps1_astar0a0 = est_info[[1]][[1]]$tmle_info$eps1,
                      eps1_astar1a1 = est_info[[2]][[2]]$tmle_info$eps1),
                      c(eps2_astar0a1 = est_info[[1]][[2]]$tmle_info$eps2,
                         eps2_astar0a0 = est_info[[1]][[1]]$tmle_info$eps2,
                         eps2_astar1a1 = est_info[[2]][[2]]$tmle_info$eps2))
  rownames(fluctuation_eps) = c("1st regression", "2nd regression")
  colnames(fluctuation_eps) = c("astar=0,a=1", "astar=a=0", "astar=a=1")
  
  gcomp_ests = c(SDE_gcomp = SDE_ests[4], SIE_gcomp = SIE_ests[4])
  names(gcomp_ests) = rownames(CI_tmle) = rownames(CI_1s) = rownames(CI_iptw) = 
    c("SDE", "SIE")
  
  CI = list(CI_tmle = CI_tmle, CI_1step = CI_1s, CI_iptw = CI_iptw)
  
  IC = data.frame(IC_SDE_tmle = D_SDE, IC_SIE_tmle = D_SIE,
                  IC_SDE_1step = D_SDE_1s, IC_SIE_1step = D_SIE_1s,
                  IC_SDE_iptw = D_SDE_iptw, IC_SIE_iptw = D_SIE_iptw)
  
  return(list(CI = CI, IC = IC, gcomp_ests = gcomp_ests, fluctuation_eps = fluctuation_eps))
  
} 



