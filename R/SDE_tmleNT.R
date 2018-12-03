#' @title SDE_tmleNT
#' @description Called within SDE_tmle to compute stochastic direct and indirect effects when no
#' transporting is desired
#' @export
SDE_tmleNT = function(data, sl, V=10, covariates, truncate = list(lower =.0001, upper = .9999), 
                    alpha = .05) 
{

  # n = 1e3
  # data = gendata.SDEtransport(n, f_W = f_W, f_S = f_S, f_A = f_A, f_Z = f_Z, f_M = f_M, f_Y = f_Y)
  # 
  
  L = truncate$lower
  U = truncate$upper
  
  # get the stochastic dist of M and true params if you want 
  gstar_info = get_gstarM(data = data, sl, V=10, covariates = covariates, transport = FALSE)
  gstarM_astar = list(gstarM_astar0 = gstar_info$gstarM_astar0, gstarM_astar1 = gstar_info$gstarM_astar1)
  
  # perform initial fits for the first regression
  init_info = get.mediation.initdataNT(data = data, covariates = covariates, sl = sl)
  
  Y_preds = init_info$Y_preds
  est_info = lapply(0:1, FUN = function(astar) {
    lapply(0:1, FUN = function(a) {
      # get tmle info
      # get iptw here while I'm at it
      update = mediation.step1(initdata = init_info$initdata, init_info$Y_preds, data = data, 
                               gstarM_astar[[astar+1]], a, transport = FALSE)
      iptw_info = list(update$IC_iptw, update$est_iptw)
      
      Y_Mg = get.stochasticM(gstarM_astar[[astar+1]], Y_preds[[2]], Y_preds[[3]]) 
      A_ps = init_info$initdata$A_ps
      EE_gcomp_info = mediation.step2(data = data, sl = sl, Qstar_M = Y_preds[[1]], 
                                      Qstar_Mg = Y_Mg, covariates$covariates_QZ, Hm = update$Hm, A_ps = A_ps, 
                                      a = a, tmle = FALSE,
                                      EE = TRUE, transport = FALSE)
      
      Qstar_Mg = get.stochasticM(gstarM_astar[[astar+1]], update$Qstar_M1, update$Qstar_M0) 
      # compute Qstar_Mg here
      
      tmle_info = mediation.step2(data = data, sl = sl, Qstar_M = update$Qstar_M, 
                                  Qstar_Mg = Qstar_Mg, covariates$covariates_QZ, Hm = update$Hm, A_ps = A_ps, 
                                  a = a, tmle = TRUE,
                                  EE = FALSE, transport = FALSE)
      
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



