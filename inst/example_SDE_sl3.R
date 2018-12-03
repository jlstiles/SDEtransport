data(data_example)
head(data_example)

covariates = list(covariates_S = c("W1","W2"),
                  covariates_A = c("S","W1", "W2"),
                  covariates_Z = c("S","A","W1","W2"),
                  covariates_M = c("Z","W1","W2"),
                  covariates_Y = c("M","Z","W1","W2"),
                  covariates_QZ = c("S","W1","W2"))

# make the superlearner using sl3
lglm = make_learner(Lrnr_glm)
lmean = make_learner(Lrnr_mean)
# lxgboost = make_learner(Lrnr_xgboost)
lrnr_stack = make_learner(Stack, list(lglm, lmean))
metalearner = make_learner(Lrnr_nnls)

# define the superlearner
sl <- Lrnr_sl$new(learners = lrnr_stack,
                  metalearner = metalearner)

# takes about 9 seconds
time = proc.time()
result = SDE_tmle(data_example, sl = sl, covariates = covariates, 
                   truncate = list(lower =.0001, upper = .9999), alpha = .05, transport = TRUE)
proc.time() - time 

result$CI
head(result$IC)
result$gcomp_ests
result$fluctuation_eps

# If not transporting:
covariatesNT = list(
                  covariates_A = c("W1", "W2"),
                  covariates_Z = c("A","W1","W2"),
                  covariates_M = c("Z","W1","W2"),
                  covariates_Y = c("M","Z","W1","W2"),
                  covariates_QZ = c("W1","W2"))
# debug(SDE_tmleNT)
resultNT = SDE_tmle(data_example, sl = sl, covariates = covariatesNT, 
                    truncate = list(lower =.0001, upper = .9999), alpha = .05, transport = FALSE)

resultNT$CI
head(resultNT$IC)
resultNT$gcomp_ests
resultNT$fluctuation_eps

