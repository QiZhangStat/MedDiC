library(Rcpp)
library(RcppArmadillo)

sourceCpp('meddic.cpp')
source('meddic.R')

# We demonstrate the use with the following simulation example
load('sim_data.rdata')

ls()
names(dt.sim)

Y = dt.sim$Y # the response vector
X0 = dt.sim$X0 # the confounder matrix without the intercept, can be NULL
Z = dt.sim$Z # the exposure matrix
M = dt.sim$M # the mediator matrix
dim(X0)
dim(Z)

# The whole meddic pipeline include the following four steps

# (1) prepare the data for meddic analysis
dt.meddic <- meddic_data_prep(Y,X0,Z,M,thr.lse=0.2,no.penalty.exposure=0,no.penalty.confounder=1)
# (2) calculate the scores
score.meddic <- meddic_get_score(dt.meddic)
# The score calculation for the debiased lasso is the most time consuming step, especially when p and q are large. We recommend to save the results in a file.
  save(score.meddic,file='score_sim1.rdata')
# (3)Estimating the direct, indirect and total effects
out.meddic <- meddic_get_results(dt.meddic,score.meddic)
# (4) Perform inference: calculating the raw p-values without adjusting multiple testing and 95% CIs.
inference_out <- inference.elem(out.meddic,conf.level = 0.95)

# In this simulation example, id.monly, id.donly and id.b are the IDs of the exposures with only mediation effect, only direct effect, and both effect types, respectively. All of them have non-zero total effects.
# The following evaluate the accuracy of detecting these effects at FDR=0.05
# mediation effect
which(p.adjust(inference_out$pval.asymp$indirect,method='fdr')<0.05)
c(dt.sim$id.b,dt.sim$id.monly)
# direct effect
which(p.adjust(inference_out$pval.asymp$direct,method='fdr')<0.05)
c(dt.sim$id.b,dt.sim$id.donly)
# total effect
which(p.adjust(inference_out$pval.asymp$tot,method='fdr')<0.05)
c(dt.sim$id.b,dt.sim$id.monly, dt.sim$id.donly)


# meddic_Wrapper is a wrapper function that run the whole procedure. NOT RUN
#meddic_out <- meddic_Wrapper(Y,X0,Z,M,thr.lse=0.2,no.penalty.exposure=0,no.penalty.confounder=1,f_score = 'score_sim.rdata')
#inference_out <- meddic_out$results_inference_alpha05


