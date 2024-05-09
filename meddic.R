


meddic_data_prep <- function(Y,X0,Z,M,thr.lse=0.2,no.penalty.exposure=0,no.penalty.confounder=0){
  # A wrapper function that prepare the data for Meddic
  # Inputs:
  # Y: the outcome vector
  # X0: the confounder matrix
  # Z: the matrix of exposures
  # M: the matrix of mediators
  # thr.lse: a threshold of (p+q)/n or q/n to decide when to use OLS and when to use the debiased alsso. The default is 0.2, meaning that the adaptive lasso is used if the ratio between the numbers of covariates and the sample size is above 0.2.
  # no.penalty.exposure, no.penalty.confounder: if not zero, the exposure or the confounders will not be penalized in the scaled lasso regression. For example, when the exposure is low dimensional and the confounders is high dimensional, it may be helpful to set no.penalty.exposure to be 1 and no.penalty.confounder to be 0.
  # Outputs:
  # A list with the data shaped for meddic analysis


  if(is.null(X0)){
    p0=0

  }else{
    p0 = ncol(X0)
  }
  # pre-process the data by standization and computing the inner products.
  Xtilde = stan_col_Rcpp(cbind(X0,Z))
  X = stan_col_Rcpp(cbind(X0,Z,M))
  y.stan = stan_col_Rcpp(Y)

  scale.coef = sd(Y)/apply(Z,2,sd)

  gram_Xtilde = inner_prod_XX_Rcpp(Xtilde)

  gram_X = inner_prod_XX_Rcpp(X)

  xy_tilde = inner_prod_XY_Rcpp(Xtilde,y.stan)
  xy = inner_prod_XY_Rcpp(X,y.stan)
  # determine whether to use debiased lasso or OLS
  n = length(Y)
  p.tilde = ncol(Xtilde)
  p = ncol(X)
  w_lambda_tilde = rep(1,p.tilde)
  w_lambda = rep(1,p)
  if(no.penalty.confounder[1]!=0){
    w_lambda_tilde[no.penalty.confounder]=0
    w_lambda[no.penalty.confounder]=0
  }
  if(no.penalty.exposure[1]!=0){
    w_lambda_tilde[no.penalty.exposure+p0]=0
    w_lambda[no.penalty.exposure+p0]=0
  }

  l1_penalty = ifelse(p>n*thr.lse,1,0)
  l1_penalty_tilde = ifelse(p.tilde>n*thr.lse,1,0)

  out.meddic <- list()
  out.meddic$y_stan = y.stan
  out.meddic$x_stan = X
  out.meddic$ipy = xy
  out.meddic$ipx = gram_X
  out.meddic$x_stan_tilde = Xtilde
  out.meddic$ipy_tilde = xy_tilde
  out.meddic$ipx_tilde = gram_Xtilde
  out.meddic$l1_penalty = l1_penalty
  out.meddic$l1_penalty_tilde = l1_penalty_tilde
  out.meddic$w_lambda_tilde=w_lambda_tilde
  out.meddic$w_lambda=w_lambda
  out.meddic$p0=p0
  out.meddic$scale.coef = scale.coef

  return(out.meddic)
}


meddic_get_score <- function(out.meddic,tol=1e-5, max_iter=20, max_iter_lasso=100){
  # A wrapper function that calculate the score matrices for Meddic
  # Inputs:
  # The output of meddic_data_prep
  # Outputs:
  # A list of the two scores for the marginal outcome model and the conditional outcome model
  y_stan=out.meddic$y_stan
  x_stan=out.meddic$x_stan
  ipy=out.meddic$ipy
  ipx=out.meddic$ipx
  x_stan_tilde=out.meddic$x_stan_tilde
  ipy_tilde=out.meddic$ipy_tilde
  ipx_tilde=out.meddic$ipx_tilde
  l1_penalty=out.meddic$l1_penalty
  l1_penalty_tilde=out.meddic$l1_penalty_tilde
  w_lambda_tilde=out.meddic$w_lambda_tilde
  w_lambda=out.meddic$w_lambda
  p0=out.meddic$p0

  if(l1_penalty_tilde>0){
    score_nomed = getScore_Rcpp_one_lambda(ipx_tilde, x_stan_tilde, 1, tol, max_iter,max_iter_lasso)
  }else{
    score_nomed = getPseudoScoreLD_Rcpp(ipx_tilde, x_stan_tilde)
  }
  if(l1_penalty>0){
    score_wmed = getScore_Rcpp_one_lambda(ipx, x_stan, 1, tol, max_iter,max_iter_lasso)
  }else{
    score_nomed = getPseudoScoreLD_Rcpp(ipx, x_stan)
  }

  return(list(no_med=score_nomed,w_med=score_wmed))
}



meddic_get_results <- function(out.meddic,scores,tol=1e-5, max_iter=20, max_iter_lasso=100){
  # A wrapper function that execute Meddic model using the input data and scores.
  # Inputs:
  # out.meddic: the output of meddic_data_prep
  # scores: a list with two elements no_med and w_med,which are for the scores for the marginal outcome model and the conditional outcome model, respectively.
  # Outputs:
  # The output of meddic
  y_stan=out.meddic$y_stan
  x_stan=out.meddic$x_stan
  ipy=out.meddic$ipy
  ipx=out.meddic$ipx
  x_stan_tilde=out.meddic$x_stan_tilde
  ipy_tilde=out.meddic$ipy_tilde
  ipx_tilde=out.meddic$ipx_tilde
  l1_penalty=out.meddic$l1_penalty
  l1_penalty_tilde=out.meddic$l1_penalty_tilde
  w_lambda_tilde=out.meddic$w_lambda_tilde
  w_lambda=out.meddic$w_lambda
  p0=out.meddic$p0
  scale.coef <- out.meddic$scale.coef
  score_nomed = scores$no_med
  score_wmed =scores$w_med
  # run the main algorithm
  out.meddic2 <- meddic_Rcpp(y_stan = y_stan,x_stan = x_stan,ipy = ipy,ipx = ipx,x_stan_tilde = x_stan_tilde,ipy_tilde = ipy_tilde,ipx_tilde = ipx_tilde, l1_penalty = l1_penalty,l1_penalty_tilde = l1_penalty_tilde,score_nomed=score_nomed,score_wmed=score_wmed,w_lambda = w_lambda,w_lambda_tilde = w_lambda_tilde, tol=tol, max_iter=max_iter, max_iter_lasso=max_iter_lasso)
  # adding more intermediate products to the output
  out.meddic <- c(out.meddic,out.meddic2)
  out.meddic$tol=tol
  out.meddic$max_iter=max_iter
  out.meddic$max_iter_lasso=max_iter_lasso
  # calculate z scores for the exposures, excluding confounders.
  zscores.element.wise = data.frame(tot=out.meddic$coef_tot/sqrt(diag(out.meddic$var_tot)),direct=out.meddic$coef_direct/sqrt(diag(out.meddic$var_direct)),indirect=out.meddic$coef_indirect/sqrt(diag(out.meddic$var_indirect)))
  if(p0>0){
    zscores.element.wise = zscores.element.wise[-(1:p0),]
  }
  out.meddic$z.elem = zscores.element.wise
  # calculate the coefficients and the standard error for the covariates before standardization.
  if(p0>0){
    out.meddic$coef.est = data.frame(tot=scale.coef*out.meddic$coef_tot[-(1:p0)],direct=scale.coef*out.meddic$coef_direct[-(1:p0)],indirect=scale.coef*out.meddic$coef_indirect[-(1:p0)])
    out.meddic$coef.se = data.frame(tot=scale.coef*sqrt(diag(out.meddic$var_tot)[-(1:p0)]),direct=scale.coef*sqrt(diag(out.meddic$var_direct)[-(1:p0)]),indirect=scale.coef*sqrt(diag(out.meddic$var_indirect)[-(1:p0)]))
  }else{
    out.meddic$coef.est = data.frame(tot=scale.coef*out.meddic$coef_tot,direct=scale.coef*out.meddic$coef_direct,indirect=scale.coef*out.meddic$coef_indirect)
    out.meddic$coef.se = data.frame(tot=scale.coef*sqrt(diag(out.meddic$var_tot)),direct=scale.coef*sqrt(diag(out.meddic$var_direct)),indirect=scale.coef*sqrt(diag(out.meddic$var_indirect)))
  }
  return(out.meddic)
}




inference.elem <- function(meddic,conf.level=0.95){
  # This function reports the p-values and the confidence intervals for the meddic output.
  z.asymp = meddic$z.elem
  q = nrow(z.asymp)
  pval.asymp = data.frame(do.call('cbind',lapply(z.asymp, function(pp) 2*(1-pnorm(abs(pp))))))
  qthr=qnorm((1+conf.level)/2)
  ci.asymp = data.frame(totLB=meddic$coef.est$tot-qthr*meddic$coef.se$tot,totUB=meddic$coef.est$tot+qthr*meddic$coef.se$tot,directLB=meddic$coef.est$direct-qthr*meddic$coef.se$direct,directUB=meddic$coef.est$direct+qthr*meddic$coef.se$direct,indirectLB=meddic$coef.est$indirect-qthr*meddic$coef.se$indirect,indirectUB=meddic$coef.est$indirect+qthr*meddic$coef.se$indirect)
  output = list(conf.level=conf.level,pval.asymp=pval.asymp,ci.asymp=ci.asymp)
  return(output)

}


meddic_Wrapper <- function(Y,X0,Z,M,thr.lse=0.2,no.penalty.exposure=0,no.penalty.confounder=0,tol=1e-5, max_iter=20, max_iter_lasso=100,f_score=NULL){
  # A wrapper function for Meddic that combine all steps.
  # Inputs:
  # Y: the outcome vector
  # X0: the confounder matrix
  # Z: the matrix of exposures
  # M: the matrix of mediators
  # thr.lse: a threshold of (p+q)/n or q/n to decide when to use OLS and when to use the debiased alsso. The default is 0.2, meaning that the adaptive lasso is used if the ratio between the numbers of covariates and the sample size is above 0.2.
  # score_nomed, score_wmed: the score objects for the debiased lasso of the marginal outcome model and the conditional outcome model, if they have been pre-computed. The default value is 0, meaning that they will be computed within this function call.
  # no.penalty.exposure, no.penalty.confounder: if not zero, the exposure or the confounders will not be penalized in the scaled lasso regression. For example, when the exposure is low dimensional and the confounders is high dimensional, it may be helpful to set no.penalty.exposure to be 1 and no.penalty.confounder to be 0.
  # f_score: the file name that the scores are saved in. If NULL, the scores are not saved
  # Outputs:
  # The meddic output

  # prepare the data for meddict analysis
  dt.meddic <- meddic_data_prep(Y,X0,Z,M,thr.lse,no.penalty.exposure,no.penalty.confounder)
  # calculate the scores
  score.meddic <- meddic_get_score(dt.meddic,tol, max_iter, max_iter_lasso)
  # The score calculation for debaised lasso is the most time consuming step, especially when p and q are large. We recommend to save the results in a file.
  if(!is.null(f_score)){
    save(score.meddic,file=f_score)
    print(paste0('The scores for the debiased lasso are saved in file ',f_score))
  }
  # Estimating the direct, indirect and total effects
  out.meddic <- meddic_get_results(dt.meddic,score.meddic,tol, max_iter, max_iter_lasso)
  # Performa inference at teh level alpha=0.05, without adjusting for multiple testing.
  out.inference <- inference.elem(out.meddic)
  out.meddic$results_inference_alpha05 = out.inference
  return(out.meddic)
}


