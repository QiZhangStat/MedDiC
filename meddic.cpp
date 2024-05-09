#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
#include <Rcpp.h>
#include <cmath>
#include <algorithm>
using namespace std;
using namespace Rcpp;




// [[Rcpp::export]]
arma::mat stan_col_Rcpp(arma::mat x){
  // column-wise standardization for a matrix
  // Inputs:
  // x: a matrix
  //Output:
  // The standardized matrix for x.
  int n=x.n_rows;
  int p=x.n_cols;
  arma::mat x_stan=arma::zeros(n,p);
  for(size_t i = 0; i < p; i++){
    arma::vec x_cent_col_i = x.col(i) - mean(x.col(i));
    NumericVector tmpx_i = wrap(x_cent_col_i);
    double x_l2norm_col_i = sqrt(sum(pow(tmpx_i,2)));
    if(x_l2norm_col_i == 0){
      stop("some variables have 0 variance, please remove them at first.");
    }
    x_stan.col(i)=sqrt(n)*x_cent_col_i/x_l2norm_col_i;
  }
  return x_stan;
}

// [[Rcpp::export]]
arma::mat inner_prod_XX_Rcpp(arma::mat x_stan){
  // Inner product of a matrix.
  // Inputs:
  // x_stan: a matrix
  //Output:
  // t(x_stan)%*%x_stan
  arma::mat IP_XX=x_stan.t()*x_stan;
  return IP_XX;
}

// [[Rcpp::export]]
arma::mat inner_prod_XY_Rcpp(arma::mat x_stan, arma::mat y_stan){
  // Inner product of two matrices
  // Inputs:
  // x_stan: a matrix
  //y_stan: a second matrix
  //Output:
  // t(x_stan)%*%y_stan.
  arma::mat IP_XX=x_stan.t()*y_stan;
  return IP_XX;
}


Rcpp::IntegerVector randN(int N, int nsamp) {
  // Randomly sample from integers 1 to N with replacement
  // Inputs:
  // N: Largest integer to sample from.
  // nsamp: number of samples from [1,2,...,N] with replacement to obtain.
  // Output:
  // sampsinteger: nsamp-length vector of samples from [1,2,...,N] with replacement to obtain.
  Rcpp::NumericVector samps = ceiling(runif(nsamp)*N);
  Rcpp::IntegerVector sampsinteger = as<Rcpp::IntegerVector>(samps);
  return sampsinteger;
}


arma::vec rowSums(const arma::mat & X){
  // RCPP version of rowsums()
  int nRows = X.n_rows;
  arma::vec out(nRows);
  for(int i = 0; i < nRows; i++){
    out(i) = sum(X.row(i));
  }
  return out;
}


arma::vec colSums(const arma::mat & X){
  // RCPP version of colSums()
  int nCols = X.n_cols;
  arma::vec out(nCols);
  for(int i = 0; i < nCols; i++){
    out(i) = sum(X.col(i));
  }
  return out;
}

int vecminInd(NumericVector x) {
  // RCPP version of R function which.min()
  // Rcpp supports STL-style iterators
  NumericVector::iterator it = std::min_element(x.begin(), x.end());
  // we want the value so dereference
  return it - x.begin();
//  return 0;
}

int vecfirsttrueInd(LogicalVector x) {
  // Find the first true value of a logical vector. RCPP version of which()[1]
  int output = x.length()-1;
  for(int i = 0; i < x.length(); i++){
    if(x[i]){
      output= i;
      break;
    }
  }
  return output;
}



Rcpp::List lm_Rcpp(arma::vec ipy, arma::mat ipx, arma::vec y_stan,arma::mat x_stan){
  // RCPP verstion of lm(), if the intercept is needed, it has to be part of x_stan. For this project, all covariates have been standardized.
  // Inputs:
  // ipy: the inner product of x_stan and y_stan
  // ipx: the inner product of x_stan
  // y_stan: the vector of response variable
  // x_stan: the matrix of covariates
  //Output:
  // A list that includes the coefficients, the variance of the estimates, the estimate of the noise level, the predicted responses and the residuals
  arma::vec est = solve(ipx,ipy);
  arma::vec ypred = x_stan*est;
  arma::vec epsi = y_stan-ypred;
  Rcpp::NumericVector tmp_epsi = wrap(epsi);
  double sigma = sqrt(sum(pow(tmp_epsi,2))/(x_stan.n_rows-x_stan.n_cols));
  arma::mat var = pow(sigma,2)*inv(ipx);
  Rcpp::List output;
  output.push_back(est,"coef");
  output.push_back(var,"var");
  output.push_back(ypred,"ypred");
  output.push_back(epsi,"residual");
  output.push_back(sigma,"sigma");

  return output;
}


/* ****** fast Lasso regression from SILGGM package. The same function has been used in FastGGM package.******* */
arma::vec FastLasso_Rcpp(arma::vec ipy, arma::mat ipx, double lambda, size_t N,double tol=1e-5,int max_iter_lasso=200){
  // This function modified from the fast lasso regression function from SILGGM package. The same function has been used in FastGGM package. We have modified it so utilize RcppArmadillo.
  // Inputs:
  // ipy: the inner product of x_stan and y_stan
  // ipx: the inner product of x_stan
  // lambda: the penalization parameter
  // N: sample size
  //Output:
  // the regression coefficient vector
  //double tol = 1e-3; //threshold for stopping
  size_t p = ipy.n_elem;
  double dbm;  //maximum of beta difference for while loop
  arma::vec beta=arma::zeros(p); //initialize output beta vector with length p and filled with 0.0
  arma::vec gc=arma::zeros(p); //initialize grandient components
  //update beta vector with coordinate descent, covariance updates
  size_t iter = 0;
  do{
    dbm = 0;
    for(size_t j = 0; j < p; j++){
      double z = (ipy(j) - gc(j))/N + beta(j);
      double beta_tmp = max(0.0, z - lambda) - max(0.0, -z - lambda);
      double diffBeta = beta_tmp - beta(j);
      double diffabs = abs(diffBeta);
      if (diffabs > 0){
        beta(j) = beta_tmp;
        gc = gc+ipx.col(j)*diffBeta;
        dbm = max(dbm, diffabs);
      }
    }
  iter++;
  }
  while(dbm >= tol&&iter<max_iter_lasso);

  return beta;
}


/* ****** fast Lasso regression allowing weighted penalty.******* */
arma::vec FastLasso_W_Rcpp(arma::vec ipy, arma::mat ipx, arma::vec lambda, size_t N,double tol=1e-5,int max_iter_lasso=200){
  // This function is the same as FastLasso_Rcpp except that it allows weighted penalty through setting the input lambda as a vector.
  // Inputs:
  // ipy: the inner product of x_stan and y_stan
  // ipx: the inner product of x_stan
  // lambda: the penalization parameter vector.
  // N: sample size
  //Output:
  // The estimated regression coefficient vector
  //double tol = 1e-3; //threshold for stopping
  size_t p = ipy.n_elem;
  double dbm;  //maximum of beta difference for while loop
  arma::vec beta=arma::zeros(p); //initialize output beta vector with length p and filled with 0.0
  arma::vec gc=arma::zeros(p); //initialize grandient components

  //update beta vector with coordinate descent, covariance updates

  size_t iter = 0;
  do{
    dbm = 0;
    for(size_t j = 0; j < p; j++){
      double z = (ipy(j) - gc(j))/N + beta(j);
      double beta_tmp = max(0.0, z - lambda(j)) - max(0.0, -z - lambda(j));
      double diffBeta = beta_tmp - beta(j);
      double diffabs = abs(diffBeta);
      if (diffabs > 0){
        beta(j) = beta_tmp;
        gc = gc+ipx.col(j)*diffBeta;
        dbm = max(dbm, diffabs);
      }
    }
    iter++;
  }
  while(dbm >= tol&&iter<max_iter_lasso);

  return beta;
}



Rcpp::List scaledLasso_Rcpp(arma::vec ipy, arma::mat ipx, arma::vec y_stan, arma::mat x_stan, double lambda_scale=1, double tol=1e-5, int max_iter=10, int max_iter_lasso=100){
  // This is our implementation of the scaled lasso
  // Inputs:
  // ipy: the inner product of x_stan and y_stan
  // ipx: the inner product of x_stan
  // y_stan: the vector of response variable
  // x_stan: the matrix of covariates
  // lambda_scale: a multiplier to the default tuning parameter sqrt(2 * log(p)/n)
  // max_iter: the maximal number of iterations of the scaled lasso update
  // max_iter_lasso: the max iteration for the lasso algorithm performed within each iteration of the scaled lasso update.
  //Output:
  // A list that includes the coefficients, the variance of the estimates, the estimate of the noise level, the predicted responses and the residuals, and the number of iterations used.

  int p=ipy.n_elem;
  int n=y_stan.n_elem;
  //penalty parameter for L1 norm of betas
  //  if(NumericVector::is_na(lambda)){
  double lambda = lambda_scale*sqrt(2 * log(p)/n);
//  cout << "Use lambda = lambda_scale*sqrt(2*log(p)/n)" << endl;
  //  }
//  cout << "In this case, lambda = " << lambda << endl;

  double sigma = 1.0;
  arma::vec beta=arma::zeros(p);
  arma::vec ypred=arma::zeros(n);
  arma::vec epsi=arma::zeros(n);
  double reg, diffBeta, sigma_tmp, diffSigma;

  //Scaled_lasso
  size_t iter = 0;
  do {
    //Update beta when fix sigma
    reg = sigma * lambda;
    arma::vec beta_tmp = FastLasso_Rcpp(ipy, ipx, reg, n,tol,max_iter_lasso);
    diffBeta = max(abs(beta_tmp - beta));
    beta = beta_tmp;
    //Update sigma when fix beta
    ypred=x_stan*beta;
    epsi = y_stan - ypred;
    Rcpp::NumericVector tmp_epsi = wrap(epsi);
    sigma_tmp = sqrt(sum(pow(tmp_epsi,2))/n);
    diffSigma = abs(sigma_tmp - sigma);
    sigma = sigma_tmp;
    iter++;
  }
  while((diffSigma >= tol || diffBeta >= tol) && iter < max_iter); //Stop iteration

  Rcpp::List result;

  result.push_back(beta,"coef");
  result.push_back(ypred,"ypred");
  result.push_back(epsi,"residual");
  result.push_back(sigma,"sigma");
  result.push_back(iter,"iter");

  return result;

}



// [[Rcpp::export]]
Rcpp::List scaledLasso_W_Rcpp(arma::vec ipy, arma::mat ipx, arma::vec y_stan, arma::mat x_stan, arma::vec w_lambda = 0, double lambda_scale=1, double tol=1e-5, int max_iter=10, int max_iter_lasso=100){
  // This is our implementation of the scaled lasso allowing different covariates to be penalized differently.
  // Inputs:
  // ipy: the inner product of x_stan and y_stan
  // ipx: the inner product of x_stan
  // y_stan: the vector of response variable
  // x_stan: the matrix of covariates
  // w_lambda: a length p multiplier vector to the tuning parameter lambda so that different covariates are penalized differently.
  // lambda_scale: a multiplier to the default tuning parameter sqrt(2 * log(p)/n)
  // max_iter: the maximal number of iterations of the scaled lasso update
  // max_iter_lasso: the max iteration for the lasso algorithm performed within each iteration of the scaled lasso update.
  //Output:
  // A list that includes the coefficients, the variance of the estimates, the estimate of the noise level, the predicted responses and the residuals, and the number of iterations used.
  int p=ipy.n_elem;
  int n=y_stan.n_elem;
  if(w_lambda.n_elem<2){
    w_lambda = arma::ones(p);
  }
  //penalty parameter for L1 norm of betas
  double lambda = lambda_scale*sqrt(2 * log(p)/n);

  double sigma = 1.0;
  arma::vec beta=arma::zeros(p);
  arma::vec ypred=arma::zeros(n);
  arma::vec epsi=arma::zeros(n);
  double reg, diffBeta, sigma_tmp, diffSigma;

  //Scaled_lasso
  size_t iter = 0;
  do {
    //Update beta when fix sigma
    reg = sigma * lambda;
    arma::vec beta_tmp = FastLasso_W_Rcpp(ipy, ipx, reg*w_lambda, n,tol,max_iter_lasso);
    diffBeta = max(abs(beta_tmp - beta));
    beta = beta_tmp;
    //Update sigma when fix beta
    ypred=x_stan*beta;
    epsi = y_stan - ypred;
    Rcpp::NumericVector tmp_epsi = wrap(epsi);
    sigma_tmp = sqrt(sum(pow(tmp_epsi,2))/n);
    diffSigma = abs(sigma_tmp - sigma);
    sigma = sigma_tmp;
    iter++;
  }
  while((diffSigma >= tol || diffBeta >= tol) && iter < max_iter); //Stop iteration
  //cout << "sigma_while" << sigma<< endl;
  NumericVector betawrap = wrap(beta);
  Rcpp::LogicalVector idnzlgc = (betawrap!=0);
  double pnz = sum(idnzlgc);
  if(pnz<n*0.5){
    arma::rowvec idnz = as<arma::rowvec>(idnzlgc);
    arma::mat tmpipx0=ipx.cols(find(idnz));
    arma::mat tmpipx=tmpipx0.rows(find(idnz));
    Rcpp::List lm_temp= lm_Rcpp(ipy(find(idnz)), tmpipx, y_stan, x_stan.cols(find(idnz)));
    sigma = lm_temp["sigma"];
    // use LS after scaled lasso for sigma estimation because it is less biased upward.
//    cout << "Sparsity of beta = " << pnz << " post-selection-LS is used in estimating sigma" << endl;
  }
//  cout << "sigma_ls = " << sigma<< endl;

  Rcpp::List result;
  double sigmaout = sigma;
//  cout << "sigma_output = " << sigma<< endl;
  result.push_back(beta,"coef");
  result.push_back(ypred,"ypred");
  result.push_back(epsi,"residual");
  result.push_back(sigmaout,"sigma");
  result.push_back(iter,"iter");

  return result;

}


// [[Rcpp::export]]
Rcpp::List getScore_Rcpp_one_lambda(arma::mat ipx, arma::mat x_stan, double lambda_scale=1, double tol=1e-5, int max_iter=20, int max_iter_lasso=200){
  // Calculate the score matrix using the scaled lasso for the debiased lasso
  // Inputs:
  // ipx: the inner product of x_stan
  // x_stan: the matrix of covariates
  // lambda_scale: a multiplier to the default tuning parameter sqrt(2 * log(p)/n)
  // max_iter: the maximal number of iterations of the scaled lasso update
  //Output:
  // A list that includes the p-by-p score matrix (W matrix in Zhang and Zhang 2014) and the n-by-p R matrix.
  int p=x_stan.n_cols;
  int n=x_stan.n_rows;
  arma::mat Rmtx=arma::zeros(n,p);
  Rcpp::IntegerVector idall = seq(0,p-1);
  Rcpp::List output;
//  double lambda = lambda_scale*sqrt(2 * log(p)/n);
  for(size_t i = 0; i < p; i++){
    // some technical maneuvers to extract the response and the covariate matrices for the scaled lasso regression.
    Rcpp::LogicalVector idlogical= (idall!=i);
    arma::colvec idcol=as<arma::colvec>(idlogical);
    arma::rowvec idrow=as<arma::rowvec>(idlogical);
    arma::mat tmpxx0=ipx.cols(find(idrow));
    arma::mat tmpxx=tmpxx0.rows(find(idrow));
    arma::mat tmpyx0=ipx.row(i);
    arma::mat tmpyx=tmpyx0.cols(find(idcol));
    arma::mat tmp_x=x_stan.cols(find(idcol));
    // regression one predictor on all other predictors using the scaled lasso
    Rcpp::List scaledlasso_ii = scaledLasso_Rcpp(arma::vectorise(tmpyx), tmpxx, arma::vectorise(x_stan.col(i)),tmp_x, lambda_scale, tol, max_iter, max_iter_lasso);
    arma::colvec resi =as<arma::colvec>(scaledlasso_ii["residual"]);
    resi = resi/inner_prod_XY_Rcpp(resi,x_stan.col(i))(0,0);
    Rmtx.col(i)=resi;
  }
  arma::mat score = Rmtx.t()*x_stan;
  score.diag(0) = arma::zeros(p);
  output.push_back(score,"score");
  output.push_back(Rmtx,"resid");
  return output;
}




// [[Rcpp::export]]
Rcpp::List debiasedlasso_Rcpp(arma::vec y_stan, Rcpp::List score, Rcpp::List init_est){
  // Debiased lasso as in Zhang and Zhang (2014) after the scores and the initial estimates are calculated.
  // Inputs:
  // y_stan: the response vector
  // score: the output of function getScore_Rcpp_one_lambda()
  // init_est: the output of an regression model such as the scaled lasso.
  //Output:
  // A list that includes the estimated coefficients, and the p-by-p variance-covariance matrix.
  arma::mat Rmtx = score["resid"];
  arma::mat Wmtx = score["score"];
  arma::vec beta_init = init_est["coef"];
  double sigma = init_est["sigma"];
  arma::vec est = Rmtx.t()*y_stan - Wmtx*beta_init;
  arma::mat var = pow(sigma,2)*Rmtx.t()*Rmtx;
  Rcpp::List output;
  output.push_back(est,"coef");
  output.push_back(var,"var");

  return output;
}

// [[Rcpp::export]]
Rcpp::List getPseudoScoreLD_Rcpp(arma::mat ipx, arma::mat x_stan){
  // This function calculates the score matrix for the "debiased" estimator when the dimension is low and when the initial estimator is based on OLS. It allows the debiased lasso function keep the same input/output structure when the dimension is low, and return OLS results as the output.
  // Inputs:
  // ipx: the inner product of x_stan
  // x_stan: the matrix of covariates
  //Output:
  // A list that includes of OLS version of the p-by-p score matrix (W matrix in Zhang and Zhang 2014) and the n-by-p R matrix.
  int p=x_stan.n_cols;
  Rcpp::List output;
  arma::mat Rmtx = solve(ipx,x_stan.t()).t();
  arma::mat score = arma::zeros(p,p);
  output.push_back(Rmtx,"resid");
  output.push_back(score,"score");
  return output;
}





// [[Rcpp::export]]
Rcpp::List meddic_Rcpp(arma::vec y_stan, arma::mat x_stan, arma::vec ipy, arma::mat ipx, arma::mat x_stan_tilde, arma::vec ipy_tilde,arma::mat ipx_tilde, double l1_penalty = 1,double l1_penalty_tilde = 1, Rcpp::List score_nomed=0,Rcpp::List score_wmed=0, arma::vec w_lambda_tilde=0, arma::vec w_lambda=0, double tol=1e-5, int max_iter=10, int max_iter_lasso=100){
  // This is the main function for MedDiC
  // Inputs:
  // y_stan: the vector of outcome variable
  // x_stan: the matrix of the confounders, exposures and mediators
  // ipy: the inner product of x_stan and y_stan
  // ipx: the inner product of x_stan
  // x_stan_tilde: the matrix of confounders and exposures.
  // ipy_tilde: the inner product of x_stan_tilde and y_stan
  // ipx_tilde: the inner product of x_stan_tilde
  // l1_penalty: 1 if the scaled lasso will be used for the conditional outcome model (regressing the outcome on the confoudners, exposures and mediators). If it is not 1, OLS will be used.
  // l1_penalty_tilde: 1 if the scaled lasso will be used for the conditional outcome model (regressing the outcome on the confoudners and exposures).  If it is not 1, OLS will be used.
  //socre_nomed, score_wmed: the debiased lasso scores for the marginal outcome model and  the conditional outcome model, usually the output of the function getScore_Rcpp_one_lambda. If they have not been pre-computed, please use input 0, and they will be calculated.
  // w_lambda, w_lambda_tilde: the weight vectors for the conditional outcome model and the marginal outcome model so that different covariates are penalized differently in the initial fit of adaptive scaled lasso. For example, if the exposure is high dimensional, and the confounders is low dimensional, one may let the weights for the confounders to be 0 so that they are not penalized.
  // max_iter: the maximal number of iterations of the scaled lasso update
  // max_iter_lasso: the max iteration for the lasso algorithm performed within each iteration of the scaled lasso update.
  //Output:
  // A list that includes the estimates and the variance for the total, direct and the indirect effects and the M-->Y link, the noise level, and the initial estimates, score matrices and he weight vector used.

  int p_tilde = ipy_tilde.n_elem;
  int p = ipy.n_elem;
  int n = y_stan.n_elem;
  // Define the default weights used for the initial fit of the adaptive scaled lasso.
  if(w_lambda_tilde.n_elem<2){
    w_lambda_tilde = arma::ones(p_tilde);
  }
  if(w_lambda.n_elem<2){
    w_lambda = arma::ones(p);
  }
  // A small number to be added to the initial fit of the scaled lasso to define the weights for the adaptive scaled lasso.
  double thr_nomed = stddev(y_stan)/sqrt(n);
  double thr_wmed = stddev(y_stan)/sqrt(n);

  Rcpp::List init_nomed0;
  Rcpp::List init_nomed;
  Rcpp::List debiased_nomed;

  Rcpp::List init_wmed0;
  Rcpp::List init_wmed;
  Rcpp::List debiased_wmed;
  if(l1_penalty_tilde!=1){
    // If L1 penalty is not needed, OLS replaces the debaised lasso for the marginal outcome model
    cout << "Least squares is used for fitting the model without mediators" << endl;
    score_nomed = getPseudoScoreLD_Rcpp(ipx_tilde, x_stan_tilde);
    debiased_nomed = lm_Rcpp(ipy_tilde, ipx_tilde, y_stan, x_stan_tilde);
    init_nomed = debiased_nomed;
  } else {
    cout << "Debiased lasso is used for fitting the model without mediators" << endl;
    cout << "Scaled lasso with universal tuning parameter sqrt(2*log(p)/n) is used for the initial estimate of the model without mediators" << endl;
    // Initial fit using the scaled lasso
    init_nomed0 = scaledLasso_W_Rcpp(ipy_tilde, ipx_tilde, y_stan, x_stan_tilde,w_lambda_tilde, 1, tol, max_iter, max_iter_lasso);
    arma::vec beta_init0 = init_nomed0["coef"];
    // Define the weights for the final fit of the adaptive scaled lasso.
    beta_init0 = abs(beta_init0)+thr_nomed;
    arma::vec w_lambda_tilde_adpt = w_lambda_tilde/beta_init0;
    w_lambda_tilde_adpt = w_lambda_tilde_adpt*sum(w_lambda_tilde)/sum(w_lambda_tilde_adpt);
    // fitting the adaptive scaled lasso to get the initial fit for the coefficients in the marginal outcome model
    init_nomed = scaledLasso_W_Rcpp(ipy_tilde, ipx_tilde, y_stan, x_stan_tilde,w_lambda_tilde_adpt, 1, tol, max_iter, max_iter_lasso);

    // If the score has not been pre-computed, get_Score_Rcpp_one_lambda function is used here.
    if(score_nomed.length()<2){
      cout << "Need to calculate scores for debiased lasso for the model without mediators. It may take awhile" << endl;
      cout << "Scaled lasso with universal tuning parameter sqrt(2*log(p)/n) is used for calculating the scores of the model without mediators" << endl;
      score_nomed = getScore_Rcpp_one_lambda(ipx_tilde, x_stan_tilde, 1, tol, max_iter,max_iter_lasso);
    }
    // debiased lasso estimate for the marginal outcome model. It includes the total effects.
    debiased_nomed = debiasedlasso_Rcpp(y_stan,score_nomed,init_nomed);
  }

  if(l1_penalty!=1){
    // If L1 penalty is not needed, OLS replaces the debaised lasso for the conditional outcome model
    cout << "Least squares is used for fitting the model with mediators" << endl;
    score_wmed = getPseudoScoreLD_Rcpp(ipx, x_stan);
    debiased_wmed = lm_Rcpp(ipy, ipx, y_stan, x_stan);
    init_wmed = debiased_wmed;
  } else {
    cout << "Debiased lasso is used for fitting the model with mediators" << endl;
    cout << "Scaled lasso with universal tuning parameter sqrt(2*log(p)/n) is used for the initial estimate of the model with mediators" << endl;
    // Initial fit using the scaled lasso
    init_wmed0 = scaledLasso_W_Rcpp(ipy, ipx, y_stan, x_stan,w_lambda,1, tol, max_iter, max_iter_lasso);
    arma::vec beta_init0 = init_wmed0["coef"];
    // Define the weights for the final fit of the adaptive scaled lasso.
    beta_init0 = abs(beta_init0)+thr_wmed;
    arma::vec w_lambda_adpt = w_lambda/beta_init0;
    w_lambda_adpt = w_lambda_adpt*sum(w_lambda)/sum(w_lambda_adpt);
    // fitting the adaptive scaled lasso to get the initial fit for the coefficients in the conditional outcome model
    init_wmed = scaledLasso_W_Rcpp(ipy, ipx, y_stan, x_stan,w_lambda_adpt,1, tol, max_iter, max_iter_lasso);
    // If the score has not been pre-computed, get_Score_Rcpp_one_lambda function is used here.
    if(score_wmed.length()<2){
      cout << "Need to calculate scores for debiased lasso for the model with mediators. It may take awhile." << endl;
      cout << "Scaled lasso with universal tuning parameter sqrt(2*log(p)/n) is used for calculating the scores of the model with mediators" << endl;
      score_wmed = getScore_Rcpp_one_lambda(ipx, x_stan, 1, tol, max_iter,max_iter_lasso);
    }
    // debiased lasso estimate for the conditional outcome model. It includes the direct effects.
    debiased_wmed = debiasedlasso_Rcpp(y_stan,score_wmed,init_wmed);
  }
  // extract the components need for computing the indirect effect
  arma::vec coef_nomed = debiased_nomed["coef"];
  arma::vec coef_wmed = debiased_wmed["coef"];
  arma::mat var_nomed = debiased_nomed["var"];
  arma::mat var_wmed = debiased_wmed["var"];

  arma::mat Rmtx_nomed = score_nomed["resid"];
  arma::mat Rmtx_wmed = score_wmed["resid"];
  arma::vec epsi_nomed = init_nomed["residual"];
  arma::vec epsi_wmed = init_wmed["residual"];
  double sigma_nomed = init_nomed["sigma"];
  double sigma_wmed = init_wmed["sigma"];

  int q=Rmtx_nomed.n_cols;
  Rcpp::IntegerVector idall = seq(0,Rmtx_wmed.n_cols-1);
  Rcpp::LogicalVector idq= (idall<q);
  arma::colvec idq_col=as<arma::colvec>(idq);
  arma::rowvec idq_row=as<arma::rowvec>(idq);
  arma::vec idq_vec = as<arma::vec>(idq);
  arma::mat Rmtx_wmed_q = Rmtx_wmed.cols(find(idq_col));
  arma::vec beta_wmed = coef_wmed(find(idq_vec));
  arma::mat var_direct0 = var_wmed.cols(find(idq_col));
  arma::mat var_direct = var_direct0.rows(find(idq_row));

  Rcpp::LogicalVector idm= (idall>=q);
  arma::colvec idm_col=as<arma::colvec>(idm);
  arma::rowvec idm_row=as<arma::rowvec>(idm);
  arma::vec idm_vec = as<arma::vec>(idm);
  arma::vec coef_med = coef_wmed(find(idm_vec));
  arma::mat var_med0 = var_wmed.cols(find(idm_col));
  arma::mat var_med = var_med0.rows(find(idm_row));

  Rcpp::NumericVector cross_prod_epsi = wrap(epsi_nomed%epsi_wmed);
  Rcpp::NumericVector prod_epsi_nomed = wrap(epsi_nomed%epsi_nomed);
  Rcpp::NumericVector prod_epsi_wmed = wrap(epsi_wmed%epsi_wmed);
  double corr_epsi = sum(cross_prod_epsi)/sqrt(sum(prod_epsi_nomed)*sum(prod_epsi_wmed));
  Rcpp::NumericVector vars = NumericVector::create(pow(sigma_nomed,2),pow(sigma_wmed,2),corr_epsi*sigma_wmed*sigma_nomed);
  // The following are the indirect effect and its variance matrix.
  arma::vec coef_indirect = coef_nomed-beta_wmed;
  arma::mat var_indirect = var_nomed+var_direct-vars[2]*(Rmtx_nomed.t()*Rmtx_wmed_q+Rmtx_wmed_q.t()*Rmtx_nomed);

  Rcpp::List output;
  output.push_back(coef_nomed,"coef_tot");
  output.push_back(var_nomed,"var_tot");
  output.push_back(beta_wmed,"coef_direct");
  output.push_back(var_direct,"var_direct");
  output.push_back(coef_indirect,"coef_indirect");
  output.push_back(var_indirect,"var_indirect");
  output.push_back(coef_med,"coef_mediator");
  output.push_back(var_med,"var_mediator");
  output.push_back(vars,"var_noise");
  output.push_back(init_nomed,"init_est_no_med");
  output.push_back(init_wmed,"init_est_w_med");
  output.push_back(score_nomed,"score_no_med");
  output.push_back(score_wmed,"score_w_med");
  output.push_back(w_lambda,"w_lambda_no_med");
  output.push_back(w_lambda,"w_lambda_w_med");
  return output;
}



