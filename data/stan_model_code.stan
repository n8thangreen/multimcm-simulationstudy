functions {
#include /include/distributions.stan
}

data {
 int<lower=1> nTx;
int<lower=1, upper=3> cf_model;         // cure fraction
int<lower=0> N_1;
array[nTx] int<lower=0> n_1;
int<lower=0> H_1;
vector<lower=0>[N_1] t_1;
vector<lower=0, upper=1>[N_1] d_1;
matrix[N_1, H_1] X_1;
vector[H_1] mu_S_1;
vector<lower=0>[H_1] sigma_S_1;
real<lower=0> a_shape_1;
real<lower=0> b_shape_1;
 int<lower=0> N_2;
array[nTx] int<lower=0> n_2;
int<lower=0> H_2;
vector<lower=0>[N_2] t_2;
vector<lower=0, upper=1>[N_2] d_2;
matrix[N_2, H_2] X_2;
vector[H_2] mu_S_2;
vector<lower=0>[H_2] sigma_S_2;
real<lower=0> a_shape_2;
real<lower=0> b_shape_2;
 int<lower=0> N_3;
array[nTx] int<lower=0> n_3;
int<lower=0> H_3;
vector<lower=0>[N_3] t_3;
vector<lower=0, upper=1>[N_3] d_3;
matrix[N_3, H_3] X_3;
vector[H_3] mu_S_3;
vector<lower=0>[H_3] sigma_S_3;
real<lower=0> a_shape_3;
real<lower=0> b_shape_3;
 int<lower=0> N_4;
array[nTx] int<lower=0> n_4;
int<lower=0> H_4;
vector<lower=0>[N_4] t_4;
vector<lower=0, upper=1>[N_4] d_4;
matrix[N_4, H_4] X_4;
vector[H_4] mu_S_4;
vector<lower=0>[H_4] sigma_S_4;
real<lower=0> a_shape_4;
real<lower=0> b_shape_4;
 int<lower=0> N_5;
array[nTx] int<lower=0> n_5;
int<lower=0> H_5;
vector<lower=0>[N_5] t_5;
vector<lower=0, upper=1>[N_5] d_5;
matrix[N_5, H_5] X_5;
vector[H_5] mu_S_5;
vector<lower=0>[H_5] sigma_S_5;
real<lower=0> a_shape_5;
real<lower=0> b_shape_5;
 int<lower=0> N_6;
array[nTx] int<lower=0> n_6;
int<lower=0> H_6;
vector<lower=0>[N_6] t_6;
vector<lower=0, upper=1>[N_6] d_6;
matrix[N_6, H_6] X_6;
vector[H_6] mu_S_6;
vector<lower=0>[H_6] sigma_S_6;
real<lower=0> a_shape_6;
real<lower=0> b_shape_6;
 int<lower=0> N_7;
array[nTx] int<lower=0> n_7;
int<lower=0> H_7;
vector<lower=0>[N_7] t_7;
vector<lower=0, upper=1>[N_7] d_7;
matrix[N_7, H_7] X_7;
vector[H_7] mu_S_7;
vector<lower=0>[H_7] sigma_S_7;
real<lower=0> a_shape_7;
real<lower=0> b_shape_7;
 int<lower=0> N_8;
array[nTx] int<lower=0> n_8;
int<lower=0> H_8;
vector<lower=0>[N_8] t_8;
vector<lower=0, upper=1>[N_8] d_8;
matrix[N_8, H_8] X_8;
vector[H_8] mu_S_8;
vector<lower=0>[H_8] sigma_S_8;
real<lower=0> a_shape_8;
real<lower=0> b_shape_8;
 int<lower=0> N_9;
array[nTx] int<lower=0> n_9;
int<lower=0> H_9;
vector<lower=0>[N_9] t_9;
vector<lower=0, upper=1>[N_9] d_9;
matrix[N_9, H_9] X_9;
vector[H_9] mu_S_9;
vector<lower=0>[H_9] sigma_S_9;
real<lower=0> a_shape_9;
real<lower=0> b_shape_9;
 int<lower=0> N_10;
array[nTx] int<lower=0> n_10;
int<lower=0> H_10;
vector<lower=0>[N_10] t_10;
vector<lower=0, upper=1>[N_10] d_10;
matrix[N_10, H_10] X_10;
vector[H_10] mu_S_10;
vector<lower=0>[H_10] sigma_S_10;
real<lower=0> a_shape_10;
real<lower=0> b_shape_10;
int<lower=1, upper=2> bg_model;
 vector[bg_model == 1 ? H_1 : 0] mu_bg;
 vector<lower=0>[bg_model == 1 ? H_1 : 0] sigma_bg;
vector<lower=0>[bg_model == 2 ? N_1 : 0] h_bg_1;
vector<lower=0>[bg_model == 2 ? N_2 : 0] h_bg_2;
vector<lower=0>[bg_model == 2 ? N_3 : 0] h_bg_3;
vector<lower=0>[bg_model == 2 ? N_4 : 0] h_bg_4;
vector<lower=0>[bg_model == 2 ? N_5 : 0] h_bg_5;
vector<lower=0>[bg_model == 2 ? N_6 : 0] h_bg_6;
vector<lower=0>[bg_model == 2 ? N_7 : 0] h_bg_7;
vector<lower=0>[bg_model == 2 ? N_8 : 0] h_bg_8;
vector<lower=0>[bg_model == 2 ? N_9 : 0] h_bg_9;
vector<lower=0>[bg_model == 2 ? N_10 : 0] h_bg_10;
 matrix[nTx, nTx] Tx_dmat;         // treatment design matrix
 vector[cf_model == 3 ? nTx : 0] mu_alpha;             // treatment regression coefficients
 vector<lower=0>[cf_model == 3 ? nTx : 0] sigma_alpha;
 int<lower=0> t_max;
vector[cf_model == 2 ? nTx : 0] mu_alpha_1;
vector<lower=0>[cf_model == 2 ? nTx : 0] sigma_alpha_1;
vector[cf_model == 2 ? nTx : 0] mu_alpha_2;
vector<lower=0>[cf_model == 2 ? nTx : 0] sigma_alpha_2;
vector[cf_model == 2 ? nTx : 0] mu_alpha_3;
vector<lower=0>[cf_model == 2 ? nTx : 0] sigma_alpha_3;
vector[cf_model == 2 ? nTx : 0] mu_alpha_4;
vector<lower=0>[cf_model == 2 ? nTx : 0] sigma_alpha_4;
vector[cf_model == 2 ? nTx : 0] mu_alpha_5;
vector<lower=0>[cf_model == 2 ? nTx : 0] sigma_alpha_5;
vector[cf_model == 2 ? nTx : 0] mu_alpha_6;
vector<lower=0>[cf_model == 2 ? nTx : 0] sigma_alpha_6;
vector[cf_model == 2 ? nTx : 0] mu_alpha_7;
vector<lower=0>[cf_model == 2 ? nTx : 0] sigma_alpha_7;
vector[cf_model == 2 ? nTx : 0] mu_alpha_8;
vector<lower=0>[cf_model == 2 ? nTx : 0] sigma_alpha_8;
vector[cf_model == 2 ? nTx : 0] mu_alpha_9;
vector<lower=0>[cf_model == 2 ? nTx : 0] sigma_alpha_9;
vector[cf_model == 2 ? nTx : 0] mu_alpha_10;
vector<lower=0>[cf_model == 2 ? nTx : 0] sigma_alpha_10;
 array[cf_model == 1 ? 1 : 0] real<lower=0> a_cf;

 array[cf_model == 1 ? 1 : 0] real<lower=0> b_cf;
vector[cf_model == 3 ? nTx : 0] mu_sd_cf;
vector<lower=0>[cf_model == 3 ? nTx : 0] sigma_sd_cf;

}

parameters {
real<lower=0> shape_1;
 real<lower=0> shape_2;
 real<lower=0> shape_3;
 real<lower=0> shape_4;
 real<lower=0> shape_5;
 real<lower=0> shape_6;
 real<lower=0> shape_7;
 real<lower=0> shape_8;
 real<lower=0> shape_9;
 real<lower=0> shape_10;// coefficients in linear predictor (including intercept)
 vector[bg_model == 1 ? H_1 : 0] beta_bg;
 vector[cf_model != 2 ? nTx : 0] alpha;
vector[cf_model == 2 ? nTx : 0] alpha_1;
vector[H_1] beta_1;
vector[cf_model == 2 ? nTx : 0] alpha_2;
vector[H_2] beta_2;
vector[cf_model == 2 ? nTx : 0] alpha_3;
vector[H_3] beta_3;
vector[cf_model == 2 ? nTx : 0] alpha_4;
vector[H_4] beta_4;
vector[cf_model == 2 ? nTx : 0] alpha_5;
vector[H_5] beta_5;
vector[cf_model == 2 ? nTx : 0] alpha_6;
vector[H_6] beta_6;
vector[cf_model == 2 ? nTx : 0] alpha_7;
vector[H_7] beta_7;
vector[cf_model == 2 ? nTx : 0] alpha_8;
vector[H_8] beta_8;
vector[cf_model == 2 ? nTx : 0] alpha_9;
vector[H_9] beta_9;
vector[cf_model == 2 ? nTx : 0] alpha_10;
vector[H_10] beta_10;
 vector<lower=0, upper=1>[cf_model == 1 ? nTx : 0] cf_pooled;
vector[cf_model == 3 ? nTx : 0] lp_cf_1;
vector[cf_model == 3 ? nTx : 0] lp_cf_2;
vector[cf_model == 3 ? nTx : 0] lp_cf_3;
vector[cf_model == 3 ? nTx : 0] lp_cf_4;
vector[cf_model == 3 ? nTx : 0] lp_cf_5;
vector[cf_model == 3 ? nTx : 0] lp_cf_6;
vector[cf_model == 3 ? nTx : 0] lp_cf_7;
vector[cf_model == 3 ? nTx : 0] lp_cf_8;
vector[cf_model == 3 ? nTx : 0] lp_cf_9;
vector[cf_model == 3 ? nTx : 0] lp_cf_10;
 vector<lower=0>[cf_model == 3 ? nTx : 0] sd_cf;

}

transformed parameters {
vector[N_1] lp_1_bg;

vector<lower=0>[N_1] lambda_1_bg;

vector[N_2] lp_2_bg;

vector<lower=0>[N_2] lambda_2_bg;

vector[N_3] lp_3_bg;

vector<lower=0>[N_3] lambda_3_bg;

vector[N_4] lp_4_bg;

vector<lower=0>[N_4] lambda_4_bg;

vector[N_5] lp_5_bg;

vector<lower=0>[N_5] lambda_5_bg;

vector[N_6] lp_6_bg;

vector<lower=0>[N_6] lambda_6_bg;

vector[N_7] lp_7_bg;

vector<lower=0>[N_7] lambda_7_bg;

vector[N_8] lp_8_bg;

vector<lower=0>[N_8] lambda_8_bg;

vector[N_9] lp_9_bg;

vector<lower=0>[N_9] lambda_9_bg;

vector[N_10] lp_10_bg;

vector<lower=0>[N_10] lambda_10_bg;
vector[N_1] lp_1;// rate parameters
vector[N_1] lambda_1;
 vector[N_2] lp_2;// rate parameters
vector[N_2] lambda_2;
 vector[N_3] lp_3;// rate parameters
vector[N_3] lambda_3;
 vector[N_4] lp_4;// rate parameters
vector[N_4] lambda_4;
 vector[N_5] lp_5;// rate parameters
vector[N_5] lambda_5;
 vector[N_6] lp_6;// rate parameters
vector[N_6] lambda_6;
 vector[N_7] lp_7;// rate parameters
vector[N_7] lambda_7;
 vector[N_8] lp_8;// rate parameters
vector[N_8] lambda_8;
 vector[N_9] lp_9;// rate parameters
vector[N_9] lambda_9;
 vector[N_10] lp_10;// rate parameters
vector[N_10] lambda_10;
 vector<lower=0, upper=1>[cf_model == 3 ? nTx : 0] cf_global;
vector<lower=0, upper=1>[nTx] cf_1;
vector[cf_model == 2 ? nTx : 0] tx_cf_1;
vector<lower=0, upper=1>[nTx] cf_2;
vector[cf_model == 2 ? nTx : 0] tx_cf_2;
vector<lower=0, upper=1>[nTx] cf_3;
vector[cf_model == 2 ? nTx : 0] tx_cf_3;
vector<lower=0, upper=1>[nTx] cf_4;
vector[cf_model == 2 ? nTx : 0] tx_cf_4;
vector<lower=0, upper=1>[nTx] cf_5;
vector[cf_model == 2 ? nTx : 0] tx_cf_5;
vector<lower=0, upper=1>[nTx] cf_6;
vector[cf_model == 2 ? nTx : 0] tx_cf_6;
vector<lower=0, upper=1>[nTx] cf_7;
vector[cf_model == 2 ? nTx : 0] tx_cf_7;
vector<lower=0, upper=1>[nTx] cf_8;
vector[cf_model == 2 ? nTx : 0] tx_cf_8;
vector<lower=0, upper=1>[nTx] cf_9;
vector[cf_model == 2 ? nTx : 0] tx_cf_9;
vector<lower=0, upper=1>[nTx] cf_10;
vector[cf_model == 2 ? nTx : 0] tx_cf_10;
 vector[cf_model == 3 ? nTx : 0] lp_cf_global;
// correlated event times
  lp_1 = X_1*beta_1;

// background survival with uncertainty

if (bg_model == 1) {
  lp_1_bg = X_1*beta_bg;
} else {
  lp_1_bg = log(h_bg_1);
}

lambda_1_bg = exp(lp_1_bg);

// correlated event times
  lp_2 = X_2*beta_2;

// background survival with uncertainty

if (bg_model == 1) {
  lp_2_bg = X_2*beta_bg;
} else {
  lp_2_bg = log(h_bg_2);
}

lambda_2_bg = exp(lp_2_bg);

// correlated event times
  lp_3 = X_3*beta_3;

// background survival with uncertainty

if (bg_model == 1) {
  lp_3_bg = X_3*beta_bg;
} else {
  lp_3_bg = log(h_bg_3);
}

lambda_3_bg = exp(lp_3_bg);

// correlated event times
  lp_4 = X_4*beta_4;

// background survival with uncertainty

if (bg_model == 1) {
  lp_4_bg = X_4*beta_bg;
} else {
  lp_4_bg = log(h_bg_4);
}

lambda_4_bg = exp(lp_4_bg);

// correlated event times
  lp_5 = X_5*beta_5;

// background survival with uncertainty

if (bg_model == 1) {
  lp_5_bg = X_5*beta_bg;
} else {
  lp_5_bg = log(h_bg_5);
}

lambda_5_bg = exp(lp_5_bg);

// correlated event times
  lp_6 = X_6*beta_6;

// background survival with uncertainty

if (bg_model == 1) {
  lp_6_bg = X_6*beta_bg;
} else {
  lp_6_bg = log(h_bg_6);
}

lambda_6_bg = exp(lp_6_bg);

// correlated event times
  lp_7 = X_7*beta_7;

// background survival with uncertainty

if (bg_model == 1) {
  lp_7_bg = X_7*beta_bg;
} else {
  lp_7_bg = log(h_bg_7);
}

lambda_7_bg = exp(lp_7_bg);

// correlated event times
  lp_8 = X_8*beta_8;

// background survival with uncertainty

if (bg_model == 1) {
  lp_8_bg = X_8*beta_bg;
} else {
  lp_8_bg = log(h_bg_8);
}

lambda_8_bg = exp(lp_8_bg);

// correlated event times
  lp_9 = X_9*beta_9;

// background survival with uncertainty

if (bg_model == 1) {
  lp_9_bg = X_9*beta_bg;
} else {
  lp_9_bg = log(h_bg_9);
}

lambda_9_bg = exp(lp_9_bg);

// correlated event times
  lp_10 = X_10*beta_10;

// background survival with uncertainty

if (bg_model == 1) {
  lp_10_bg = X_10*beta_bg;
} else {
  lp_10_bg = log(h_bg_10);
}

lambda_10_bg = exp(lp_10_bg);
lambda_1 = exp(lp_1);
 lambda_2 = exp(lp_2);
 lambda_3 = exp(lp_3);
 lambda_4 = exp(lp_4);
 lambda_5 = exp(lp_5);
 lambda_6 = exp(lp_6);
 lambda_7 = exp(lp_7);
 lambda_8 = exp(lp_8);
 lambda_9 = exp(lp_9);
 lambda_10 = exp(lp_10);
 if (cf_model == 1) {
	 cf_1 = cf_pooled;
	 cf_2 = cf_pooled;
	 cf_3 = cf_pooled;
	 cf_4 = cf_pooled;
	 cf_5 = cf_pooled;
	 cf_6 = cf_pooled;
	 cf_7 = cf_pooled;
	 cf_8 = cf_pooled;
	 cf_9 = cf_pooled;
	 cf_10 = cf_pooled;
}
if (cf_model == 3) {
	 lp_cf_global = Tx_dmat*alpha;
	 cf_global = inv_logit(lp_cf_global);
	 cf_1 = inv_logit(lp_cf_1);
	 cf_2 = inv_logit(lp_cf_2);
	 cf_3 = inv_logit(lp_cf_3);
	 cf_4 = inv_logit(lp_cf_4);
	 cf_5 = inv_logit(lp_cf_5);
	 cf_6 = inv_logit(lp_cf_6);
	 cf_7 = inv_logit(lp_cf_7);
	 cf_8 = inv_logit(lp_cf_8);
	 cf_9 = inv_logit(lp_cf_9);
	 cf_10 = inv_logit(lp_cf_10);
}
if (cf_model == 2) {
tx_cf_1 = Tx_dmat*alpha_1;
cf_1 = inv_logit(tx_cf_1);
tx_cf_2 = Tx_dmat*alpha_2;
cf_2 = inv_logit(tx_cf_2);
tx_cf_3 = Tx_dmat*alpha_3;
cf_3 = inv_logit(tx_cf_3);
tx_cf_4 = Tx_dmat*alpha_4;
cf_4 = inv_logit(tx_cf_4);
tx_cf_5 = Tx_dmat*alpha_5;
cf_5 = inv_logit(tx_cf_5);
tx_cf_6 = Tx_dmat*alpha_6;
cf_6 = inv_logit(tx_cf_6);
tx_cf_7 = Tx_dmat*alpha_7;
cf_7 = inv_logit(tx_cf_7);
tx_cf_8 = Tx_dmat*alpha_8;
cf_8 = inv_logit(tx_cf_8);
tx_cf_9 = Tx_dmat*alpha_9;
cf_9 = inv_logit(tx_cf_9);
tx_cf_10 = Tx_dmat*alpha_10;
cf_10 = inv_logit(tx_cf_10);
}

}

model {
int idx_1;
int idx_2;
int idx_3;
int idx_4;
int idx_5;
int idx_6;
int idx_7;
int idx_8;
int idx_9;
int idx_10;      
beta_1 ~ normal(mu_S_1, sigma_S_1);
      
beta_2 ~ normal(mu_S_2, sigma_S_2);
      
beta_3 ~ normal(mu_S_3, sigma_S_3);
      
beta_4 ~ normal(mu_S_4, sigma_S_4);
      
beta_5 ~ normal(mu_S_5, sigma_S_5);
      
beta_6 ~ normal(mu_S_6, sigma_S_6);
      
beta_7 ~ normal(mu_S_7, sigma_S_7);
      
beta_8 ~ normal(mu_S_8, sigma_S_8);
      
beta_9 ~ normal(mu_S_9, sigma_S_9);
      
beta_10 ~ normal(mu_S_10, sigma_S_10);
if (bg_model == 1) {
	 beta_bg ~ normal(mu_bg, sigma_bg);
}
shape_1 ~ gamma(a_shape_1, b_shape_1);
 shape_2 ~ gamma(a_shape_2, b_shape_2);
 shape_3 ~ gamma(a_shape_3, b_shape_3);
 shape_4 ~ gamma(a_shape_4, b_shape_4);
 shape_5 ~ gamma(a_shape_5, b_shape_5);
 shape_6 ~ gamma(a_shape_6, b_shape_6);
 shape_7 ~ gamma(a_shape_7, b_shape_7);
 shape_8 ~ gamma(a_shape_8, b_shape_8);
 shape_9 ~ gamma(a_shape_9, b_shape_9);
 shape_10 ~ gamma(a_shape_10, b_shape_10);// cure fraction 
 if (cf_model == 3) {
 	 alpha ~ normal(mu_alpha, sigma_alpha);
 sd_cf ~ cauchy(mu_sd_cf, sigma_sd_cf);  // truncated

lp_cf_1 ~ normal(lp_cf_global, sd_cf);
lp_cf_2 ~ normal(lp_cf_global, sd_cf);
lp_cf_3 ~ normal(lp_cf_global, sd_cf);
lp_cf_4 ~ normal(lp_cf_global, sd_cf);
lp_cf_5 ~ normal(lp_cf_global, sd_cf);
lp_cf_6 ~ normal(lp_cf_global, sd_cf);
lp_cf_7 ~ normal(lp_cf_global, sd_cf);
lp_cf_8 ~ normal(lp_cf_global, sd_cf);
lp_cf_9 ~ normal(lp_cf_global, sd_cf);
lp_cf_10 ~ normal(lp_cf_global, sd_cf);
} else if (cf_model == 2) {
alpha_1 ~ normal(mu_alpha_1, sigma_alpha_1);
alpha_2 ~ normal(mu_alpha_2, sigma_alpha_2);
alpha_3 ~ normal(mu_alpha_3, sigma_alpha_3);
alpha_4 ~ normal(mu_alpha_4, sigma_alpha_4);
alpha_5 ~ normal(mu_alpha_5, sigma_alpha_5);
alpha_6 ~ normal(mu_alpha_6, sigma_alpha_6);
alpha_7 ~ normal(mu_alpha_7, sigma_alpha_7);
alpha_8 ~ normal(mu_alpha_8, sigma_alpha_8);
alpha_9 ~ normal(mu_alpha_9, sigma_alpha_9);
alpha_10 ~ normal(mu_alpha_10, sigma_alpha_10);
} else {
	 cf_pooled ~ beta(a_cf, b_cf);
}
idx_1 = 1;

// likelihood
for (Tx in 1:nTx) {
  for (i in idx_1:(idx_1 + n_1[Tx] - 1)) {

     target += log_sum_exp(
    log(cf_1[Tx]) +
      exp_log_S(t_1[i], lambda_1_bg[i]),
    log1m(cf_1[Tx]) +
      exp_log_S(t_1[i], lambda_1_bg[i]) + weibull_log_S(t_1[i], shape_1, lambda_1[i]));

     target += d_1[i] * log_sum_exp(
    log(lambda_1_bg[i]),
    log1m(cf_1[Tx]) +
      weibull_lpdf(t_1[i] | shape_1, lambda_1[i]) - log(cf_1[Tx] + (1 - cf_1[Tx])*weibull_Surv(t_1[i], shape_1, lambda_1[i])));
  }

  idx_1 = idx_1 + n_1[Tx];
}
 idx_2 = 1;

// likelihood
for (Tx in 1:nTx) {
  for (i in idx_2:(idx_2 + n_2[Tx] - 1)) {

     target += log_sum_exp(
    log(cf_2[Tx]) +
      exp_log_S(t_2[i], lambda_2_bg[i]),
    log1m(cf_2[Tx]) +
      exp_log_S(t_2[i], lambda_2_bg[i]) + weibull_log_S(t_2[i], shape_2, lambda_2[i]));

     target += d_2[i] * log_sum_exp(
    log(lambda_2_bg[i]),
    log1m(cf_2[Tx]) +
      weibull_lpdf(t_2[i] | shape_2, lambda_2[i]) - log(cf_2[Tx] + (1 - cf_2[Tx])*weibull_Surv(t_2[i], shape_2, lambda_2[i])));
  }

  idx_2 = idx_2 + n_2[Tx];
}
 idx_3 = 1;

// likelihood
for (Tx in 1:nTx) {
  for (i in idx_3:(idx_3 + n_3[Tx] - 1)) {

     target += log_sum_exp(
    log(cf_3[Tx]) +
      exp_log_S(t_3[i], lambda_3_bg[i]),
    log1m(cf_3[Tx]) +
      exp_log_S(t_3[i], lambda_3_bg[i]) + weibull_log_S(t_3[i], shape_3, lambda_3[i]));

     target += d_3[i] * log_sum_exp(
    log(lambda_3_bg[i]),
    log1m(cf_3[Tx]) +
      weibull_lpdf(t_3[i] | shape_3, lambda_3[i]) - log(cf_3[Tx] + (1 - cf_3[Tx])*weibull_Surv(t_3[i], shape_3, lambda_3[i])));
  }

  idx_3 = idx_3 + n_3[Tx];
}
 idx_4 = 1;

// likelihood
for (Tx in 1:nTx) {
  for (i in idx_4:(idx_4 + n_4[Tx] - 1)) {

     target += log_sum_exp(
    log(cf_4[Tx]) +
      exp_log_S(t_4[i], lambda_4_bg[i]),
    log1m(cf_4[Tx]) +
      exp_log_S(t_4[i], lambda_4_bg[i]) + weibull_log_S(t_4[i], shape_4, lambda_4[i]));

     target += d_4[i] * log_sum_exp(
    log(lambda_4_bg[i]),
    log1m(cf_4[Tx]) +
      weibull_lpdf(t_4[i] | shape_4, lambda_4[i]) - log(cf_4[Tx] + (1 - cf_4[Tx])*weibull_Surv(t_4[i], shape_4, lambda_4[i])));
  }

  idx_4 = idx_4 + n_4[Tx];
}
 idx_5 = 1;

// likelihood
for (Tx in 1:nTx) {
  for (i in idx_5:(idx_5 + n_5[Tx] - 1)) {

     target += log_sum_exp(
    log(cf_5[Tx]) +
      exp_log_S(t_5[i], lambda_5_bg[i]),
    log1m(cf_5[Tx]) +
      exp_log_S(t_5[i], lambda_5_bg[i]) + weibull_log_S(t_5[i], shape_5, lambda_5[i]));

     target += d_5[i] * log_sum_exp(
    log(lambda_5_bg[i]),
    log1m(cf_5[Tx]) +
      weibull_lpdf(t_5[i] | shape_5, lambda_5[i]) - log(cf_5[Tx] + (1 - cf_5[Tx])*weibull_Surv(t_5[i], shape_5, lambda_5[i])));
  }

  idx_5 = idx_5 + n_5[Tx];
}
 idx_6 = 1;

// likelihood
for (Tx in 1:nTx) {
  for (i in idx_6:(idx_6 + n_6[Tx] - 1)) {

     target += log_sum_exp(
    log(cf_6[Tx]) +
      exp_log_S(t_6[i], lambda_6_bg[i]),
    log1m(cf_6[Tx]) +
      exp_log_S(t_6[i], lambda_6_bg[i]) + weibull_log_S(t_6[i], shape_6, lambda_6[i]));

     target += d_6[i] * log_sum_exp(
    log(lambda_6_bg[i]),
    log1m(cf_6[Tx]) +
      weibull_lpdf(t_6[i] | shape_6, lambda_6[i]) - log(cf_6[Tx] + (1 - cf_6[Tx])*weibull_Surv(t_6[i], shape_6, lambda_6[i])));
  }

  idx_6 = idx_6 + n_6[Tx];
}
 idx_7 = 1;

// likelihood
for (Tx in 1:nTx) {
  for (i in idx_7:(idx_7 + n_7[Tx] - 1)) {

     target += log_sum_exp(
    log(cf_7[Tx]) +
      exp_log_S(t_7[i], lambda_7_bg[i]),
    log1m(cf_7[Tx]) +
      exp_log_S(t_7[i], lambda_7_bg[i]) + weibull_log_S(t_7[i], shape_7, lambda_7[i]));

     target += d_7[i] * log_sum_exp(
    log(lambda_7_bg[i]),
    log1m(cf_7[Tx]) +
      weibull_lpdf(t_7[i] | shape_7, lambda_7[i]) - log(cf_7[Tx] + (1 - cf_7[Tx])*weibull_Surv(t_7[i], shape_7, lambda_7[i])));
  }

  idx_7 = idx_7 + n_7[Tx];
}
 idx_8 = 1;

// likelihood
for (Tx in 1:nTx) {
  for (i in idx_8:(idx_8 + n_8[Tx] - 1)) {

     target += log_sum_exp(
    log(cf_8[Tx]) +
      exp_log_S(t_8[i], lambda_8_bg[i]),
    log1m(cf_8[Tx]) +
      exp_log_S(t_8[i], lambda_8_bg[i]) + weibull_log_S(t_8[i], shape_8, lambda_8[i]));

     target += d_8[i] * log_sum_exp(
    log(lambda_8_bg[i]),
    log1m(cf_8[Tx]) +
      weibull_lpdf(t_8[i] | shape_8, lambda_8[i]) - log(cf_8[Tx] + (1 - cf_8[Tx])*weibull_Surv(t_8[i], shape_8, lambda_8[i])));
  }

  idx_8 = idx_8 + n_8[Tx];
}
 idx_9 = 1;

// likelihood
for (Tx in 1:nTx) {
  for (i in idx_9:(idx_9 + n_9[Tx] - 1)) {

     target += log_sum_exp(
    log(cf_9[Tx]) +
      exp_log_S(t_9[i], lambda_9_bg[i]),
    log1m(cf_9[Tx]) +
      exp_log_S(t_9[i], lambda_9_bg[i]) + weibull_log_S(t_9[i], shape_9, lambda_9[i]));

     target += d_9[i] * log_sum_exp(
    log(lambda_9_bg[i]),
    log1m(cf_9[Tx]) +
      weibull_lpdf(t_9[i] | shape_9, lambda_9[i]) - log(cf_9[Tx] + (1 - cf_9[Tx])*weibull_Surv(t_9[i], shape_9, lambda_9[i])));
  }

  idx_9 = idx_9 + n_9[Tx];
}
 idx_10 = 1;

// likelihood
for (Tx in 1:nTx) {
  for (i in idx_10:(idx_10 + n_10[Tx] - 1)) {

     target += log_sum_exp(
    log(cf_10[Tx]) +
      exp_log_S(t_10[i], lambda_10_bg[i]),
    log1m(cf_10[Tx]) +
      exp_log_S(t_10[i], lambda_10_bg[i]) + weibull_log_S(t_10[i], shape_10, lambda_10[i]));

     target += d_10[i] * log_sum_exp(
    log(lambda_10_bg[i]),
    log1m(cf_10[Tx]) +
      weibull_lpdf(t_10[i] | shape_10, lambda_10[i]) - log(cf_10[Tx] + (1 - cf_10[Tx])*weibull_Surv(t_10[i], shape_10, lambda_10[i])));
  }

  idx_10 = idx_10 + n_10[Tx];
}
}

generated quantities {
real mean_bg;
 // real pbeta_bg;
 real log_lik = 0;
 vector[t_max] S_bg;
vector[t_max] S_1;
matrix[t_max, nTx] S_1_pred;
real mean_1;
int idx_1;
real log_lik_1;
// real pbeta_1 = normal_rng(mu_S_1[1], sigma_S_1[1]);


vector[t_max] S_2;
matrix[t_max, nTx] S_2_pred;
real mean_2;
int idx_2;
real log_lik_2;
// real pbeta_2 = normal_rng(mu_S_2[1], sigma_S_2[1]);


vector[t_max] S_3;
matrix[t_max, nTx] S_3_pred;
real mean_3;
int idx_3;
real log_lik_3;
// real pbeta_3 = normal_rng(mu_S_3[1], sigma_S_3[1]);


vector[t_max] S_4;
matrix[t_max, nTx] S_4_pred;
real mean_4;
int idx_4;
real log_lik_4;
// real pbeta_4 = normal_rng(mu_S_4[1], sigma_S_4[1]);


vector[t_max] S_5;
matrix[t_max, nTx] S_5_pred;
real mean_5;
int idx_5;
real log_lik_5;
// real pbeta_5 = normal_rng(mu_S_5[1], sigma_S_5[1]);


vector[t_max] S_6;
matrix[t_max, nTx] S_6_pred;
real mean_6;
int idx_6;
real log_lik_6;
// real pbeta_6 = normal_rng(mu_S_6[1], sigma_S_6[1]);


vector[t_max] S_7;
matrix[t_max, nTx] S_7_pred;
real mean_7;
int idx_7;
real log_lik_7;
// real pbeta_7 = normal_rng(mu_S_7[1], sigma_S_7[1]);


vector[t_max] S_8;
matrix[t_max, nTx] S_8_pred;
real mean_8;
int idx_8;
real log_lik_8;
// real pbeta_8 = normal_rng(mu_S_8[1], sigma_S_8[1]);


vector[t_max] S_9;
matrix[t_max, nTx] S_9_pred;
real mean_9;
int idx_9;
real log_lik_9;
// real pbeta_9 = normal_rng(mu_S_9[1], sigma_S_9[1]);


vector[t_max] S_10;
matrix[t_max, nTx] S_10_pred;
real mean_10;
int idx_10;
real log_lik_10;
// real pbeta_10 = normal_rng(mu_S_10[1], sigma_S_10[1]);

real pshape_1 = gamma_rng(a_shape_1, b_shape_1);
 real pshape_2 = gamma_rng(a_shape_2, b_shape_2);
 real pshape_3 = gamma_rng(a_shape_3, b_shape_3);
 real pshape_4 = gamma_rng(a_shape_4, b_shape_4);
 real pshape_5 = gamma_rng(a_shape_5, b_shape_5);
 real pshape_6 = gamma_rng(a_shape_6, b_shape_6);
 real pshape_7 = gamma_rng(a_shape_7, b_shape_7);
 real pshape_8 = gamma_rng(a_shape_8, b_shape_8);
 real pshape_9 = gamma_rng(a_shape_9, b_shape_9);
 real pshape_10 = gamma_rng(a_shape_10, b_shape_10);mean_1 = exp(beta_1[1]);
 mean_2 = exp(beta_2[1]);
 mean_3 = exp(beta_3[1]);
 mean_4 = exp(beta_4[1]);
 mean_5 = exp(beta_5[1]);
 mean_6 = exp(beta_6[1]);
 mean_7 = exp(beta_7[1]);
 mean_8 = exp(beta_8[1]);
 mean_9 = exp(beta_9[1]);
 mean_10 = exp(beta_10[1]);// background rate
if (bg_model == 1) {
	mean_bg = exp(beta_bg[1]);
} else {
// mean_bg = 0.001;
mean_bg = mean(h_bg_1);
mean_bg = mean(h_bg_2);
mean_bg = mean(h_bg_3);
mean_bg = mean(h_bg_4);
mean_bg = mean(h_bg_5);
mean_bg = mean(h_bg_6);
mean_bg = mean(h_bg_7);
mean_bg = mean(h_bg_8);
mean_bg = mean(h_bg_9);
mean_bg = mean(h_bg_10);
}
// posterior mean checks
for (j in 1:nTx) {
  for (i in 1:t_max) {
    S_bg[i] = exp_Surv(i, mean_bg);
    S_1[i] = exp_weibull_Surv(i, shape_1, mean_1, mean_bg);
    S_1_pred[i, j] = cf_1[j]*S_bg[i] + (1 - cf_1[j])*S_1[i];
  }
}

 // posterior mean checks
for (j in 1:nTx) {
  for (i in 1:t_max) {
    S_bg[i] = exp_Surv(i, mean_bg);
    S_2[i] = exp_weibull_Surv(i, shape_2, mean_2, mean_bg);
    S_2_pred[i, j] = cf_2[j]*S_bg[i] + (1 - cf_2[j])*S_2[i];
  }
}

 // posterior mean checks
for (j in 1:nTx) {
  for (i in 1:t_max) {
    S_bg[i] = exp_Surv(i, mean_bg);
    S_3[i] = exp_weibull_Surv(i, shape_3, mean_3, mean_bg);
    S_3_pred[i, j] = cf_3[j]*S_bg[i] + (1 - cf_3[j])*S_3[i];
  }
}

 // posterior mean checks
for (j in 1:nTx) {
  for (i in 1:t_max) {
    S_bg[i] = exp_Surv(i, mean_bg);
    S_4[i] = exp_weibull_Surv(i, shape_4, mean_4, mean_bg);
    S_4_pred[i, j] = cf_4[j]*S_bg[i] + (1 - cf_4[j])*S_4[i];
  }
}

 // posterior mean checks
for (j in 1:nTx) {
  for (i in 1:t_max) {
    S_bg[i] = exp_Surv(i, mean_bg);
    S_5[i] = exp_weibull_Surv(i, shape_5, mean_5, mean_bg);
    S_5_pred[i, j] = cf_5[j]*S_bg[i] + (1 - cf_5[j])*S_5[i];
  }
}

 // posterior mean checks
for (j in 1:nTx) {
  for (i in 1:t_max) {
    S_bg[i] = exp_Surv(i, mean_bg);
    S_6[i] = exp_weibull_Surv(i, shape_6, mean_6, mean_bg);
    S_6_pred[i, j] = cf_6[j]*S_bg[i] + (1 - cf_6[j])*S_6[i];
  }
}

 // posterior mean checks
for (j in 1:nTx) {
  for (i in 1:t_max) {
    S_bg[i] = exp_Surv(i, mean_bg);
    S_7[i] = exp_weibull_Surv(i, shape_7, mean_7, mean_bg);
    S_7_pred[i, j] = cf_7[j]*S_bg[i] + (1 - cf_7[j])*S_7[i];
  }
}

 // posterior mean checks
for (j in 1:nTx) {
  for (i in 1:t_max) {
    S_bg[i] = exp_Surv(i, mean_bg);
    S_8[i] = exp_weibull_Surv(i, shape_8, mean_8, mean_bg);
    S_8_pred[i, j] = cf_8[j]*S_bg[i] + (1 - cf_8[j])*S_8[i];
  }
}

 // posterior mean checks
for (j in 1:nTx) {
  for (i in 1:t_max) {
    S_bg[i] = exp_Surv(i, mean_bg);
    S_9[i] = exp_weibull_Surv(i, shape_9, mean_9, mean_bg);
    S_9_pred[i, j] = cf_9[j]*S_bg[i] + (1 - cf_9[j])*S_9[i];
  }
}

 // posterior mean checks
for (j in 1:nTx) {
  for (i in 1:t_max) {
    S_bg[i] = exp_Surv(i, mean_bg);
    S_10[i] = exp_weibull_Surv(i, shape_10, mean_10, mean_bg);
    S_10_pred[i, j] = cf_10[j]*S_bg[i] + (1 - cf_10[j])*S_10[i];
  }
}
// prior mean checks
// pmean_pfs = exp(pbeta_pfs);
// pmean_bg = exp(pbeta_bg);

// for (i in 1:t_max) {
//  pS_bg[i] = exp_Surv(i, pmean_bg);
//  pS_1[i] = exp_weibull_Surv(i, shape_1, lambda_1[j], pmean_bg);
//  S_1_prior[i] = pmean_cf_1*pS_bg[i] + (1 - pmean_cf_1)*pS_1[i,j];
// }

 // prior mean checks
// pmean_pfs = exp(pbeta_pfs);
// pmean_bg = exp(pbeta_bg);

// for (i in 1:t_max) {
//  pS_bg[i] = exp_Surv(i, pmean_bg);
//  pS_2[i] = exp_weibull_Surv(i, shape_2, lambda_2[j], pmean_bg);
//  S_2_prior[i] = pmean_cf_2*pS_bg[i] + (1 - pmean_cf_2)*pS_2[i,j];
// }

 // prior mean checks
// pmean_pfs = exp(pbeta_pfs);
// pmean_bg = exp(pbeta_bg);

// for (i in 1:t_max) {
//  pS_bg[i] = exp_Surv(i, pmean_bg);
//  pS_3[i] = exp_weibull_Surv(i, shape_3, lambda_3[j], pmean_bg);
//  S_3_prior[i] = pmean_cf_3*pS_bg[i] + (1 - pmean_cf_3)*pS_3[i,j];
// }

 // prior mean checks
// pmean_pfs = exp(pbeta_pfs);
// pmean_bg = exp(pbeta_bg);

// for (i in 1:t_max) {
//  pS_bg[i] = exp_Surv(i, pmean_bg);
//  pS_4[i] = exp_weibull_Surv(i, shape_4, lambda_4[j], pmean_bg);
//  S_4_prior[i] = pmean_cf_4*pS_bg[i] + (1 - pmean_cf_4)*pS_4[i,j];
// }

 // prior mean checks
// pmean_pfs = exp(pbeta_pfs);
// pmean_bg = exp(pbeta_bg);

// for (i in 1:t_max) {
//  pS_bg[i] = exp_Surv(i, pmean_bg);
//  pS_5[i] = exp_weibull_Surv(i, shape_5, lambda_5[j], pmean_bg);
//  S_5_prior[i] = pmean_cf_5*pS_bg[i] + (1 - pmean_cf_5)*pS_5[i,j];
// }

 // prior mean checks
// pmean_pfs = exp(pbeta_pfs);
// pmean_bg = exp(pbeta_bg);

// for (i in 1:t_max) {
//  pS_bg[i] = exp_Surv(i, pmean_bg);
//  pS_6[i] = exp_weibull_Surv(i, shape_6, lambda_6[j], pmean_bg);
//  S_6_prior[i] = pmean_cf_6*pS_bg[i] + (1 - pmean_cf_6)*pS_6[i,j];
// }

 // prior mean checks
// pmean_pfs = exp(pbeta_pfs);
// pmean_bg = exp(pbeta_bg);

// for (i in 1:t_max) {
//  pS_bg[i] = exp_Surv(i, pmean_bg);
//  pS_7[i] = exp_weibull_Surv(i, shape_7, lambda_7[j], pmean_bg);
//  S_7_prior[i] = pmean_cf_7*pS_bg[i] + (1 - pmean_cf_7)*pS_7[i,j];
// }

 // prior mean checks
// pmean_pfs = exp(pbeta_pfs);
// pmean_bg = exp(pbeta_bg);

// for (i in 1:t_max) {
//  pS_bg[i] = exp_Surv(i, pmean_bg);
//  pS_8[i] = exp_weibull_Surv(i, shape_8, lambda_8[j], pmean_bg);
//  S_8_prior[i] = pmean_cf_8*pS_bg[i] + (1 - pmean_cf_8)*pS_8[i,j];
// }

 // prior mean checks
// pmean_pfs = exp(pbeta_pfs);
// pmean_bg = exp(pbeta_bg);

// for (i in 1:t_max) {
//  pS_bg[i] = exp_Surv(i, pmean_bg);
//  pS_9[i] = exp_weibull_Surv(i, shape_9, lambda_9[j], pmean_bg);
//  S_9_prior[i] = pmean_cf_9*pS_bg[i] + (1 - pmean_cf_9)*pS_9[i,j];
// }

 // prior mean checks
// pmean_pfs = exp(pbeta_pfs);
// pmean_bg = exp(pbeta_bg);

// for (i in 1:t_max) {
//  pS_bg[i] = exp_Surv(i, pmean_bg);
//  pS_10[i] = exp_weibull_Surv(i, shape_10, lambda_10[j], pmean_bg);
//  S_10_prior[i] = pmean_cf_10*pS_bg[i] + (1 - pmean_cf_10)*pS_10[i,j];
// }
// likelihood
  idx_1 = 1;

  for (Tx in 1:nTx) {
    for (i in idx_1:(idx_1 + n_1[Tx] - 1)) {
      log_lik_1 += log_sum_exp(
      log(cf_1[Tx]) +
      surv_exp_lpdf(t_1[i] | d_1[i], lambda_1_bg[i]),
      log1m(cf_1[Tx]) +
      joint_exp_weibull_lpdf(t_1[i] | d_1[i], shape_1, lambda_1[i], lambda_1_bg[i]));
    }

    idx_1 = idx_1 + n_1[Tx];
  }
  log_lik = log_lik + log_lik_1;

 // likelihood
  idx_2 = 1;

  for (Tx in 1:nTx) {
    for (i in idx_2:(idx_2 + n_2[Tx] - 1)) {
      log_lik_2 += log_sum_exp(
      log(cf_2[Tx]) +
      surv_exp_lpdf(t_2[i] | d_2[i], lambda_2_bg[i]),
      log1m(cf_2[Tx]) +
      joint_exp_weibull_lpdf(t_2[i] | d_2[i], shape_2, lambda_2[i], lambda_2_bg[i]));
    }

    idx_2 = idx_2 + n_2[Tx];
  }
  log_lik = log_lik + log_lik_2;

 // likelihood
  idx_3 = 1;

  for (Tx in 1:nTx) {
    for (i in idx_3:(idx_3 + n_3[Tx] - 1)) {
      log_lik_3 += log_sum_exp(
      log(cf_3[Tx]) +
      surv_exp_lpdf(t_3[i] | d_3[i], lambda_3_bg[i]),
      log1m(cf_3[Tx]) +
      joint_exp_weibull_lpdf(t_3[i] | d_3[i], shape_3, lambda_3[i], lambda_3_bg[i]));
    }

    idx_3 = idx_3 + n_3[Tx];
  }
  log_lik = log_lik + log_lik_3;

 // likelihood
  idx_4 = 1;

  for (Tx in 1:nTx) {
    for (i in idx_4:(idx_4 + n_4[Tx] - 1)) {
      log_lik_4 += log_sum_exp(
      log(cf_4[Tx]) +
      surv_exp_lpdf(t_4[i] | d_4[i], lambda_4_bg[i]),
      log1m(cf_4[Tx]) +
      joint_exp_weibull_lpdf(t_4[i] | d_4[i], shape_4, lambda_4[i], lambda_4_bg[i]));
    }

    idx_4 = idx_4 + n_4[Tx];
  }
  log_lik = log_lik + log_lik_4;

 // likelihood
  idx_5 = 1;

  for (Tx in 1:nTx) {
    for (i in idx_5:(idx_5 + n_5[Tx] - 1)) {
      log_lik_5 += log_sum_exp(
      log(cf_5[Tx]) +
      surv_exp_lpdf(t_5[i] | d_5[i], lambda_5_bg[i]),
      log1m(cf_5[Tx]) +
      joint_exp_weibull_lpdf(t_5[i] | d_5[i], shape_5, lambda_5[i], lambda_5_bg[i]));
    }

    idx_5 = idx_5 + n_5[Tx];
  }
  log_lik = log_lik + log_lik_5;

 // likelihood
  idx_6 = 1;

  for (Tx in 1:nTx) {
    for (i in idx_6:(idx_6 + n_6[Tx] - 1)) {
      log_lik_6 += log_sum_exp(
      log(cf_6[Tx]) +
      surv_exp_lpdf(t_6[i] | d_6[i], lambda_6_bg[i]),
      log1m(cf_6[Tx]) +
      joint_exp_weibull_lpdf(t_6[i] | d_6[i], shape_6, lambda_6[i], lambda_6_bg[i]));
    }

    idx_6 = idx_6 + n_6[Tx];
  }
  log_lik = log_lik + log_lik_6;

 // likelihood
  idx_7 = 1;

  for (Tx in 1:nTx) {
    for (i in idx_7:(idx_7 + n_7[Tx] - 1)) {
      log_lik_7 += log_sum_exp(
      log(cf_7[Tx]) +
      surv_exp_lpdf(t_7[i] | d_7[i], lambda_7_bg[i]),
      log1m(cf_7[Tx]) +
      joint_exp_weibull_lpdf(t_7[i] | d_7[i], shape_7, lambda_7[i], lambda_7_bg[i]));
    }

    idx_7 = idx_7 + n_7[Tx];
  }
  log_lik = log_lik + log_lik_7;

 // likelihood
  idx_8 = 1;

  for (Tx in 1:nTx) {
    for (i in idx_8:(idx_8 + n_8[Tx] - 1)) {
      log_lik_8 += log_sum_exp(
      log(cf_8[Tx]) +
      surv_exp_lpdf(t_8[i] | d_8[i], lambda_8_bg[i]),
      log1m(cf_8[Tx]) +
      joint_exp_weibull_lpdf(t_8[i] | d_8[i], shape_8, lambda_8[i], lambda_8_bg[i]));
    }

    idx_8 = idx_8 + n_8[Tx];
  }
  log_lik = log_lik + log_lik_8;

 // likelihood
  idx_9 = 1;

  for (Tx in 1:nTx) {
    for (i in idx_9:(idx_9 + n_9[Tx] - 1)) {
      log_lik_9 += log_sum_exp(
      log(cf_9[Tx]) +
      surv_exp_lpdf(t_9[i] | d_9[i], lambda_9_bg[i]),
      log1m(cf_9[Tx]) +
      joint_exp_weibull_lpdf(t_9[i] | d_9[i], shape_9, lambda_9[i], lambda_9_bg[i]));
    }

    idx_9 = idx_9 + n_9[Tx];
  }
  log_lik = log_lik + log_lik_9;

 // likelihood
  idx_10 = 1;

  for (Tx in 1:nTx) {
    for (i in idx_10:(idx_10 + n_10[Tx] - 1)) {
      log_lik_10 += log_sum_exp(
      log(cf_10[Tx]) +
      surv_exp_lpdf(t_10[i] | d_10[i], lambda_10_bg[i]),
      log1m(cf_10[Tx]) +
      joint_exp_weibull_lpdf(t_10[i] | d_10[i], shape_10, lambda_10[i], lambda_10_bg[i]));
    }

    idx_10 = idx_10 + n_10[Tx];
  }
  log_lik = log_lik + log_lik_10;
}

