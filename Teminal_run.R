library(dplyr)
library(rstan)

#read data, mask identifies missingness
mask <- as.matrix(read.csv(mask))             #replace with mask file
final_data <- as.data.frame(read.csv(dat.sim)) #replace with the simulated dataset
final_data$x2_c <- final_data$x2 - mean(final_data$x2, na.rm=TRUE)

time_points <- as.matrix(final_data[, grep("time_", names(final_data))])
qol_values <- as.matrix(final_data[, grep("qol_", names(final_data))])

dat_death <- subset(final_data, status==1)
observed_times_death <- as.matrix(dat_death[, grep("time_", names(dat_death))])

#Calculate backward times for each subject at each time point
backward_time_list <- lapply(1:nrow(observed_times_death), function(i) {
  backward_times <- dat_death$observed_time[i] - observed_times_death[i, ]
  backward_times[backward_times < 0] <- NA
  return(backward_times)
})

backward_time_vector <- unlist(backward_time_list)
backward_time_vector <- backward_time_vector[!is.na(backward_time_vector)]

#define spline parameters
num_knots <- 3
degree <- 1
knots <- unname(quantile(backward_time_vector, probs = seq(0, 1, length.out = num_knots + 2)[-c(1, num_knots + 2)]))

time_points[is.na(time_points)] <- 0
qol_values[is.na(qol_values)] <- 0

#prepare data for Stan
stan_data <- list(
  N = nrow(final_data),
  K = length(unique(final_data$cluster)),
  T = length(grep("time_", names(final_data))),
  cluster = final_data$cluster,
  treatment = final_data$treatment,
  time_points = time_points,
  x1 = final_data$x1,
  x2 = final_data$x2_c,
  Y = qol_values,
  MASK = mask,
  survival_time = final_data$observed_time,
  status = final_data$status,
  num_knots = num_knots,
  knots = knots,
  spline_degree = degree
)

#compile the Stan model
stan_model_code <- "
// Function for spline using De Boor's algorithm
functions {
  vector build_b_spline(real[] t, real[] ext_knots, int ind, int order) {
    vector[size(t)] b_spline;
    vector[size(t)] w1 = rep_vector(0, size(t));
    vector[size(t)] w2 = rep_vector(0, size(t));
    
    if (order == 1) {
      for (i in 1:size(t)) {
        b_spline[i] = (ext_knots[ind] <= t[i]) && (t[i] <= ext_knots[ind + 1]);
      }
    } else {
      if (ext_knots[ind] != ext_knots[ind + order - 1]) {
        w1 = (to_vector(t) - rep_vector(ext_knots[ind], size(t))) / (ext_knots[ind + order - 1] - ext_knots[ind]);
      }
      
      if (ext_knots[ind + 1] != ext_knots[ind + order]) {
        w2 = 1 - (to_vector(t) - rep_vector(ext_knots[ind + 1], size(t))) / (ext_knots[ind + order] - ext_knots[ind + 1]);
      }
      
      b_spline = w1 .* build_b_spline(t, ext_knots, ind, order - 1) + w2 .* build_b_spline(t, ext_knots, ind + 1, order - 1);
    }
    
    return b_spline;
  }
}

data {
  int<lower=0> N;                   // total number of subjects
  int<lower=0> K;                   // number of clusters
  int<lower=0> T;                   // number of time points
  int<lower=1,upper=K> cluster[N];  // cluster indicator
  vector[N] x1;                     // binary covariate
  vector[N] x2;                     // continuous covariate
  vector[N] treatment;
  matrix[N, T] time_points;
  matrix[N, T] Y;                   // outcome measurements
  matrix[N, T] MASK;                // mask for missing values
  vector[N] survival_time;          // observed survival or censoring times
  int<lower=0,upper=1> status[N];   // censoring indicator
  int num_knots;                    // number of internal knots for the splines
  vector[num_knots] knots;          // sequence of internal knots
  int spline_degree;                // degree of the splines (order - 1)
}

transformed data {
  matrix[N, T] X1 = rep_matrix(x1, T);
  matrix[N, T] X2 = rep_matrix(x2, T);
  matrix[N, T] TRT = rep_matrix(treatment, T);
  
  int num_basis = num_knots + spline_degree + 1;  // total number of B-splines
  
  int non_missing_count = 0;
  
  for (i in 1:N) {
    for (t in 1:T) {
      if (MASK[i, t] > 0) {
        non_missing_count += 1;
      }
    }
  }
  
  vector[non_missing_count] Y_non_missing;  
  int non_missing_indices[non_missing_count, 2];  
  int index = 1;  

  for (i in 1:N) {
    for (t in 1:T) {
      if (MASK[i, t] > 0) {
        Y_non_missing[index] = Y[i, t];  
        non_missing_indices[index, 1] = i;  
        non_missing_indices[index, 2] = t;  
        index += 1;
      }
    }
  }
}

parameters {
  real alpha1;
  real alpha2;
  real phi0;
  real phi1;
  real phi2;

  real v1;
  real v2;

  real<lower=0> sigma_b;
  real<lower=0> sigma_u;
  real<lower=0> sigma_e;
  real<lower=0> gamma;


  row_vector[num_basis] a_backward;  
  row_vector[num_basis] a_treatment;
  real<lower=0> sigma_rw2_backward; 
  real<lower=0> sigma_rw2_treatment; 
  
  
  vector[N] z_b;
  vector[K] z_u;
  vector<lower=0, upper=1>[N] U;
}

transformed parameters {
  vector[N] lambda;
  vector[N] death_time;
  vector[N] b_i = z_b * sigma_b;
  vector[K] u_i = z_u * sigma_u;
  b_i = b_i - mean(b_i);
  u_i = u_i - mean(u_i);

  vector[non_missing_count] backward_times_non_missing;  
  matrix[num_basis, non_missing_count] B;  
  vector[non_missing_count] spline_contribution_backward;
  vector[non_missing_count] spline_contribution_treatment;

  matrix[N, T] BI = rep_matrix(b_i, T);
  matrix[N, T] UI;
  matrix[N, T] MU;

  lambda =  exp(phi0 + phi1 * x1 + phi2 * x2 + v1 * b_i + v2 * u_i[cluster]);

  for (i in 1:N) {
    if (status[i] == 1) {
      death_time[i] = survival_time[i]; 
    } else {
      death_time[i] = pow(-log(exp(-lambda[i] * pow(survival_time[i], gamma)) - U[i] * exp(-lambda[i] * pow(survival_time[i], gamma))) / lambda[i], 1 / gamma);
    }
  }

  vector[spline_degree + num_knots + 1] ext_knots_temp;
  vector[2 * spline_degree + num_knots + 2] ext_knots;

  for (ind in 1:non_missing_count) {
    int subj = non_missing_indices[ind, 1];
    int time_idx = non_missing_indices[ind, 2];
    backward_times_non_missing[ind] = death_time[subj] - time_points[subj, time_idx];
  }
  
  ext_knots_temp = append_row(rep_vector(min(backward_times_non_missing), spline_degree+1), knots);
  ext_knots = append_row(ext_knots_temp, rep_vector(max(backward_times_non_missing), spline_degree+1));
  
  for (ind in 1:num_basis) {
    B[ind, :] = to_row_vector(build_b_spline(to_array_1d(backward_times_non_missing), to_array_1d(ext_knots), ind, spline_degree + 1));
  }
  
  spline_contribution_backward = to_vector(a_backward * B);
  spline_contribution_treatment = to_vector(a_treatment * B);
  
  MU = rep_matrix(0, N, T);
  UI = rep_matrix(u_i[cluster], T);
  
  for (ind in 1:non_missing_count) {
    int subj = non_missing_indices[ind, 1];
    int time_idx = non_missing_indices[ind, 2];
    MU[subj, time_idx] = X1[subj, time_idx] * alpha1 + X2[subj, time_idx] * alpha2 + spline_contribution_backward[ind] + TRT[subj, time_idx] * spline_contribution_treatment[ind] + BI[subj, time_idx] + UI[subj, time_idx];
  }
}

model {
  // define priors
  alpha1 ~ normal(0, 10);
  alpha2 ~ normal(0, 10);
  phi0 ~ normal(0, 10);
  phi1 ~ normal(0, 10);
  phi2 ~ normal(0, 10);

  v1 ~ normal(0, 1);
  v2 ~ normal(0, 1);

  sigma_b ~ normal(0, 10);
  sigma_u ~ normal(0, 10);
  sigma_e ~ normal(0, 10);
  gamma ~ gamma(0.5, 0.5);
  
  a_backward[1] ~ normal(0, 10);
  a_backward[2] ~ normal(0, 10);
  a_treatment[1] ~ normal(0, 10);
  a_treatment[2] ~ normal(0, 10);
  sigma_rw2_backward ~ inv_gamma(0.01, 0.01);
  sigma_rw2_treatment ~ inv_gamma(0.01, 0.01);


  z_b ~ normal(0, 1);
  z_u ~ normal(0, 1);
  U ~ uniform(0, 1);
  
  for (i in 3:num_basis) {
    a_backward[i] ~ normal(2 * a_backward[i - 1] - a_backward[i - 2], sqrt(sigma_rw2_backward));
    a_treatment[i] ~ normal(2 * a_treatment[i - 1] - a_treatment[i - 2], sqrt(sigma_rw2_treatment));
  }
  
  // Longitudinal model
   for (i in 1:N) {
    for (t in 1:T) {
      if (MASK[i, t] > 0) {
        Y[i, t] ~ normal(MU[i, t], sigma_e);
      }
    }
   }
  
  // Survival likelihood
  for (i in 1:N) {
    if (status[i] == 1) {
      // Observed death times
      target += weibull_lpdf(survival_time[i] | gamma, pow(lambda[i], -1 / gamma));
    } else {
      // Censored observations
      target += weibull_lccdf(survival_time[i] | gamma, pow(lambda[i], -1 / gamma));
    }
  }
}
"

stan_model <- stan_model(model_code = stan_model_code)

init_fn <- function() {
  list(alpha1 = 0, alpha2 = 0, phi0 = 0, phi1 = 0, phi2 = 0, v1 = 0, v2 = 0, sigma_u = 1, sigma_b = 1, sigma_e = 1, gamma = 1,
       z_b = rnorm(nrow(final_data), 0, 1),
       z_u = rnorm(length(unique(final_data$cluster)), 0, 1),  
       U = runif(nrow(final_data), 0, 1), 
       a_backward = rnorm(length(knots)+degree+1, 0, 10),  
       a_treatment = rnorm(length(knots)+degree+1, 0, 10),
       sigma_rw2_backward = 1,
       sigma_rw2_treatment = 1)
}

# Compile and sample from the Stan model
fit <- sampling(stan_model, data = stan_data, init = init_fn, iter = 2500, warmup = 1000, chains = 2, control = list(adapt_delta = 0.99, max_treedepth = 15), cores = 2, refresh=100)

rest_model_sum <- as.data.frame(summary(fit)$summary)
rest_model_sum$parameter <- row.names(rest_model_sum)


