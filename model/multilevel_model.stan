data {
  int<lower=1> N;
  vector[N] y;
  vector<lower=0>[N] se;
  vector[N] time_cen_yr;
  vector[N] time_sqrt;
  int<lower=1> n_site;
  int<lower=1> n_length;
  int<lower=1> n_graft;
  int<lower=1> n_study;
  array[N] int<lower=1, upper=n_site> site_id;
  array[N] int<lower=1, upper=n_length> length_id;
  array[N] int<lower=1, upper=n_graft> graft_id;
  array[N] int<lower=1, upper=n_study> study_id;
}

parameters {
  real beta0;
  real beta_time_cen;
  real beta_time_sqrt;

  vector[n_site] site_raw;
  vector[n_length] length_raw;
  vector[n_graft] graft_raw;
  vector[n_study] study_raw;

  vector[n_site] site_time_cen_raw;
  vector[n_length] length_time_cen_raw;
  vector[n_graft] graft_time_cen_raw;

  vector[n_site] site_time_sqrt_raw;
  vector[n_length] length_time_sqrt_raw;
  vector[n_graft] graft_time_sqrt_raw;

  real<lower=0> sigma_site;
  real<lower=0> sigma_length;
  real<lower=0> sigma_graft;
  real<lower=0> sigma_study;

  real<lower=0> sigma_site_time_cen;
  real<lower=0> sigma_length_time_cen;
  real<lower=0> sigma_graft_time_cen;

  real<lower=0> sigma_site_time_sqrt;
  real<lower=0> sigma_length_time_sqrt;
  real<lower=0> sigma_graft_time_sqrt;
}

transformed parameters {
  vector[n_site] site_effect          = sigma_site * site_raw;
  vector[n_length] length_effect      = sigma_length * length_raw;
  vector[n_graft] graft_effect        = sigma_graft * graft_raw;
  vector[n_study] study_effect        = sigma_study * study_raw;

  vector[n_site] site_time_cen        = sigma_site_time_cen * site_time_cen_raw;
  vector[n_length] length_time_cen    = sigma_length_time_cen * length_time_cen_raw;
  vector[n_graft] graft_time_cen      = sigma_graft_time_cen * graft_time_cen_raw;

  vector[n_site] site_time_sqrt       = sigma_site_time_sqrt * site_time_sqrt_raw;
  vector[n_length] length_time_sqrt   = sigma_length_time_sqrt * length_time_sqrt_raw;
  vector[n_graft] graft_time_sqrt     = sigma_graft_time_sqrt * graft_time_sqrt_raw;
}

model {
  // Priors
  beta0          ~ normal(0, 5);
  beta_time_cen  ~ normal(0, 1);
  beta_time_sqrt ~ normal(0, 1);

  site_raw             ~ normal(0, 1);
  length_raw           ~ normal(0, 1);
  graft_raw            ~ normal(0, 1);
  study_raw            ~ normal(0, 1);

  site_time_cen_raw      ~ normal(0, 1);
  length_time_cen_raw    ~ normal(0, 1);
  graft_time_cen_raw     ~ normal(0, 1);

  site_time_sqrt_raw     ~ normal(0, 1);
  length_time_sqrt_raw   ~ normal(0, 1);
  graft_time_sqrt_raw    ~ normal(0, 1);

  sigma_site             ~ normal(0, 1);
  sigma_length           ~ normal(0, 1);
  sigma_graft            ~ normal(0, 1);
  sigma_study            ~ normal(0, 1);

  sigma_site_time_cen    ~ normal(0, 0.5);
  sigma_length_time_cen  ~ normal(0, 0.5);
  sigma_graft_time_cen   ~ normal(0, 0.5);

  sigma_site_time_sqrt   ~ normal(0, 0.5);
  sigma_length_time_sqrt ~ normal(0, 0.5);
  sigma_graft_time_sqrt  ~ normal(0, 0.5);

  // Likelihood
  for (n in 1:N) {
    real mu_n =
      beta0
      + site_effect[site_id[n]]
      + length_effect[length_id[n]]
      + graft_effect[graft_id[n]]
      + study_effect[study_id[n]]
      + (beta_time_cen
         + site_time_cen[site_id[n]]
         + length_time_cen[length_id[n]]
         + graft_time_cen[graft_id[n]]) * time_cen_yr[n]
      + (beta_time_sqrt
         + site_time_sqrt[site_id[n]]
         + length_time_sqrt[length_id[n]]
         + graft_time_sqrt[graft_id[n]]) * time_sqrt[n];

    y[n] ~ normal(mu_n, se[n]);
  }
}

generated quantities {
  vector[N] mu;
  vector[N] log_lik;
  for (n in 1:N) {
    mu[n] =
      beta0
      + site_effect[site_id[n]]
      + length_effect[length_id[n]]
      + graft_effect[graft_id[n]]
      + study_effect[study_id[n]]
      + (beta_time_cen
         + site_time_cen[site_id[n]]
         + length_time_cen[length_id[n]]
         + graft_time_cen[graft_id[n]]) * time_cen_yr[n]
      + (beta_time_sqrt
         + site_time_sqrt[site_id[n]]
         + length_time_sqrt[length_id[n]]
         + graft_time_sqrt[graft_id[n]]) * time_sqrt[n];

    log_lik[n] = normal_lpdf(y[n] | mu[n], se[n]);
  }
}
