data {
  int<lower=1> N;
  int<lower=1> p;
  matrix[N, N] Omega;
  matrix[N, p] X;
  vector[N] y;
}
transformed data {
  real delta = 1e-9;
}
parameters {
  real<lower=0> sigma;
  vector[p] beta;
}
transformed parameters {
  real sigma2;
  vector[N] mu;
  matrix[N, N] Sigma;
  matrix[N, N] L_Sigma;

  // sigma
  sigma2 = sigma * sigma;
  
  // Sigma
  Sigma = sigma2 * Omega;
  
  // diagonal elements
  for (n in 1:N)
    Sigma[n, n] = Sigma[n, n] + delta;   
  
  L_Sigma = cholesky_decompose(Sigma);
  
  // mean
  mu = X * beta;  
  
}
model {
  sigma ~ inv_gamma(.1, .1);
  beta ~ normal(0, 100);
  
  y ~ multi_normal_cholesky(mu, L_Sigma);

}

