data {
  int<lower=1> N;
  int<lower=1> p;
  matrix[N, p] X;
  vector[N] y;
  vector[2] coords[N];
}
transformed data {
  real delta = 1e-9;
}
parameters {
  real<lower=0> sigma;
  real<lower=0> phi;
  vector[p] beta;
}
transformed parameters {
  vector[N] mu;
  matrix[N, N] L_Sigma;
  
  {
    matrix[N, N] Sigma;
    // Sigma
    Sigma = cov_exp_quad(coords, sigma, phi);
    
    // diagonal elements
    for (n in 1:N)
      Sigma[n, n] = Sigma[n, n] + delta;  

    L_Sigma = cholesky_decompose(Sigma); 
  }
  
  // mean
  mu = X * beta;  
  
}
model {
  sigma ~ inv_gamma(1, 1);
  phi ~ inv_gamma(4, 11.2);
  beta ~ normal(0, 100);
  
  y ~ multi_normal_cholesky(mu, L_Sigma);
  
}

