data {
  int<lower=1> N;
  vector[2] coords[N];
  vector[N] y;
}
transformed data {
  real delta = 1e-9;
}
parameters {
  real<lower=0> rho;
  real<lower=0> sigma;
  real<lower=0> tau;
  vector[N] eta;
}
transformed parameters {
  vector[N] f;
  
  {
    matrix[N, N] L_K;
    matrix[N, N] K = cov_exp_quad(coords, sigma, rho);
    
    // diagonal elements
    for (n in 1:N)
      K[n, n] = K[n, n] + delta;
    
    L_K = cholesky_decompose(K);
    f = L_K * eta;
  }  
  
  
}
model {
  rho ~ inv_gamma(6.7, 28.4);
  sigma ~ std_normal();
  tau ~ std_normal();
  eta ~ std_normal();
  
  y ~ normal(f, tau);
}

