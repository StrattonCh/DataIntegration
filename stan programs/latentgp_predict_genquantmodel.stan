data {
  int<lower=1> N;
  int<lower=1> n;
  int<lower=1> m;
  vector[2] coords[n];
  vector[2] coords_all[N];
  vector[n] y;
  vector[m] z;
}
transformed data {
  real delta = 1e-4;
}
parameters {
  real<lower=0> rho;
  real<lower=0> sigma;
  real<lower=0> tau;
  vector[n] eta;
  vector[m] eta2;
  real beta;
}
transformed parameters {
  vector[n] f;
  matrix[n, n] Sigma11;
  
  {
    matrix[n, n] L_Sigma11;
    Sigma11 = cov_exp_quad(coords, sigma, rho);
        
    // diagonal elements
    for (k in 1:n)
      Sigma11[k, k] = Sigma11[k, k] + delta;
    
    L_Sigma11 = cholesky_decompose(Sigma11);
    f = L_Sigma11 * eta;
  }  

  // storage
  vector[m] f_pred;
  matrix[N, N] Sigma_total;
  matrix[m, m] Sigma22;
  matrix[n, m] Sigma12;
  matrix[n, n] Sigma11_inv;
  matrix[n, m] prod;
  matrix[m, n] prod_t;
  matrix[m, m] Sigma2_1;
  matrix[m, m] Sigma2_1_fixed;
  vector[m] mu2_1;
  matrix[m, m] L_Sigma2_1;
  real linear_pred[m];
  real mu[m];

  // assignment
  Sigma_total = cov_exp_quad(coords_all, sigma, rho);
  Sigma22 = Sigma_total[(n+1):N,(n+1):N];
  Sigma12 = Sigma_total[1:n, (n+1):N];
  Sigma11_inv = inverse(Sigma11);

  // joint MVN prediction
  prod = Sigma11_inv * Sigma12;
  prod_t = prod';
  Sigma2_1 = Sigma22 - prod_t * Sigma12;
  mu2_1 = prod_t * f;

  // resolve symmetry issues due to rounding
  for(row in 1:(m-1))
    for(col in (row+1):m)
      Sigma2_1[col, row] = Sigma2_1[row, col];

  // resolve p.d. issues due to rounding
  for(k in 1:m)
    Sigma2_1[k, k] = Sigma2_1[k, k] + delta;
  
  // generate MVN RV's
  L_Sigma2_1 = cholesky_decompose(Sigma2_1);
  f_pred = L_Sigma2_1 * eta2 + mu2_1;

  // poisson process
  for(i in 1:m){
    linear_pred[i] = f_pred[i] * beta;
    mu[i] = exp(linear_pred[i]);
  }
}
model {
  rho ~ inv_gamma(6.7, 28.4);
  sigma ~ std_normal();
  tau ~ std_normal();
  eta ~ std_normal();
  eta2 ~ std_normal();
  beta ~ normal(0, 100);

  y ~ normal(f, tau);


}

