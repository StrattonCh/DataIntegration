library(tidyverse)

mh_binom <- function(seed, num_mcmc, warmup = num_mcmc/2, y, a = 1, b = 1, prop_sd = .01){
  # seed
  set.seed(seed)
  
  # convenience
  n <- length(y)
  
  # storage
  pi_mcmc <- matrix(NA, num_mcmc, 1)
  accept_ratio <- matrix(0, num_mcmc, 1)
  
  # initialize
  pi <- rbeta(1, 1, 1)
  
  message(paste0("Beginning sampling at "), Sys.time())
  pb <- txtProgressBar(min = 0, max = num_mcmc, style = 3, width = 50, char = "=")
  for(iter in 2:num_mcmc){
    # pi - metropolis
    # proposal
    pi_s <- -1
    while(pi_s <= 0 | pi_s >= 1) pi_s <- pi + rnorm(1, 0, prop_sd)

    ## evaluate proposal
    log_post_current <- sum(dbinom(y, 1, pi, log = T)) + dbeta(pi, 1, 1, log = T)
    log_post_s <- sum(dbinom(y, 1, pi_s, log = T)) + dbeta(pi_s, 1, 1, log = T)
    log_r <- log_post_s - log_post_current
    if(log(runif(1)) < log_r){
      pi <- pi_s
      accept_ratio[iter,] <- 1
    }
    
    # storage
    pi_mcmc[iter,] <- pi

    # progress
    setTxtProgressBar(pb, iter)
  }
  close(pb)
  message(paste0("Ending sampling at "), Sys.time())
  
  samples <- cbind(
    pi_mcmc[(warmup+1):num_mcmc,]
  )
  
  colnames(samples) <- c("pi")
  
  return(
    list(
      samples = samples,
      accept_ratio = accept_ratio
    )
  )
  
  
  
}

df <- tibble(
  y = rbinom(500, 1, .5)
)

samples <- mh_binom(
  seed = 1,
  num_mcmc = 5000,
  y = df$y,
  prop_sd = .05
)

mean(samples$accept_ratio)
mean(samples$samples)
quantile(samples$samples, c(.025, .975))
plot(samples$samples[,1], type = "l")

# seed = 1
# num_mcmc = 5000
# warmup = 2500
# y = df$y
# a = 1
# b = 1
# initial_prop_var = .1
# adapt_period = .1*num_mcmc
# sd = 2.4^2
# epsilon = 1e-9

mh_binom_adapt <- function(seed, num_mcmc, warmup = num_mcmc/2, y, a = 1, b = 1,
                           initial_prop_var = .01, adapt_period = .1*num_mcmc, sd = 2.4^2, epsilon = 1e-9){
  # seed
  set.seed(seed)
  
  # convenience
  n <- length(y)
  
  # storage
  pi_mcmc <- matrix(NA, num_mcmc, 1)
  Ct_mcmc <- matrix(NA, num_mcmc, 1)
  accept_ratio <- matrix(0, num_mcmc, 1)
  
  # initialize
  pi <- rbeta(1, 1, 1);pi_mcmc[1,] <- pi
  
  # adaptive
  Ct <- initial_prop_var; Ct_mcmc[1,] <- Ct
  
  message(paste0("Beginning sampling at "), Sys.time())
  pb <- txtProgressBar(min = 0, max = num_mcmc, style = 3, width = 50, char = "=")
  for(iter in 2:num_mcmc){
    # pi - metropolis
    # proposal
    pi_s <- -1
    if(iter <= adapt_period){
      Ct <- initial_prop_var
    } else if(iter > adapt_period){
      Ct <- sd * cov(pi_mcmc[1:(iter-1),,drop = F]) + sd * epsilon
      Ct_mcmc[iter,] <- Ct
    }
    while(pi_s <= 0 | pi_s >= 1) pi_s <- pi + rnorm(1, 0, sqrt(Ct))
    
    ## evaluate proposal
    log_post_current <- sum(dbinom(y, 1, pi, log = T)) + dbeta(pi, 1, 1, log = T)
    log_post_s <- sum(dbinom(y, 1, pi_s, log = T)) + dbeta(pi_s, 1, 1, log = T)
    log_r <- log_post_s - log_post_current
    if(log(runif(1)) < log_r){
      pi <- pi_s
      accept_ratio[iter,] <- 1
    }
    
    # storage
    pi_mcmc[iter,] <- pi
    
    # progress
    setTxtProgressBar(pb, iter)
  }
  close(pb)
  message(paste0("Ending sampling at "), Sys.time())
  
  samples <- cbind(
    pi_mcmc[(warmup+1):num_mcmc,],
    Ct_mcmc[(warmup+1):num_mcmc,]
  )
  
  colnames(samples) <- c("pi", "Ct")
  
  return(
    list(
      samples = samples,
      accept_ratio = accept_ratio
    )
  )
}

samples <- mh_binom_adapt(
  seed = 1,
  num_mcmc = 5000,
  y = df$y,
  initial_prop_var = .05^2
)

mean(samples$accept_ratio)
mean(samples$samples[,1])
quantile(samples$samples[1], c(.025, .975))
plot(samples$samples[,2], type = "l")
