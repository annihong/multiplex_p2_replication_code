library(dplyr)
library(doSNOW)
library(rstan)
#library(foreach)
source("./helper.R")
source("./simulation_helper.R")
# source("../stan_helper.R")
if (!require("foreach")) install.packages("foreach")
#PARALLEL INFO
NUM_CORE = parallel::detectCores() - 1
TOTAL_ITER = 1000
EXPERIMENT = "simulation"
OUTPUT_PATH = "" ## path to save the results

# RSTAN CONFIG:
CHAINS = 4
WARMUP = 2000
THIN = 1

# NETWORK INFO: 
n = 30
t = 2
H = t*(t - 1)/2
N = n*(n - 1)/2 # number of dyads
K = 2^(2*t)

# true mean of the priors
#Sigma = diag(2*t)
mu_0 = rep(0,t)
rho_0 = rep(0,t)
cross_mu_0 = rep(0,H)
cross_rho_0 = rep(0,H)
eta = 2
Sigma_Dist <- "LJK" 
nu = 2*t + 1
ig_shape = 3
ig_scale = 50

######## BEGIN SIMULATING NETWORK FROM GROUND TRUTH #########

prior_res <- list()

for (i in 1:TOTAL_ITER) {
  Sigma <- sample_Sigma()
  sampled_params <- sample_prior(n, t, mu_0, rho_0, cross_mu_0, cross_rho_0, Sigma)
  simulated_network <- simulate_network(n, t, sampled_params)
  prior_res[[i]] <- list(sampled_params = sampled_params, simulated_network = simulated_network)
}



########  END SIMULATING NETWORK FROM GROUND TRUTH #########
cl <- parallel::makeCluster(NUM_CORE, outfile=paste0(OUTPUT_PATH, "log/", EXPERIMENT, "LogSim.txt"))  # Create a cluster with 4 workers
doSNOW::registerDoSNOW(cl)  # Register the cluster for use with foreach

#results <- foreach(iter=1:TOTAL_ITER, .packages="rstan") %dopar% {
results <- foreach(i=1:length(iter_ids), .packages = "multip2") %dopar% {
    iter = iter_ids[i]
    sampled_params <- prior_res[[iter]]$sampled_params
    sampled_params$iter = iter
    M <- prior_res[[iter]]$simulated_network #network to fit:
    m_fit <- multip2::Mp2Model(M)
    m_fit <-  multip2::fit(m_fit, chains = CHAINS, warmup = WARMUP, thin = THIN, iter = WARMUP*2, network_sim = TRUE, seed = iter, par = "x_beta", include = FALSE, stan_file = "multiplex_p2_revert.stan")
    # m_fit <-  multip2::fit(m_fit, chains = 1, warmup =10, thin = 1, iter = 1000, network_sim = FALSE, prior_sim = TRUE, stan_file = "multiplex_p2_low_mem_prior_sim_test.stan")


    sim_result <- list(sampled_params=sampled_params,Mp2_fit = m_fit)
    saveRDS(sim_result,file = paste0(OUTPUT_PATH, iter, "_", EXPERIMENT, "draws = ", CHAINS * WARMUP, ".Rds"))
    rm(sim_result)
    rm(m_fit)
}

stopCluster(cl)
