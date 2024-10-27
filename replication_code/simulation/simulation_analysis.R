set.seed(1)
library(multip2)
path <- ""
MODE = "PS_new"
simulations <- list.files(path)
length(simulations)
library(ggplot2)
library(dplyr)

sampled_params_list <- list()
posterior_draws_list <- list()


zscore_helper <- function(true_param, vector) {
    post_mean = mean(vector)
    post_sd = sd(vector)
    zscore = (post_mean-true_param)/post_sd
    #print(paste0("post_mean: ", post_mean, "post_sd", post_sd))
    return(zscore)
}

post_contract_helper <- function(to_rank, vector, variance = 100) {
    post_sd = sd(vector)
    post_contract = 1 - (post_sd^2)/variance
    return(post_contract)
}


rank_helper <- function(to_rank, vector){
    subsample <- sample(vector,999, replace = FALSE)
    #subsample <- vector[-length(vector)]
    # print(length(subsample))
    rank <- rank(c(to_rank, subsample))[1]
    return(rank)
}
# given the parameter name, 
calc_stats <- function(sampled_params, posterior_draws, param, stat_func){
    prior_params <- sampled_params[[param]]
    draws <- posterior_draws[[param]]
    # print(param)
    # print(nrow(draws))
    post_draws <- as.matrix(draws)
    stats <- c()

    for (t in 1:length(prior_params)) {
        prior_param = prior_params[[t]]
        if (MODE == "PS" & identical(stat_func, post_contract_helper) & param == "mu") {
            n = nrow(sampled_params$C)
            A_idx = 1 + 2 * (t - 1)
            B_idx = 2 + 2 * (t - 1)
            v_mu = 100 
            v_A = sampled_params$Sigma[A_idx, A_idx]
            v_B = sampled_params$Sigma[B_idx, B_idx]
            cov_AB = sampled_params$Sigma[A_idx, B_idx]
            v_PS_mu = v_mu + (1/n) * (v_A + v_B + 2 * cov_AB)
            stat <- stat_func(prior_param, post_draws[,t], variance = v_PS_mu)
        } else {
            stat <- stat_func(prior_param, post_draws[,t])
        }
        stats <- c(stats, stat)
    }
    names(stats) <- paste(param, 1:length(prior_params), sep = "_")
    return(stats)
}

calc_stats_file <- function(simulation_file){
    print(simulation_file)
    file_path <- paste0(path, folder, simulation_file)
    #print(simulation_file)
    # e1 <- new.env(parent = baseenv())
    # load(file_path, envir = e1)

    sim_result <- readRDS(file_path)
    stan_fit <- sim_result$Mp2_fit$fit_res$stan_fit
    rstan::check_hmc_diagnostics(stan_fit)
    sampled_params <- sim_result$sampled_params
    if (MODE == "OG") {
        posterior_draws <- rstan::extract(stan_fit, c("mu", "rho", "cross_mu", "cross_rho"))
    } else if (MODE == "PS"){
        posterior_draws <- rstan::extract(stan_fit, c("mu", "rho", "cross_mu", "cross_rho", "A_bar", "B_bar"))
        PS_mu = posterior_draws$mu + posterior_draws$A_bar + posterior_draws$B_bar
        posterior_draws$mu <- PS_mu
        posterior_draws <- posterior_draws[c("mu", "rho", "cross_mu", "cross_rho")]
        t_total = ncol(PS_mu)
        A_bar_prior <- colMeans(sampled_params$C[,1 + 2*(1:t_total - 1)])
        B_bar_prior <- colMeans(sampled_params$C[,2 + 2*(1:t_total - 1)])
        sampled_params$mu <- sampled_params$mu + A_bar_prior + B_bar_prior
    } else if (MODE == "PS_new") {
        posterior_draws <- rstan::extract(stan_fit, c("PS_mu", "rho", "cross_mu", "cross_rho"))
        posterior_draws$mu <- posterior_draws$PS_mu
        posterior_draws <- posterior_draws[c("mu", "rho", "cross_mu", "cross_rho")]
        t_total = ncol(posterior_draws$mu)
        A_bar_prior <- colMeans(sampled_params$C[,1 + 2*(1:t_total - 1)])
        B_bar_prior <- colMeans(sampled_params$C[,2 + 2*(1:t_total - 1)])
        sampled_params$mu <- sampled_params$mu + A_bar_prior + B_bar_prior
    }
    

    #posterior_draws <- sim_result$posterior_draws

    params <- names(sampled_params)[1:4] # no sigma
    names(posterior_draws) <- params # rename PS_mu as mu
    res <- list(ranks = c(), zscores = c(), posts = c())
    stat_funcs <- c(rank_helper, zscore_helper, post_contract_helper)
    for (param in params) {
        for (i in 1:length(res)) {
            res[[i]] = c(res[[i]], calc_stats(sampled_params, posterior_draws, param, stat_funcs[[i]]))
        }
    }
    res$sampled_params <- sampled_params
    res$posterior_draws <- posterior_draws
    return(res)
}


res <- lapply(simulations, calc_stats_file)

rank_df <- data.frame()
zscore_df <- data.frame()
post_contract_df <- data.frame()
for (file_res in res) {
    rank_df <- rbind(rank_df,file_res$ranks)
    zscore_df <- rbind(zscore_df,file_res$zscores)
    post_contract_df <- rbind(post_contract_df, file_res$posts)
}

colnames(rank_df) <- colnames(zscore_df) <- colnames(post_contract_df) <- names(res[[1]]$ranks)

## This script is used to generate the plots for the simulation results modified from the SBC r package: https://github.com/hyunjimoon/SBC/blob/master/R/plot.R
library(tidyverse)
library(SBC)
source("./plot_helper.R")
rank_df <- res$rank_df
zscore_df=res$zscore_df
post_contract_df=res$post_contract_df


plot_long <- rank_df %>% gather(key = "variable", value = "rank")
plot_long$variable <- factor(plot_long$variable, levels = colnames(rank_df))
max_rank =1000


g <- plot_rank_hist(plot_long, max_rank = max_rank, bins=30)
ggsave("rank_1000.png", plot=g, width = 10, height = 5, dpi = 300)

#ecdf_data <- data_for_ecdf_plots.matrix(rank_df)
g <- plot_ecdf(rank_df, max_rank = max_rank)
ggsave("ecdf_1000.png", plot=g, width = 10, height = 5, dpi = 300)


#### Model Sensitivity ####

zscore_long <- tidyr::gather(zscore_df, value = "zscore", key = "params")
scatter_df <- tidyr::gather(post_contract_df, value = "post_contract", key = "params")
scatter_df[,"z_score"] <- zscore_long$zscore

scatter <- ggplot(scatter_df,aes(x=post_contract,y=z_score))+ geom_point(aes(color=params), alpha = 0.4, size = 1) + facet_wrap( ~ params, ncol=6, scales = "free_x") + xlab("Posterior contraction") + ylab("Z-score") + 
ylim(-5,5) +
theme_bw() + 
scale_color_manual(values=c(rep("#fb8500", 1),  rep("#ffb703", 1), rep("#219ebc", 2), rep("#5AB1BB", 2))) + theme(legend.position="none", axis.text.x = element_text(size = 5), axis.text.y = element_text(size = 5),
      axis.title = element_text(size = 15))

ggsave("model_sensitivity_1000.png", plot=scatter, width = 8, height = 3)


