file_path <- "/Users/annihong/Documents/Rprojects/multiplex-social-networks/revision/study1/"


##LIBRARIES##
require(tidyverse)
library(multip2)
# require(ggmcmc)
# require(igraph)
# require(network)
# require(ergm)
# library(viridis)
# require(statnet.common)


##READ IN DATA##
m_1 <- readRDS(paste0(file_path, "study_1_models_m1.Rds"))$m_1
m_2 <- readRDS(paste0(file_path, "study_1_models_m2.Rds"))$m_2



##LOAD IN THE STAN FIT OBJECTS##
fit1 <- m_1$fit_res$stan_fit
fit2 <- m_2$fit_res$stan_fit

##MODEL ESTIMATES##

s1 <- multip2::summary.Mp2Model(m_1)

s2 <- multip2::summary.Mp2Model(m_2)

##CONVERGENCE CHECK## NEEDS REVISION (figure out PS or HC, change the stan model used to estimation)

#S <- ggs(fit2, family = "PS|rho|Corr|fixed")
S <- ggs(fit2, family = "PS")
S <- ggs(fit2, family = "^Corr\[(\d+),(?=\d+\]$)(\d+)\]$")
ggmcmc(S, file = paste0("diagnostics",".pdf"),  family = "PS|rho", param_page=6) #plot=c("traceplot", "running", "geweke"),
ptrace <- ggs_traceplot(S, family = "^Corr") +
      facet_wrap(~ Parameter, ncol = 6, scales = "free")
p <- ggs_grb(S, family = "^mu") + facet_wrap(~ Parameter, ncol = 3, scales = "free")
p <- ggs_grb(S, family = "PS_mu") + facet_wrap(~ Parameter, ncol = 3, scales = "free")
ggs_geweke(S, family = "mu|rho") + facet_wrap(~ Parameter, ncol = 3, scales = "free")
ggs_running(S, family = "Sigma") + facet_wrap(~ Parameter, ncol = 6, scales = "free")
ggs_density(ggs(radon$s.radon, par_labels=P, family="sigma"))
ptrace <- ggs_traceplot(S, family = "mu")
ggs_running(S, family = "^PS_mu")
ggs_Rhat(S, family = "^mu")
ggs_geweke(S, family = "^mu")
ggsave(ptrace)




##### RANDOM EFFECTS M2
t = m_1$t
multi_result_corr <- s2$correlation

# labs <- rownames(multi_result_corr)
# sender_labs <- grep("^sender", labs, value = T)
# receiver_labs <- grep("^receiver", labs, value = T)
labs <- c(paste("sender", c("PO", "SC", "PE"), sep =":"), paste("receiver", c("PO", "SC", "PE"), sep =":"))
corr_M <- matrix(, ncol = 2*t, nrow = 2*t, dimnames = list(labs, labs))

fill_corr_M <- function(corr_M, res) {
      res_M <- corr_M 
      for (col in colnames(corr_M)){
            for (row in rownames(corr_M)){
                  lab <- paste0(col, "_", row)
                  if (lab %in% names(res)) {
                        res_M[col, row] = res[lab]
                  } else {
                        lab <- paste0(row, "_", col)
                        res_M[col, row] = res[lab]
                  }
            }
      }
      return(res_M)
}


#sig_dat <- fill_corr_M(corr_M, multi_result_corr[, '2.5%'] * multi_result_corr[, '97.5%'] > 0 )
# val_dat <- fill_corr_M(corr_M, multi_result_corr[, "mean"])
# upper_dat <- fill_corr_M(corr_M, multi_result_corr[, '97.5%'])
# lower_dat <-fill_corr_M(corr_M, multi_result_corr[, '2.5%'])
df_corr <- expand.grid(Var1 = rownames(corr_M), Var2 = colnames(corr_M))
df_corr$value <- as.vector(fill_corr_M(corr_M, multi_result_corr[, "mean"]))
df_corr$upper <- as.vector(fill_corr_M(corr_M, multi_result_corr[, '97.5%']))
df_corr$lower <- as.vector(fill_corr_M(corr_M, multi_result_corr[, '2.5%']))
df_corr$CI_low <- paste0("(", round(df_corr$lower , 2))
df_corr$CI_high <- paste0(round(df_corr$upper, 2), ")")
df_corr$CI <- paste0(df_corr$CI_low, ", ", df_corr$CI_high)

val_dat_tri <- df_corr 
# %>%
#   rowwise() %>%
#   mutate(pair = sort(c(Var1, Var2)) %>% paste(collapse = ",")) %>%
#   group_by(pair) %>%
#   distinct(pair, .keep_all = T)


val_dat_tri$Var1 <- factor(val_dat_tri$Var1, levels = unique(val_dat_tri$Var1))
val_dat_tri$Var2 <- factor(val_dat_tri$Var2, levels = unique(val_dat_tri$Var2))



p_sigma <- ggplot(val_dat_tri, aes(Var1, Var2)) +
      geom_tile(aes(fill = value), show.legend = FALSE) +
      geom_text(aes(label = CI), size = 5, nudge_y=-0.2) +
      scale_fill_gradient2(low = "#219ebc", high = "#ffb703", name = "Estimated Covariance") + 
      #scale_color_gradient2() + 
      geom_text(aes(label = round(value, 2)), size=6, nudge_y=0.1) +
      xlab("") + ylab("") +
      theme_minimal() + 
      theme(axis.text = element_text(size = 16, color = "black"))+
      scale_y_discrete(position = "left")

ggsave("/Users/annihong/Documents/Rprojects/multiplex-social-networks/revision/study1/model2_corr.png", plot=p_sigma, width = 10, height = 8.5)


##POSTERIOR PREDICTIVE CHECKS##
sim_nets <- extract_network_draws(results$m_2, 1000)
dep_net <- results$m_2$data$dep_net


png("/Users/annihong/Documents/Rprojects/multiplex-social-networks/revision/study1/triad_census.png", height = 400, width = 1200, res = 100)
simulated_network_checks(sim_nets, dep_net, "Triad_census") 
dev.off()

png("/Users/annihong/Documents/Rprojects/multiplex-social-networks/revision/study1/id.png", height = 400, width = 1200, res = 100)
simulated_network_checks(sim_nets, dep_net, "Indegree_distribution")
dev.off()

png("/Users/annihong/Documents/Rprojects/multiplex-social-networks/revision/study1/od.png", height = 400, width = 1200, res = 100)
simulated_network_checks(sim_nets, dep_net, "Outdegree_distribution")
dev.off()

png("/Users/annihong/Documents/Rprojects/multiplex-social-networks/revision/study1/baselineGOF.png", height = 600, width = 1000, res = 100)
simulated_network_checks(sim_nets, dep_net, "multiplex_gof_baseline") + labs(title = "")
dev.off()


png("/Users/annihong/Documents/Rprojects/multiplex-social-networks/revision/study1/randomGOF.png", height = 600, width = 800, res = 100)
simulated_network_checks(sim_nets, dep_net, "multiplex_gof_random")  + labs(title = "")
dev.off()

######old way#####
t = 3
n = 30
M_net <- lapply(dep_net, igraph::graph_from_adjacency_matrix)
observed_stats <- as.data.frame(descriptive_stats(M_net))
observed_stats_df <- data.frame("var" = rownames(observed_stats), "sim_stats" = observed_stats[,1])
sim_stats <- observed_stats
# post processing the generated results
sim_sample_size = 4000
fit = fit1 # use model 2 for gof
#sim_nws <- tail(rstan::extract(fit)$y_tilde, sim_sample_size)
sim_nws <-rstan::extract(fit)$y_tilde[sample(1:8000, sim_sample_size),]
sim_nws_igraphs <- dyads_to_matrix_list(sim_nws, n, t, "igraph") #convert the simulated dyad data to igraph objects 

descriptive_labels <- c("Pol Info Density", 
                "Sci Info Density", 
                "Influ Density",
                "Pol Info Reciprocity", 
                "Sci Info Reciprocity", 
                "Influ Reciprocity",
                "Pol X Sci Jaccard",
                "Pol X Influ Jaccard",
                "Influ X Sci Jaccard",
                "Pol X Sci Cross Jaccard",
                "Pol X Influ Cross Jaccard",
                "Influ X Sci Cross Jaccard")

#calculate the descriptive stats for all the simulated data 
stats <- descriptive_stats_list(sim_nws_igraphs, avg=F)#[,c(1,2,3,7,8,9,4,5,6,10,11,12)]
basic_stats <- gather(data.frame(stats), key="var", value="sim_stats")
#varcov_stats <- gather(data.frame(stats[,7:22]), key="var", value="sim_stats")
#basic_stats$var <- factor(basic_stats$var, colnames(stats)[c(1,2,3,7,8,9,4,5,6,10,11,12)])
#levels(basic_stats$var) <- rev(levels(basic_stats$var))
p3 <- ggplot(basic_stats,aes(x = var, y=sim_stats, fill=var),show.legend = FALSE) +
 stat_boxplot(geom ='errorbar') +
    geom_boxplot() +
    scale_fill_manual(values=c(rep("#219ebc", 3), rep("#5AB1BB", 3), rep("#fb8500", 3), rep("#ffb703", 3)))+
    geom_point(data=observed_stats_df, mapping = aes(var,sim_stats), shape = 21, colour = "black", fill = "white", size=3) +
    theme_bw() + 
    theme(
      legend.position="none",
      plot.title = element_text(size=11),
      axis.text.y = element_text(size=15, color = "black"),
      axis.text.x = element_text(size=15, color = "black")) +
    #scale_x_discrete(labels= descriptive_labels) + 
    xlab("") +ylab("") + coord_flip()
    #ggtitle("Descriptive Statistics Observed vs Simulated basic statistics") +


ggsave(paste0(file_path, "multiplex_PPC_basic.png"), plot=p3,  width = 10, height = 5)





