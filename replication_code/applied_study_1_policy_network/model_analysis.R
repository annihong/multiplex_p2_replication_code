file_path <- ""


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


png("triad_census.png", height = 400, width = 1200, res = 100)
simulated_network_checks(sim_nets, dep_net, "Triad_census") 
dev.off()

png("id.png", height = 400, width = 1200, res = 100)
simulated_network_checks(sim_nets, dep_net, "Indegree_distribution")
dev.off()

png("od.png", height = 400, width = 1200, res = 100)
simulated_network_checks(sim_nets, dep_net, "Outdegree_distribution")
dev.off()

png("baselineGOF.png", height = 600, width = 1000, res = 100)
simulated_network_checks(sim_nets, dep_net, "multiplex_gof_baseline") + labs(title = "")
dev.off()


png("randomGOF.png", height = 600, width = 800, res = 100)
simulated_network_checks(sim_nets, dep_net, "multiplex_gof_random")  + labs(title = "")
dev.off()




