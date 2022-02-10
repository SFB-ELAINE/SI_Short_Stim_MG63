# Script for combining ROS and beta-catenin measurements and +++++++++++++++
# calculating correlation factor
# Author: Kai Budde
# Created: 2022/02/09
# Last changed: 2022/02/09

rm(list=ls())
# Sys.setlocale("LC_TIME","English")

# General function parameters ##############################################
options(stringsAsFactors = FALSE, warn=-1)

library("tidyverse")
library("ggplot2")
library("ggpubr")

# Input parameters #########################################################

ROS_beta_cat_data <- "data/ROS_and_Beta_Catenin_data.csv"
output_dir <- "plots"

# Data input, cleaning, and additional columns #############################
df_complete <- read.csv(file = ROS_beta_cat_data)


# Plot ratios of mean fluorescence intensities #############################
df_final_beta_catenin <-  df_complete %>%
  group_by(TimePoint, Experiment) %>%
  summarize(ratio_mean_fluorescence = FL2.H.mean[Group == "stim"]/
              FL2.H.mean[Group == "control"])

df_final_beta_catenin <- as.data.frame(df_final_beta_catenin)

for(i in 0:(max(df_final_beta_catenin$TimePoint)+1)){
  if(!(i %in% df_final_beta_catenin$TimePoint)){
    df_final_beta_catenin <- rbind(df_final_beta_catenin,
                                   c(i, rep(NA, length(names(df_final_beta_catenin))-1) ))
  }
}

df_final_beta_catenin$TimePoint <- as.factor(df_final_beta_catenin$TimePoint)

plot_ratios <- ggplot(df_final_beta_catenin, aes(x=TimePoint, y=ratio_mean_fluorescence)) +
  stat_boxplot(geom ='errorbar', width = 0.3, position = position_dodge(width = 0.75)) +
  geom_boxplot(alpha = 1, position = position_dodge2(width = 0.9, preserve = "single"), outlier.shape = 1) +
  geom_jitter(color="black", size=0.5, alpha=0.9) +
  # geom_point(position = position_jitterdodge(), size=1, alpha=0.9) +
  #ylim(0,20) +
  theme_bw(base_size = 18) +
  theme(#axis.title.y=element_text(size=12),
    #axis.text.x = element_blank(), 
    axis.ticks.x = element_blank()) +
  ylab("Beta-catenin ratio of mean fluorescence (stim vs. control)") +
  xlab("Time in h")

# file_name <- paste(output_dir, "/Beta_catenin_ratio.png", sep="")
# ggsave(filename = file_name, width = (297), height = (210), units = "mm")

# Plot correlation #########################################################

# plot_correlation <- ggplot(df_complete, aes(x=mean_fluorescence, y=FL2.H.mean)) +
#   geom_point() +
#   #ylim(0,20) +
#   theme_bw(base_size = 18) +
#   theme(#axis.title.y=element_text(size=12),
#     #axis.text.x = element_blank(), 
#     axis.ticks.x = element_blank()) +
#   ylab("Mean ROS fluorescence in arb. units") +
#   xlab("Mean Beta-catenin fluorescence in arb. units") +
#   stat_cor(method = "pearson")


plot_correlation <- ggscatter(df_complete, x = "FL2.H.mean", y = "mean_fluorescence",
                              add = "reg.line",  # Add regressin line
                              add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
                              conf.int = TRUE # Add confidence interval
                              ) +
  stat_cor(method = "pearson") +
  xlab("Mean ROS fluorescence in arb. units") +
  ylab("Mean Beta-catenin fluorescence in arb. units")

file_name <- paste(output_dir, "/correlation_bcat_ROS1.png", sep="")
ggsave(filename = file_name, width = (297), height = (210), units = "mm")

plot_correlation2 <- ggscatter(df_complete, x = "FL2.H.mean", y = "mean_fluorescence", color = "Group",
                              add = "reg.line",  # Add regressin line
                              palette = "jco",
                              conf.int = TRUE # Add confidence interval
                              ) +
  stat_cor(aes(color = Group)) +
  xlab("Mean ROS fluorescence in arb. units") +
  ylab("Mean Beta-catenin fluorescence in arb. units")

file_name <- paste(output_dir, "/correlation_bcat_ROS2.png", sep="")
ggsave(filename = file_name, width = (297), height = (210), units = "mm")


# plot_correlation <- ggplot(df_complete, aes(x=mean_fluorescence, y=FL2.H.mean, color = Group)) +
#   geom_point() +
#   #ylim(0,20) +
#   theme_bw(base_size = 18) +
#   theme(#axis.title.y=element_text(size=12),
#     #axis.text.x = element_blank(), 
#     axis.ticks.x = element_blank()) +
#   ylab("Mean ROS fluorescence in arb. units") +
#   xlab("Mean Beta-catenin fluorescence in arb. units")

