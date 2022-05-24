# Script for plotting histograms of calculated H2O2 concentrations +++++++++
# Author: Kai Budde
# Created: 2022/03/21
# Last changed: 2022/03/21

rm(list=ls())

# Input parameters #########################################################

ROS_data <- "data/H2O2/H2O2_in_stimulated_medium.csv"
output_dir <- "plots/H2O2"
# filter_out_timepoint <- "4h"
# filter_out_experiment_nr <- 2


# General function parameters ##############################################
options(stringsAsFactors = FALSE, warn=-1)

library(tidyverse)
library(ggplot2)
# library(cowplot)

input_dir <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(dir = input_dir)
setwd(dir = "..")

# Data input, cleaning, and additional columns #############################

df_H2O2 <- read_csv(file = ROS_data, name_repair = "universal", show_col_types = FALSE)

# Date column as date
df_H2O2$Date <- lubridate::dmy(df_H2O2$Date)

# Remove specified rows
# df_H2O2 <- df_H2O2 %>%
#   dplyr::filter(!(TimePoint == filter_out_timepoint & Experiment == filter_out_experiment_nr))

# Save tibble with ratios  #################################################

df_final_1 <- df_H2O2 %>%
  group_by(TimePoint_in_hours, Group) %>%
  summarize(mean_H2O2_Concentration_in_uM = mean(H2O2_Concentration_in_uM, na.rm = TRUE),
            sd_H2O2_Concentration_in_uM = sd(H2O2_Concentration_in_uM, na.rm = TRUE))

df_final_2 <- df_final_1 %>%
  group_by(TimePoint_in_hours) %>%
  summarize(ratio_mean_H2O2_Concentration = mean_H2O2_Concentration_in_uM[Group == "stim"]/
              mean_H2O2_Concentration_in_uM[Group == "control"],
            sum_sd_H2O2_Concentration = sd_H2O2_Concentration_in_uM[Group == "stim"] +
              sd_H2O2_Concentration_in_uM[Group == "control"])

# Plotting #################################################################

dir.create(path = output_dir, showWarnings = FALSE)

# Plot ratio of H2O2 over time

plot_H2O2_ratio <- ggplot(df_final_2, aes(x=TimePoint_in_hours,ratio_mean_H2O2_Concentration )) +
  geom_point(size = 2) +
  geom_line() +
  geom_vline(xintercept = 0, linetype="dotted", color = "blue", size=1) +
  geom_text(aes(x=0, label="stim on", y=15.5), hjust=-0.1, colour="blue", text=element_text(size=11)) +
  geom_vline(xintercept = 2, linetype="dotted", color = "blue", size=1) +
  geom_text(aes(x=2, label="stim off", y=7.5), hjust=-0.1, colour="blue", text=element_text(size=11)) +
  geom_errorbar(aes(ymin=ratio_mean_H2O2_Concentration-sum_sd_H2O2_Concentration,
                    ymax=ratio_mean_H2O2_Concentration+sum_sd_H2O2_Concentration), width=.2) +
  labs(title="H2O2 conentration ratio (stim vs. control)",
       x="Time in h", y = "Ratio") +
  scale_x_continuous(expand = c(0, 0), limits = c(-1, 37), minor_breaks = c(-1:37),
                     breaks = c(0,2,4,8,12,24,36)) + 
  theme_bw()

# print(plot_ROS_ratio)

file_name <- paste(output_dir, "/H2O2_concentration_ratio.png", sep="")
ggsave(filename = file_name, width = 297, height = 210, units = "mm")


# Boxplots -----------------------------------------------------------------

# Add empty timepoints
for(i in 0:(max(df_H2O2$TimePoint_in_hours)+1)){
  if(!(i %in% df_H2O2$TimePoint_in_hours)){
    df_H2O2 <- rbind( df_H2O2, c(NA, i, rep(NA, length( names(df_H2O2) )-2) ) )
  }
}

# Save time points as factors
df_H2O2$TimePoint_in_hours <- as.factor(df_H2O2$TimePoint_in_hours)

plot_H2O2_boxplots <- ggplot(df_H2O2, aes(x = TimePoint_in_hours, y = H2O2_Concentration_in_uM, fill = Group)) +
  stat_boxplot(geom ='errorbar', width = 0.3, position = position_dodge(width = 0.75)) +
  geom_boxplot(alpha = 1, position = position_dodge2(width = 0.9, preserve = "single"), outlier.shape = 1) +
  geom_jitter(color="black", size=0.5, alpha=0.9, position = position_dodge2(width = 0.9, preserve = "single")) +
  theme_bw(base_size = 18) +
  theme(#axis.title.y=element_text(size=12),
    #axis.text.x = element_blank(), 
    axis.ticks.x = element_blank()) +
  ylab("H2O2 concentration in uM") +
  xlab("Time in h")

file_name <- paste(output_dir, "/mean_sd_plot_FL2H_log.png", sep="")
ggsave(filename = file_name, width = (297), height = (210), units = "mm")
