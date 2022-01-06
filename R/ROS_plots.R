# Script for plotting histograms of ROS readout using TECAN reader +++++++++
# Author: Kai Budde
# Created: 2022/01/04
# Last changed: 2022/01/04

rm(list=ls())

# Input parameters #########################################################

ROS_data <- "data/ROS_raw_data.csv"
output_dir <- "plots"

# General function parameters ##############################################
options(stringsAsFactors = FALSE, warn=-1)

library(tidyverse)
library(ggplot2)

input_dir <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(dir = input_dir)
setwd(dir = "..")

# Data input, cleaning, and additional columns #############################

df_ROS <- read_csv(file = ROS_data, name_repair = "universal", show_col_types = FALSE)

# Rename columns of the image acquisition positions (first number is col,
# second number is row)

names(df_ROS) <- gsub(pattern = "\\.\\.([0-9]{1,2})\\.([0-9]{1,2})",
                      replacement = "pos_\\1_\\2", x = names(df_ROS))
positions <- grep(pattern = "pos", x = names(df_ROS), value = TRUE)

# Change NAs in fluorescence data to max possible value
max_possible_fluorescence_value <- 2^16-1
df_ROS[, positions][is.na(df_ROS[, positions])] <- max_possible_fluorescence_value

# Add column with mean and std (the "Mean" and "StDev" come from the
# TECAN software)
df_ROS <- df_ROS %>%
  mutate(calc_mean = rowMeans(.[grep(pattern = "pos", x = names(.))], na.rm = TRUE))
df_ROS <- df_ROS %>% relocate(calc_mean, .after = Mean)

df_ROS$calc_sd <- apply(df_ROS[,grep(pattern = "pos", x = names(df_ROS))], 1, sd, na.rm = TRUE)
df_ROS <- df_ROS %>% relocate(calc_sd, .after = StDev)

# Include copy of control group at 0h as stim data
df_ROS_copy <- df_ROS[df_ROS$TimePoint == "0h" & df_ROS$Group == "control",]
df_ROS_copy$Group <- "stim"
df_ROS <- rbind(df_ROS, df_ROS_copy)
rm(df_ROS_copy)

# Lines for testing purpose
# df_ROS$test <- df_ROS$calc_sd - df_ROS$StDev
# df_ROS <- df_ROS %>% relocate(test, .after = StDev)
# max(abs(df_ROS$test))
# any(is.na(df_ROS$test))

# Save tidy tibble #########################################################
df_tidy_ROS <- df_ROS %>%
  select(c("TimePoint", "Experiment", "Group", "Well", positions))

df_tidy_ROS <- df_tidy_ROS %>%
  gather(key = position, value = fluorescence, -c(1:4))
df_tidy_ROS <- df_tidy_ROS[!is.na(df_tidy_ROS$fluorescence),]
# any(is.na(df_tidy_ROS$fluorescence))

df_tidy_ROS$fluorescence <- df_tidy_ROS$fluorescence/(max_possible_fluorescence_value)
df_tidy_ROS$Experiment <- as.factor(df_tidy_ROS$Experiment)
df_tidy_ROS$Well <- as.factor(df_tidy_ROS$Well)

time_points <- unique(df_tidy_ROS$TimePoint)
experiment_groups <- unique(df_tidy_ROS$Group)
experiments <- unique(df_tidy_ROS$Experiment)
wells <- unique(df_tidy_ROS$Well)

# Save final tibble ########################################################

df_final_1 <- df_tidy_ROS %>%
  filter(Group != "pos") %>%
  group_by(TimePoint, Group) %>%
  summarize(mean_fluorescence = mean(fluorescence, na.rm = TRUE),
            sd_fluorescence = sd(fluorescence, na.rm = TRUE))

df_final_2 <- df_final_1 %>%
  group_by(TimePoint) %>%
  summarize(ratio_mean_fluorescence = mean_fluorescence[Group == "stim"]/
              mean_fluorescence[Group == "control"],
            sum_sd_fluorescence = sd_fluorescence[Group == "stim"]+
              sd_fluorescence[Group == "control"])

df_final_2$TimePoint <- gsub(pattern = "h", replacement = "", x = df_final_2$TimePoint)
df_final_2$TimePoint <- as.numeric(df_final_2$TimePoint)

# Plotting #################################################################

dir.create(path = output_dir, showWarnings = FALSE)


# Plot for every time point and experiment class and experiment a histogram
# showing the three wells
for(i in time_points){
  for(j in experiment_groups){
    for(k in experiments){
      
      df_dummy <- df_tidy_ROS %>%
        filter(TimePoint == i) %>%
        filter(Group == j) %>%
        filter(Experiment == k)
      
      df_group_mean <- df_dummy %>%
        group_by(Well) %>%
        summarize(grp_mean = mean(fluorescence, na.rm=TRUE))
      
      
      plot_histogram <- ggplot(df_dummy, aes(x=fluorescence,
                                             fill = Well, color = Well)) +
        geom_histogram(binwidth = 0.003, alpha=0.2, position="identity") +
        geom_vline(data = df_group_mean,
                   mapping = aes(xintercept=grp_mean, color=Well),
                   linetype="dashed", size=1) +
        scale_x_continuous(expand = c(0, 0), limits = c(0, 1)) + 
        scale_y_continuous(expand = c(0, 0), limits = c(0, 50)) +
        labs(title=paste("Time point: ", i, ", Group: ", j, ", Exp: ", k, sep=""),
             x="Fluorescence intensity (arb. unit)", y = "Count")+
        theme_bw()
      
      # print(plot_histogram)
      
      file_name <- paste(output_dir, "/histogram_", i, "_", j, "_", k, ".png", sep="")
      ggsave(filename = file_name, width = 297, height = 210, units = "mm")
      
    }
  }
}

# Plot for every time point and experiment class a histogram showing the
# two independent experiments
for(i in time_points){
  for(j in experiment_groups){
    df_dummy <- df_tidy_ROS %>%
      filter(TimePoint == i) %>%
      filter(Group == j)
    
    df_group_mean <- df_dummy %>%
      group_by(Experiment) %>%
      summarize(grp_mean = mean(fluorescence, na.rm=TRUE))
    
    
    plot_histogram <- ggplot(df_dummy, aes(x=fluorescence,
                                           fill = Experiment, color = Experiment)) +
      geom_histogram(binwidth = 0.003, alpha=0.2, position="identity") +
      geom_vline(data = df_group_mean,
                 mapping = aes(xintercept=grp_mean, color=Experiment),
                 linetype="dashed", size=1) +
      scale_x_continuous(expand = c(0, 0), limits = c(0, 1)) + 
      scale_y_continuous(expand = c(0, 0), limits = c(0, 100)) +
      labs(title=paste("Time point: ", i, ", Group: ", j, sep=""),
           x="Fluorescence intensity (arb. unit)", y = "Count")+
      theme_bw()
    
    # print(plot_histogram)
    
    file_name <- paste(output_dir, "/histogram_", i, "_", j, ".png", sep="")
    ggsave(filename = file_name, width = 297, height = 210, units = "mm")
    
  }
}


# Plot ratio of ROS over time

plot_ROS_ratio <- ggplot(df_final_2, aes(x=TimePoint,ratio_mean_fluorescence)) +
  geom_point(size = 2) +
  geom_line() +
  geom_vline(xintercept = 0, linetype="dotted", color = "blue", size=1) +
  geom_text(aes(x=0, label="stim on", y=1.5), hjust=-0.1, colour="blue", text=element_text(size=11)) +
  geom_vline(xintercept = 2, linetype="dotted", color = "blue", size=1) +
  geom_text(aes(x=2, label="stim off", y=1.5), hjust=-0.1, colour="blue", text=element_text(size=11)) +
  geom_errorbar(aes(ymin=ratio_mean_fluorescence-sum_sd_fluorescence,
                    ymax=ratio_mean_fluorescence+sum_sd_fluorescence), width=.2) +
  labs(title="ROS ratio (stim vs. control)",
       x="time in h", y = "Ratio") +
  scale_x_continuous(expand = c(0, 0), limits = c(-1, 25), minor_breaks = c(-1:25)) + 
  theme_bw()

# print(plot_ROS_ratio)
  
file_name <- paste(output_dir, "/ROS_ratio.png", sep="")
ggsave(filename = file_name, width = 297, height = 210, units = "mm")