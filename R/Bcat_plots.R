# Script for plotting histograms of bcat readout using FACS ++++++++++++++++
# Author: Kai Budde
# Created: 2022/01/05
# Last changed: 2022/05/24

rm(list=ls())
Sys.setlocale("LC_TIME","English")

# General function parameters ##############################################
options(stringsAsFactors = FALSE, warn=-1)


# Load packages (Using Groundhog does not work properly because it cannot
# load from Bioconductor)

# # Set groundhog day for reproducibility (see https://groundhogr.com)
# groundhog.day <- "2022-03-01"
# 
# if(!any(grepl(pattern = "groundhog", x = installed.packages(), ignore.case = TRUE))){
#   install.packages("groundhog")
# }
# library(groundhog)
# pkgs <- c("tidyr","BiocManager", "rstatix", "vctrs", "flowCore")
# groundhog.library(pkgs, groundhog.day)

list.of.packages <- c("BiocManager", "tidyverse", "rstatix")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

require("BiocManager")
require("tidyverse")
require("rstatix")

# Install from BioConductor
if (!requireNamespace("flowCore", quietly = TRUE)){
  BiocManager::install("flowCore")
}
if (!requireNamespace("ggcyto", quietly = TRUE)){
  BiocManager::install("ggcyto")
}
if (!requireNamespace("flowViz", quietly = TRUE)){
  BiocManager::install("flowViz")
}

require("flowCore")
require("ggcyto")
require("flowViz")

# Siehe https://www.bioconductor.org/packages/devel/bioc/vignettes/flowViz/inst/doc/filters.html


# Input parameters #########################################################

facs_data_dir <- "data/beta_catenin_raw_data/"
facs_metadata <- "data/beta_catenin_data_assignments.csv"
output_dir <- "plots/bcat"
save_scatter_plots <- FALSE
save_all_scatter_plots_in_one <- FALSE
save_histograms <- FALSE
use_log_values <- FALSE

# Filter
polyGate <- flowCore::polygonGate(data.frame(
  "FSC.H" = c(400, 800, 800,  400,  100),
  "SSC.H" = c(0,   200, 1024, 1024, 400)))

# Data input, cleaning, and additional columns #############################

input_dir <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(dir = input_dir)
setwd(dir = "..")

# Read assignments to FACS data
df_exp_metadata <- readr::read_csv(file = facs_metadata, name_repair = "universal")

# Read all FACS data
#description(read.FCS("data/beta_catenin_raw_data/b catenin.001", emptyValue=FALSE))
FACS_data_bcat_complete <- flowCore::read.flowSet(path = facs_data_dir, alter.names = TRUE)

# Rename the frames
sample_names <- sampleNames(FACS_data_bcat_complete)
sample_names <- df_exp_metadata$Experiment[match(x = sample_names, table = df_exp_metadata$Frame_name)]
sampleNames(FACS_data_bcat_complete) <- sample_names

FACS_data_bcat_complete@phenoData@data$name <- sample_names

dates <- flowCore::keyword(object = FACS_data_bcat_complete, keyword = "$DATE")
times <- flowCore::keyword(object = FACS_data_bcat_complete, keyword = "$ETIM")


# Plot forward vs. side scatter data for every image #######################

if(save_scatter_plots){
  dir_scatter_plots <- paste(output_dir, "/scatter", sep="")
  dir.create(path = dir_scatter_plots, showWarnings = FALSE)
  for(i in sample_names){
    plot <- flowViz::xyplot(x = `SSC.H` ~ `FSC.H`,
                            data = FACS_data_bcat_complete[[i]],
                            filter = polyGate, smooth = FALSE,
                            main=paste("Experiment: ", i, sep=""),
                            ylab="Side scatter (H)",
                            xlab="Forward scatter (H)")
    
    file_name <- paste(dir_scatter_plots, "/scatter_plot_exp_", i, ".png", sep="")
    
    png(filename=file_name, type="cairo", units="mm",
        width=297, height=210, pointsize=12, res=300)
    print(plot)
    dev.off()
    
  }
}

# alternative:
if(save_all_scatter_plots_in_one){
  plot <- ggcyto::autoplot(FACS_data_bcat_complete, x= "FSC.H", y = "SSC.H") +
    ggcyto::ggcyto_par_set(limits = list(x = c(0, 1024), y= c(0, 1024)))
  
  file_name <- paste(output_dir, "/scatter_plot_FSCH_SSCH_all.pdf", sep="")
  
  ggsave(filename = file_name, width = (2*297), height = (2*210), units = "mm")
}

# Filter with polyGate #####################################################

FACS_data_bcat_filtered <- flowCore::Subset(x = FACS_data_bcat_complete,
                                            subset = polyGate)

# Transform (some) data to log10 base
log10transform <- logTransform(transformationId="log10-transformation", logbase=10, r=1, d=1)

FACS_data_bcat_filtered_transformed <- flowCore::transform(
  FACS_data_bcat_filtered, transformList(c("FL2.H", "FSC.H", "SSC.H"),
                                         list(log10transform, log10transform, log10transform)))

# Plot histograms ##########################################################
# FL1.H -> DCFDA-ROS
# FL2.H -> beta-catenin with Cy3 as secondary AB
# Forward scatter (FSC.H): Size of cell
# Sideward scatter (SSC.H: Granularity of cell

if(save_histograms){
  
  # FL2.H (Beta-Catenin)
  plot <- ggcyto::autoplot(FACS_data_bcat_filtered, x= "FL2.H", arrange = FALSE) +
    ggplot2::ggtitle("Fluorescence intensity (height) of Beta-Catenin with Cy3 as secondary AB")
  
  file_name <- paste(output_dir, "/density_plot_FL2H_all.pdf", sep="")
  ggsave(filename = file_name, width = (2*297), height = (2*210), units = "mm")
  
  plot <- ggcyto::autoplot(FACS_data_bcat_filtered, x= "FL2.H") +
    scale_x_log10() +
    ggplot2::ggtitle("Fluorescence intensity (height) of Beta-Catenin with Cy3 as secondary AB (log axis)")

  file_name <- paste(output_dir, "/density_plot_FL2H_all_log_axis.pdf", sep="")
  ggsave(filename = file_name, width = (2*297), height = (2*210), units = "mm")
  
  plot <- ggcyto::autoplot(FACS_data_bcat_filtered_transformed, x= "FL2.H") +
    ggplot2::ggtitle("Fluorescence intensity (height) of Beta-Catenin with Cy3 as secondary AB (log data)")
  
  file_name <- paste(output_dir, "/density_plot_FL2H_all_log_data.pdf", sep="")
  ggsave(filename = file_name, width = (2*297), height = (2*210), units = "mm")
  
  
  dir_density_plos <- paste(output_dir, "/density", sep="")
  dir.create(path = dir_density_plos, showWarnings = FALSE)
  
  for(i in sample_names){
    
    plot <- suppressMessages(ggcyto::autoplot(FACS_data_bcat_filtered_transformed[[i]], x= "FL2.H") +
                               ggtitle("Fluorescence intensity (height) of Beta-Catenin with Cy3 as secondary AB (log data)") +
                               coord_cartesian(xlim = c(0, 3), ylim = c(0.0, 4)) )

    file_name <- paste(dir_density_plos, "/FL2H_log_density_exp_", i, ".png", sep="")
    
    png(filename=file_name, type="cairo", units="mm",
        width=297, height=210, pointsize=12, res=300)
    print(plot)
    dev.off()
    
  }
  
  # FSC.H (Forward scatter)
  plot <- ggcyto::autoplot(FACS_data_bcat_filtered, x= "FSC.H") +
    ggplot2::ggtitle("Forward scatter intensity (height)")
  
  file_name <- paste(output_dir, "/density_plot_FSCH_all.pdf", sep="")
  ggsave(filename = file_name, width = (2*297), height = (2*210), units = "mm")
  
  
  for(i in sample_names){

    plot <- suppressMessages(ggcyto::autoplot(FACS_data_bcat_filtered[[i]], x= "FSC.H") +
      coord_cartesian(xlim = c(0, 1024), ylim = c(0.0, 5e-3)))
    
    file_name <- paste(dir_density_plos, "/FSCH_density_exp_", i, ".png", sep="")
    
    png(filename=file_name, type="cairo", units="mm",
        width=297, height=210, pointsize=12, res=300)
    print(plot)
    dev.off()
    
  }
  
  plot <- ggcyto::autoplot(FACS_data_bcat_filtered, x= "FSC.H") +
    scale_x_log10() +
    ggplot2::ggtitle("Forward scatter intensity (height) (log axis)")
  
  file_name <- paste(output_dir, "/density_plot_FSCH_all_log_axis.pdf", sep="")
  ggsave(filename = file_name, width = (2*297), height = (2*210), units = "mm")
  
  plot <- ggcyto::autoplot(FACS_data_bcat_filtered_transformed, x= "FSC.H") +
    ggplot2::ggtitle("Forward scatter intensity (height) (log data)")
  
  file_name <- paste(output_dir, "/density_plot_FSCH_all_log_data.pdf", sep="")
  ggsave(filename = file_name, width = (2*297), height = (2*210), units = "mm")
  
  
  # SSC.H (Sideward scatter)
  plot <- ggcyto::autoplot(FACS_data_bcat_filtered, x= "SSC.H") +
    ggplot2::ggtitle("Sideward scatter intensity (height)")
  
  file_name <- paste(output_dir, "/density_plot_SSCH_all.pdf", sep="")
  ggsave(filename = file_name, width = (2*297), height = (2*210), units = "mm")
  
  for(i in sample_names){
    
    plot <- suppressMessages(ggcyto::autoplot(FACS_data_bcat_filtered[[i]], x= "SSC.H") +
                               coord_cartesian(xlim = c(0, 1024), ylim = c(0.0, 5e-3)))
    
    file_name <- paste(dir_density_plos, "/SSCH_density_exp_", i, ".png", sep="")
    
    png(filename=file_name, type="cairo", units="mm",
        width=297, height=210, pointsize=12, res=300)
    print(plot)
    dev.off()
    
  }
  
  plot <- ggcyto::autoplot(FACS_data_bcat_filtered, x= "SSC.H") +
    scale_x_log10() +
    ggplot2::ggtitle("Sideward scatter intensity (height) (log axis)")
  
  file_name <- paste(output_dir, "/density_plot_SSCH_all_log_axis.pdf", sep="")
  ggsave(filename = file_name, width = (2*297), height = (2*210), units = "mm")
  
  plot <- ggcyto::autoplot(FACS_data_bcat_filtered_transformed, x= "SSC.H") +
    ggplot2::ggtitle("Sideward scatter intensity (height) (log data)")
  
  file_name <- paste(output_dir, "/density_plot_SSCH_all_log_data.pdf", sep="")
  ggsave(filename = file_name, width = (2*297), height = (2*210), units = "mm")
  
}


# calculate mean and sd. of fluorescence ###################################
# FL1.H -> DCFDA-ROS
# FL2.H -> beta-catenin with Cy3 as secondary AB

# Calculate mean and standard deviation (sd)
if(use_log_values){
  dummy_mean_FACS_filtered <- flowCore::fsApply(FACS_data_bcat_filtered_transformed, each_col, mean)
}else{
  dummy_mean_FACS_filtered <- flowCore::fsApply(FACS_data_bcat_filtered, each_col, mean)
}
# dummy_mean_FACS_filtered <- flowCore::fsApply(FACS_data_bcat_filtered_transformed, each_col, mean)
df_mean_FACS_filtered <- as.tibble(dummy_mean_FACS_filtered)
df_mean_FACS_filtered$ID <- rownames(dummy_mean_FACS_filtered)
df_mean_FACS_filtered <- df_mean_FACS_filtered %>%
  relocate(ID)  
rm(dummy_mean_FACS_filtered)

if(use_log_values){
  dummy_sd_FACS_filtered <- fsApply(FACS_data_bcat_filtered_transformed, each_col, sd)
  
}else{
  dummy_sd_FACS_filtered <- fsApply(FACS_data_bcat_filtered, each_col, sd)
}
# dummy_sd_FACS_filtered <- fsApply(FACS_data_bcat_filtered_transformed, each_col, sd)
df_sd_FACS_filtered <- as.tibble(dummy_sd_FACS_filtered)
df_sd_FACS_filtered$ID <- rownames(dummy_sd_FACS_filtered)
df_sd_FACS_filtered <- df_sd_FACS_filtered %>%
  relocate(ID)  
rm(dummy_sd_FACS_filtered)

# Rename columns
names(df_mean_FACS_filtered)[-1] <- paste(names(df_mean_FACS_filtered)[-1],
                                          ".mean", sep="")
names(df_sd_FACS_filtered)[-1] <- paste(names(df_sd_FACS_filtered)[-1],
                                        ".sd", sep="")

# Combine both tibbles
df_FACS_filtered <- cbind(df_mean_FACS_filtered, df_sd_FACS_filtered[,-1])

# Add metadata to tibble
df_FACS_filtered$TimePoint <- df_FACS_filtered$ID
search_pattern <- "^(.{1,2})h_.+"
df_FACS_filtered$TimePoint[
  !grepl(pattern = search_pattern,
         x = df_FACS_filtered$TimePoint,
         ignore.case = TRUE)] <- NA
df_FACS_filtered$TimePoint <- gsub(pattern = search_pattern,
                                   replacement = "\\1",
                                   x = df_FACS_filtered$TimePoint,
                                   ignore.case = TRUE)

df_FACS_filtered$TimePoint <- as.numeric(df_FACS_filtered$TimePoint)


df_FACS_filtered$Group <- df_FACS_filtered$ID
search_pattern <- ".+h_([a-z]+)_.+"
df_FACS_filtered$Group[
  !grepl(pattern = search_pattern,
         x = df_FACS_filtered$Group,
         ignore.case = TRUE)] <- NA
df_FACS_filtered$Group <- gsub(pattern = search_pattern,
                                   replacement = "\\1",
                                   x = df_FACS_filtered$Group,
                                   ignore.case = TRUE)

df_FACS_filtered$Well <- df_FACS_filtered$ID
search_pattern <- ".+_([A-B][1-3])_.+"
df_FACS_filtered$Well[
  !grepl(pattern = search_pattern,
         x = df_FACS_filtered$Well,
         ignore.case = TRUE)] <- NA
df_FACS_filtered$Well <- gsub(pattern = search_pattern,
                               replacement = "\\1",
                               x = df_FACS_filtered$Well,
                               ignore.case = TRUE)
df_FACS_filtered$WellID <- as.numeric(gsub(pattern = ".+([1-3])",
                                           replacement = "\\1",
                                           x = df_FACS_filtered$Well))

df_FACS_filtered$Well[is.na(df_FACS_filtered$Well)] <- "Z"

df_FACS_filtered$Experiment <- df_FACS_filtered$ID
search_pattern <- ".+_([1-2])$"
df_FACS_filtered$Experiment[
  !grepl(pattern = search_pattern,
         x = df_FACS_filtered$Experiment,
         ignore.case = TRUE)] <- NA
df_FACS_filtered$Experiment <- gsub(pattern = search_pattern,
                              replacement = "\\1",
                              x = df_FACS_filtered$Experiment,
                              ignore.case = TRUE)
df_FACS_filtered$Experiment <- as.numeric(df_FACS_filtered$Experiment)

df_FACS_filtered$Experiment[is.na(df_FACS_filtered$Experiment)] <- 99

# Calibrate data ###########################################################

# First measurements
calib_fluorescence <- df_FACS_filtered$FL2.H.mean[df_FACS_filtered$ID == "Blank--+_1"]

df_FACS_filtered$FL2.H.mean[df_FACS_filtered$Experiment != 99 &
                              !is.na(df_FACS_filtered$Experiment) &
                              df_FACS_filtered$TimePoint != 36 &
                              !is.na(df_FACS_filtered$TimePoint)] <-
  df_FACS_filtered$FL2.H.mean[df_FACS_filtered$Experiment != 99 &
                                !is.na(df_FACS_filtered$Experiment) &
                                df_FACS_filtered$TimePoint != 36 &
                                !is.na(df_FACS_filtered$TimePoint)] - calib_fluorescence


calib_fluorescence <- df_FACS_filtered$FL2.H.mean[df_FACS_filtered$ID == "Blank--+_2"]

df_FACS_filtered$FL2.H.mean[(df_FACS_filtered$TimePoint == 36 &
                               !is.na(df_FACS_filtered$TimePoint) ) |
                              df_FACS_filtered$Experiment == 99] <- 
  df_FACS_filtered$FL2.H.mean[(df_FACS_filtered$TimePoint == 36 &
                                 !is.na(df_FACS_filtered$TimePoint) ) |
                                df_FACS_filtered$Experiment == 99] -
  calib_fluorescence


# Calculate ratios of stim and control group ###############################

df_final <-  df_FACS_filtered %>%
  group_by(TimePoint, WellID, Experiment) %>%
  summarize(ratio_mean_fluorescence = FL2.H.mean[Group == "stim"]/
              FL2.H.mean[Group == "control"],
            sum_sd_fluorescence = FL2.H.sd[Group == "stim"]+
              FL2.H.sd[Group == "control"])

df_final <- df_final[!is.na(df_final$ratio_mean_fluorescence),]


# Remove outliers from new measurements
df_final <- df_final[df_final$Experiment != 99, ]

# Add empty timepoints
df_final <- as.data.frame(df_final)

for(i in 0:(max(df_final$TimePoint)+1)){
  if(!(i %in% df_final$TimePoint)){
    df_final <- rbind(df_final, c(i, rep(NA, length(names(df_final))-1) ))
  }
}

# Save time points as factors
df_final$TimePoint <- as.factor(df_final$TimePoint)

# Save aggregates data as csv ##############################################
df_final_output <- df_FACS_filtered %>%
  select(ID, FL1.H.mean, FL2.H.mean, TimePoint, Group, Well, Experiment) %>%
  dplyr::filter(Group != "pos")

names(df_final_output)[names(df_final_output) == "FL1.H.mean"] <- "DCF_ROS_mean"
names(df_final_output)[names(df_final_output) == "FL2.H.mean"] <- "Beta_catenin_mean"

df_final_output$Time <- times[match(df_final_output$ID, dimnames(times)[[1]])]
df_final_output$Date <- dates[match(df_final_output$ID, dimnames(dates)[[1]])]
df_final_output$Date <- as.Date(df_final_output$Date, format = "%d-%b-%y")

file_name <- paste(output_dir, "/betacatenin_data.csv", sep="")
write.csv(x = df_final_output, file = file_name, row.names = FALSE)

# Plot ratios of mean fluorescence intensities #############################

plot_ratios <- ggplot(df_final, aes(x=TimePoint, y=ratio_mean_fluorescence)) +
  stat_boxplot(geom ='errorbar', width = 0.3, position = position_dodge(width = 0.75)) +
  geom_boxplot(alpha = 1, position = position_dodge2(width = 0.9, preserve = "single"), outlier.shape = 1) +
  geom_jitter(color="black", size=0.5, alpha=0.9) +
  # geom_point(position = position_jitterdodge(), size=1, alpha=0.9) +
  ylim(0,max(df_final$ratio_mean_fluorescence)) +
  theme_bw(base_size = 18) +
  geom_hline(yintercept=1, linetype="dashed", color = "black", size=1) +
  theme(#axis.title.y=element_text(size=12),
    #axis.text.x = element_blank(), 
    axis.ticks.x = element_blank()) +
  ylab("Beta-catenin ratio of mean fluorescence (stim vs. control)") +
  xlab("Time in h")

file_name <- paste(output_dir, "/Beta_catenin_FACS_ratio.png", sep="")
ggsave(filename = file_name, width = (297), height = (210), units = "mm")


# Plot mean fluorescence intensities #######################################

df_mean_data <-  df_FACS_filtered %>%
  group_by(TimePoint, WellID, Experiment)

# Keep only stim and control
df_mean_data <- as.data.frame(df_mean_data)

# Remove outliers from new measurements
df_mean_data <- df_mean_data[df_mean_data$Experiment != 99, ]

# Keep only relevant columns
names(df_mean_data)[names(df_mean_data) == "FL1.H.mean"] <- "DCF_ROS_mean"
names(df_mean_data)[names(df_mean_data) == "FL2.H.mean"] <- "Beta_catenin_mean"

df_mean_data <- df_mean_data %>%
  select(TimePoint, ID, Beta_catenin_mean, Group, Well, WellID, Experiment) %>%
  dplyr::filter(Group != "pos")

# Add time points
for(i in 0:(max(df_mean_data$TimePoint)+1)){
  if(!(i %in% df_mean_data$TimePoint)){
    df_mean_data <- rbind(df_mean_data, c(i, rep(NA, length(names(df_mean_data))-1) ))
  }
}

# Save time points as factors
# df_mean_data$TimePoint <- as.factor(df_mean_data$TimePoint)

plot_bcat_boxplots <- ggplot(df_mean_data, aes(x = as.factor(df_mean_data$TimePoint),
                                               y = Beta_catenin_mean, fill = Group)) +
  stat_boxplot(geom ='errorbar', width = 0.3, position = position_dodge(width = 0.75)) +
  geom_boxplot(alpha = 1, position = position_dodge2(width = 0.9, preserve = "single"), outlier.shape = 1) +
  geom_jitter(color="black", size=0.5, alpha=0.9, position = position_dodge2(width = 0.9, preserve = "single")) +
  theme_bw(base_size = 18) +
  theme(#axis.title.y=element_text(size=12),
    #axis.text.x = element_blank(), 
    axis.ticks.x = element_blank()) +
  ylab("Bcat FACS mean intensity") +
  xlab("Time in h")

file_name <- paste(output_dir, "/Beta_catenin_FACS_means.png", sep="")
ggsave(filename = file_name, width = (297), height = (210), units = "mm")


# Calculate paired t-test

for(t in c(2,4,8,24,36)){
  df_dummy <- df_mean_data[df_mean_data$TimePoint == t,]
  results <- rstatix::pairwise_t_test(data = df_dummy, formula = Beta_catenin_mean ~ Group)
  print(t)
  print(results$p)
  
}

# Plot means by WellID #####################################################

df_mean_data <- df_mean_data[!is.na(df_mean_data$Beta_catenin_mean),]

df_mean_data$WellID_Exp <- as.factor(paste(df_mean_data$WellID, "_", df_mean_data$Experiment, sep=""))
df_mean_data$Experiment <- as.factor(df_mean_data$Experiment)

my_colors <- c("#1b9e77", "#d95f02", "#7570b3", "#e7298a", "#66a61e", "#e6ab02")

plot_bcat_WellID <- ggplot(df_mean_data,
                           aes(x=TimePoint, y=Beta_catenin_mean,
                               color=WellID_Exp, shape = Group)) +
  geom_line(size = 1) +
  geom_point(size = 6) +
  scale_color_manual(values = my_colors) +
  ylim(0,max(df_mean_data$Beta_catenin_mean)) +
  theme_bw(base_size = 24) +
  ylab("Bcat FACS mean intensity") +
  xlab("Time in h")

file_name <- paste(output_dir, "/Beta_catenin_FACS_means_per_WellExp.png", sep="")
ggsave(filename = file_name, width = (297), height = (210), units = "mm")






# # Plot mean and standard deviations ########################################
# 
# df_FACS_filtered_for_plotting <- df_FACS_filtered %>%
#   dplyr::filter(!is.na(Group)) %>%
#   dplyr::filter(Group != "pos")
# 
# plot <- ggplot(df_FACS_filtered_for_plotting, aes(x = TimePoint, y = FL2.H.mean, color = Group)) +
#   geom_point(size = 2) +
#   # geom_line() +
#   # geom_vline(xintercept = 0, linetype="dotted", color = "blue", size=1) +
#   # geom_text(aes(x=0, label="stim on", y=1.5), hjust=-0.1, colour="blue", text=element_text(size=11)) +
#   # geom_vline(xintercept = 2, linetype="dotted", color = "blue", size=1) +
#   # geom_text(aes(x=2, label="stim off", y=1.5), hjust=-0.1, colour="blue", text=element_text(size=11)) +
#   geom_errorbar(aes(ymin=FL2.H.mean - FL2.H.sd,
#                     ymax=FL2.H.mean + FL2.H.sd), width=.2) +
#   labs(title="Fluorescence (log) (mean +- sd) of Beta-Catenin",
#        x="time in h", y = "Fluorescence (arb. units)") +
#   scale_x_continuous(expand = c(0, 0), limits = c(-1, 25), minor_breaks = c(-1:25)) + 
#   theme_bw()
# 
# file_name <- paste(output_dir, "/mean_sd_plot_FL2H_log.png", sep="")
# ggsave(filename = file_name, width = (297), height = (210), units = "mm")
