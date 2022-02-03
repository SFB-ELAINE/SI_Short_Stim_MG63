# Script for plotting histograms of bcat readout using FACS ++++++++++++++++
# Author: Kai Budde
# Created: 2022/01/05
# Last changed: 2022/02/03

rm(list=ls())
Sys.setlocale("LC_TIME","English")

# General function parameters ##############################################
options(stringsAsFactors = FALSE, warn=-1)

if (!requireNamespace("BiocManager", quietly = TRUE)){
  install.packages("BiocManager")
}
if (!requireNamespace("flowCore", quietly = TRUE)){
  BiocManager::install("flowCore")
}
if (!requireNamespace("ggcyto", quietly = TRUE)){
  install.packages("ggcyto")
}
if (!requireNamespace("flowViz", quietly = TRUE)){
  BiocManager::install("flowViz")
}

library("tidyverse")
library("ggplot2")
library("ggcyto")
library("flowCore")
library("flowViz")

# Siehe https://www.bioconductor.org/packages/devel/bioc/vignettes/flowViz/inst/doc/filters.html


# Input parameters #########################################################

facs_data_dir <- "data/beta_catenin_raw_data/"
facs_metadata <- "data/beta_catenin_data_assignments.csv"
output_dir <- "plots"
save_scatter_plots <- FALSE
save_all_scatter_plots_in_one <- FALSE
save_histograms <- FALSE

# Filter
polyGate <- flowCore::polygonGate(data.frame(
  "FSC.H" = c(400, 800, 800,  400,  100),
  "SSC.H" = c(0,   200, 1024, 1024, 400)))

# Data input, cleaning, and additional columns #############################

input_dir <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(dir = input_dir)
setwd(dir = "..")

# Read assignments to FACS data
df_exp_metadata <- read_csv(file = facs_metadata, name_repair = "universal")

# Read all FACS data
FACS_data_bcat_complete <- flowCore::read.flowSet(path = facs_data_dir, alter.names = TRUE)

# Rename the frames
sample_names <- sampleNames(FACS_data_bcat_complete)
sample_names <- df_exp_metadata$Experiment[match(x = sample_names, table = df_exp_metadata$Frame_name)]
sampleNames(FACS_data_bcat_complete) <- sample_names

FACS_data_bcat_complete@phenoData@data$name <- sample_names

# Plot forward vs. side scatter data for every image #######################

if(save_scatter_plots){
  for(i in sample_names){
    plot <- flowViz::xyplot(x = `SSC.H` ~ `FSC.H`,
                            data = FACS_data_bcat_complete[[i]],
                            filter = polyGate, smooth = FALSE,
                            main=paste("Experiment: ", i, sep=""),
                            ylab="Side scatter (H)",
                            xlab="Forward scatter (H)")
    
    file_name <- paste(output_dir, "/scatter_plot_exp_", i, ".png", sep="")
    
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
  
  file_name <- paste(output_dir, "/scatter_plot_FSCH_SSCH_complete.pdf", sep="")
  
  ggsave(filename = file_name, width = (2*297), height = (2*210), units = "mm")
}



# Filter with polyGate #####################################################

FACS_data_bcat_filtered <- flowCore::Subset(x = FACS_data_bcat_complete,
                                            subset = polyGate)

# Plot histograms ##########################################################
# FL1.H -> DCFDA-ROS
# FL2.H -> beta-catenin with Cy3 as secondary AB

FACS_data_bcat_filtered_transformed <- flowCore::transform(FACS_data_bcat_filtered, `FL2.H`=log(`FL2.H`))

if(save_histograms){
  plot <- ggcyto::autoplot(FACS_data_bcat_filtered, x= "FL2.H") +
    scale_x_log10()
  
  file_name <- paste(output_dir, "/density_plot_FL2H_complete.pdf", sep="")
  ggsave(filename = file_name, width = (2*297), height = (2*210), units = "mm")
  
  plot <- ggcyto::autoplot(FACS_data_bcat_filtered_transformed, x= "FL2.H")
  
  file_name <- paste(output_dir, "/density_plot_FL2H_log_complete.pdf", sep="")
  ggsave(filename = file_name, width = (2*297), height = (2*210), units = "mm")
}


# calculate mean and sd. of fluorescence ###################################
# FL1.H -> DCFDA-ROS
# FL2.H -> beta-catenin with Cy3 as secondary AB

# Calculate mean and standard deviation (sd)
dummy_mean_FACS_filtered <- flowCore::fsApply(FACS_data_bcat_filtered, each_col, mean)
# dummy_mean_FACS_filtered <- flowCore::fsApply(FACS_data_bcat_filtered_transformed, each_col, mean)
df_mean_FACS_filtered <- as.tibble(dummy_mean_FACS_filtered)
df_mean_FACS_filtered$ID <- rownames(dummy_mean_FACS_filtered)
df_mean_FACS_filtered <- df_mean_FACS_filtered %>%
  relocate(ID)  
rm(dummy_mean_FACS_filtered)


dummy_sd_FACS_filtered <- fsApply(FACS_data_bcat_filtered, each_col, sd)
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
search_pattern <- ".+_([A-C]+)_.+"
df_FACS_filtered$Well[
  !grepl(pattern = search_pattern,
         x = df_FACS_filtered$Well,
         ignore.case = TRUE)] <- NA
df_FACS_filtered$Well <- gsub(pattern = search_pattern,
                               replacement = "\\1",
                               x = df_FACS_filtered$Well,
                               ignore.case = TRUE)

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
  group_by(TimePoint, Well, Experiment) %>%
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

# Plot ratios of mean fluorescence intensities #############################

plot_ratios <- ggplot(df_final, aes(x=TimePoint, y=ratio_mean_fluorescence)) +
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

file_name <- paste(output_dir, "/Beta_catenin_ratio.png", sep="")
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
