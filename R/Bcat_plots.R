# Script for plotting histograms of bcat readout using FACS ++++++++++++++++
# Author: Kai Budde
# Created: 2022/01/05
# Last changed: 2022/01/09

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
sample_names <- df_exp_metadata$experiment[match(x = sample_names, table = df_exp_metadata$frame_name)]
sampleNames(FACS_data_bcat_complete) <- sample_names

FACS_data_bcat_complete@phenoData@data$name <- sample_names

# Plot forward vs. side scatter data for every image #######################

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

# alternative:
plot <- ggcyto::autoplot(FACS_data_bcat_complete, x= "FSC.H", y = "SSC.H") +
  ggcyto::ggcyto_par_set(limits = list(x = c(0, 1024), y= c(0, 1024)))

file_name <- paste(output_dir, "/scatter_plot_FSCH_SSCH_complete.pdf", sep="")

ggsave(filename = file_name, width = (2*297), height = (2*210), units = "mm")


# Filter with polyGate #####################################################

FACS_data_bcat_filtered <- flowCore::Subset(x = FACS_data_bcat_complete,
                                            subset = polyGate)

# Plot histograms ##########################################################
# FL1.H -> DCFDA-ROS
# FL2.H -> beta-catenin with Cy3 as secondary AB

plot <- ggcyto::autoplot(FACS_data_bcat_filtered, x= "FL2.H") +
  scale_x_log10()

file_name <- paste(output_dir, "/density_plot_FL2H_complete.pdf", sep="")
ggsave(filename = file_name, width = (2*297), height = (2*210), units = "mm")

FACS_data_bcat_filtered_transformed <- flowCore::transform(FACS_data_bcat_filtered, `FL2.H`=log(`FL2.H`))
plot <- ggcyto::autoplot(FACS_data_bcat_filtered_transformed, x= "FL2.H")

file_name <- paste(output_dir, "/density_plot_FL2H_log_complete.pdf", sep="")
ggsave(filename = file_name, width = (2*297), height = (2*210), units = "mm")

# calculate mean and sd. of fluorescence ###################################
# FL1.H -> DCFDA-ROS
# FL2.H -> beta-catenin with Cy3 as secondary AB

# Calculate mean and standard deviation (sd)
dummy_mean_FACS_filtered <- fsApply(FACS_data_bcat_filtered_transformed, each_col, mean)
df_mean_FACS_filtered <- as.tibble(dummy_mean_FACS_filtered)
df_mean_FACS_filtered$Exp <- rownames(dummy_mean_FACS_filtered)
df_mean_FACS_filtered <- df_mean_FACS_filtered %>%
  relocate(Exp)  
rm(dummy_mean_FACS_filtered)


dummy_sd_FACS_filtered <- fsApply(FACS_data_bcat_filtered_transformed, each_col, sd)
df_sd_FACS_filtered <- as.tibble(dummy_sd_FACS_filtered)
df_sd_FACS_filtered$Exp <- rownames(dummy_sd_FACS_filtered)
df_sd_FACS_filtered <- df_sd_FACS_filtered %>%
  relocate(Exp)  
rm(dummy_sd_FACS_filtered)

# Rename columns
names(df_mean_FACS_filtered)[-1] <- paste(names(df_mean_FACS_filtered)[-1],
                                          ".mean", sep="")
names(df_sd_FACS_filtered)[-1] <- paste(names(df_sd_FACS_filtered)[-1],
                                        ".sd", sep="")

# Combine both tibbles
df_FACS_filtered <- cbind(df_mean_FACS_filtered, df_sd_FACS_filtered[,-1])

# Add metadata to tibble
df_FACS_filtered$timepoint <- df_FACS_filtered$Exp
search_pattern <- "^(.{1,2})h_.+"
df_FACS_filtered$timepoint[
  !grepl(pattern = search_pattern,
         x = df_FACS_filtered$timepoint,
         ignore.case = TRUE)] <- NA
df_FACS_filtered$timepoint <- gsub(pattern = search_pattern,
                                   replacement = "\\1",
                                   x = df_FACS_filtered$timepoint,
                                   ignore.case = TRUE)

df_FACS_filtered$timepoint <- as.numeric(df_FACS_filtered$timepoint)


df_FACS_filtered$group <- df_FACS_filtered$Exp
search_pattern <- ".+h_([a-z]+)_.+"
df_FACS_filtered$group[
  !grepl(pattern = search_pattern,
         x = df_FACS_filtered$group,
         ignore.case = TRUE)] <- NA
df_FACS_filtered$group <- gsub(pattern = search_pattern,
                                   replacement = "\\1",
                                   x = df_FACS_filtered$group,
                                   ignore.case = TRUE)

df_FACS_filtered$well <- df_FACS_filtered$Exp
search_pattern <- ".+_([A-C]+)_.+"
df_FACS_filtered$well[
  !grepl(pattern = search_pattern,
         x = df_FACS_filtered$well,
         ignore.case = TRUE)] <- NA
df_FACS_filtered$well <- gsub(pattern = search_pattern,
                               replacement = "\\1",
                               x = df_FACS_filtered$well,
                               ignore.case = TRUE)

df_FACS_filtered$well[is.na(df_FACS_filtered$well)] <- "A"

df_FACS_filtered$replication <- df_FACS_filtered$Exp
search_pattern <- ".+_([1-2])$"
df_FACS_filtered$replication[
  !grepl(pattern = search_pattern,
         x = df_FACS_filtered$replication,
         ignore.case = TRUE)] <- NA
df_FACS_filtered$replication <- gsub(pattern = search_pattern,
                              replacement = "\\1",
                              x = df_FACS_filtered$replication,
                              ignore.case = TRUE)
df_FACS_filtered$replication <- as.numeric(df_FACS_filtered$replication)
df_FACS_filtered$replication[is.na(df_FACS_filtered$replication)] <- 1


# Plot mean and standard deviations ########################################

df_FACS_filtered_for_plotting <- df_FACS_filtered %>%
  dplyr::filter(!is.na(group)) %>%
  dplyr::filter(group != "pos")

plot <- ggplot(df_FACS_filtered_for_plotting, aes(x = timepoint, y = FL2.H.mean, color = group)) +
  geom_point(size = 2) +
  # geom_line() +
  # geom_vline(xintercept = 0, linetype="dotted", color = "blue", size=1) +
  # geom_text(aes(x=0, label="stim on", y=1.5), hjust=-0.1, colour="blue", text=element_text(size=11)) +
  # geom_vline(xintercept = 2, linetype="dotted", color = "blue", size=1) +
  # geom_text(aes(x=2, label="stim off", y=1.5), hjust=-0.1, colour="blue", text=element_text(size=11)) +
  geom_errorbar(aes(ymin=FL2.H.mean - FL2.H.sd,
                    ymax=FL2.H.mean + FL2.H.sd), width=.2) +
  labs(title="Fluorescence (log) (mean +- sd) of Beta-Catenin",
       x="time in h", y = "Fluorescence (arb. units)") +
  scale_x_continuous(expand = c(0, 0), limits = c(-1, 25), minor_breaks = c(-1:25)) + 
  theme_bw()

file_name <- paste(output_dir, "/mean_sd_plot_FL2H_log.png", sep="")
ggsave(filename = file_name, width = (297), height = (210), units = "mm")
