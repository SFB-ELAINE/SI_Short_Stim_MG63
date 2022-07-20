# Script for analyzing czi images of ROX                           +++++++++
# Author: Kai Budde
# Created: 2022/06/30
# Last changed: 2022/07/19

# Delete everything in the environment
rm(list = ls())
# close all open plots in RStudio
graphics.off()

# Please adapt the following parameters ####################################
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

# Directory of the images with CellROX measurements
input_directories_CellRox <- c(
  # "E:/PhD/Daten/ShortStim_ZellBio/CellRox/09.05.2022/",
  "E:/PhD/Daten/ShortStim_ZellBio/CellRox/Versuch 1 220712/",
  "E:/PhD/Daten/ShortStim_ZellBio/CellRox/Versuch 2 220713/",
  "E:/PhD/Daten/ShortStim_ZellBio/CellRox/Versuch 3 220715/")

output_dir <- "E:/PhD/Daten/ShortStim_ZellBio/CellRox/"

# The the following TRUE if the czi images have not yet been analyzed
analyzeCziImages <- FALSE

# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

# Load packages #############################################################

list.of.packages <- c("BiocManager", "devtools", "tidyverse")

new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

if(!("EBImage" %in% installed.packages())){
  print("Installing EBImage.")
  BiocManager::install("EBImage")
}
require(EBImage)
require(devtools)
require(tidyverse)

# Install readCzi if not already installed (used version: 0.1.??)
devtools::install_github(repo = "https://github.com/SFB-ELAINE/readCzi")#, ref = "v0.1.??")
require(readCzi)


# Install cellPixels if not already installed (used version: 0.1.??)
devtools::install_github(repo = "https://github.com/SFB-ELAINE/cellPixels")#, ref = "v0.1.??")
require(cellPixels)


## Check for CellROX -------------------------------------------------------
if(analyzeCziImages){
  input_directories <- c(input_directories_CellRox)
  
  for(i in 1:length(input_directories)){
    if(i == 1){
      df_results <- cellPixels::cellPixels(input_dir = input_directories[i],
                                           nucleus_color = "blue",
                                           protein_in_nuc_color = "green",
                                           protein_in_cytosol_color = "red",
                                           number_size_factor = 0.2,
                                           add_scale_bar = TRUE)
    }else{
      df_dummy <- cellPixels::cellPixels(input_dir = input_directories[i],
                                         nucleus_color = "blue",
                                         protein_in_nuc_color = "red",
                                         protein_in_cytosol_color = "green",
                                         number_size_factor = 0.2,
                                         add_scale_bar = TRUE)
      df_results <- rbind(df_results, df_dummy)
    }
    
  }
  rm(i)
  
  # Get and save metadata as csv
  input_directories <- c(input_directories_CellRox)
  
  for(i in 1:length(input_directories)){
    
    current_input_dir <- input_directories[i]
    # Data frame with certain parameter values
    file_names <- list.files(path = current_input_dir)
    file_names_czi <- file_names[grepl("czi", file_names)]
    file_names_czi <- paste(current_input_dir, file_names_czi, sep="/")
    number_of_czi_files <- length(file_names_czi)
    
    
    for(j in 1:number_of_czi_files){
      print(j)
      if(j==1){
        df_metadata <- readCzi::readCziMetadata(input_file = file_names_czi[j], save_metadata = FALSE)
      }else{
        df_dummy <- readCzi::readCziMetadata(input_file = file_names_czi[j], save_metadata = FALSE)
        df_metadata <- dplyr::bind_rows(df_metadata, df_dummy)
        rm(df_dummy)
      }
    }
    rm(j)
    
    current_output_dir <- paste0(current_input_dir, "output/")
    dir.create(current_output_dir, showWarnings = FALSE)
    
    write.csv(x = df_metadata,
              file = paste(output_dir,"summary_metadata.csv", sep=""),
              row.names = FALSE)
    write.csv2(x = df_metadata,
               file = paste(output_dir,"summary_metadata_de.csv", sep=""),
               row.names = FALSE)
    
    
  }
  rm(i)
  
}

## Data analysis -----------------------------------------------------------

# Import data
input_directories <- c(input_directories_CellRox)
for(i in 1:length(input_directories)){
  if(i == 1){
    df_data <- readr::read_csv(file = paste0(input_directories[i], "/output/image_analysis_summary_en.csv"), name_repair = "universal")
    df_metadata  <- readr::read_csv(file = paste0(input_directories[i], "/output/summary_metadata.csv"), name_repair = "universal")
  }else{
    df_data_dummy <- readr::read_csv(file = paste0(input_directories[i], "/output/image_analysis_summary_en.csv"), name_repair = "universal")
    df_metadata_dummy  <- readr::read_csv(file = paste0(input_directories[i], "/output/summary_metadata.csv"), name_repair = "universal")
    
    df_data <- dplyr::bind_rows(df_data, df_data_dummy)
    df_metadata <- dplyr::bind_rows(df_metadata, df_metadata_dummy)
  }
  
}
rm(list = c("df_data_dummy", "df_metadata_dummy", "i"))

# Data cleaning and add additional meta information ########################
output_dir_data <- paste0(output_dir, "data/")
dir.create(output_dir, output_dir_data = FALSE)

# Delete all rows that do not contain "NA" in manual_quality_check column anymore
df_data <- df_data[is.na(df_data$manual_quality_check),]

# Add column with information about date of experiment
df_data$date <- NA
df_data$date <- suppressWarnings(
  as.numeric(gsub(pattern = "^([0-9]{6}).+\\.czi$",
                  replacement = "\\1",
                  x = df_data$fileName)))
df_data$date <- lubridate::ymd(df_data$date)

# Add column with information of magnification of the objective of the microscope
df_data$magnification <- NA
df_data$magnification <- df_metadata$objective_magnification[match(df_metadata$fileName, df_data$fileName)]

# Add column with information about whether an image was retaken
df_data$image_retake <- FALSE
df_data$image_retake <- grepl(pattern = ".+_.+_", x = df_data$fileName, ignore.case = TRUE)

# Add column with information about which well the cells were from
df_data$well_number <- NA
# without retake
df_data$well_number[!df_data$image_retake] <- suppressWarnings(
  as.numeric(gsub(pattern = ".+([0-9]+)_[0-9]+\\.czi$",
                  replacement = "\\1",
                  x = df_data$fileName[!df_data$image_retake])))
#with retake
df_data$well_number[df_data$image_retake] <- suppressWarnings(
  as.numeric(gsub(pattern = ".+([0-9]+)_[0-9]+_[0-9]\\.czi$",
                  replacement = "\\1",
                  x = df_data$fileName[df_data$image_retake])))


# Add column with information about which image number the cells were from
df_data$image_number <- NA
# without retake
df_data$image_number[!df_data$image_retake] <- suppressWarnings(
  as.numeric(gsub(pattern = ".+_([0-9]+)\\.czi$",
                  replacement = "\\1",
                  x = df_data$fileName[!df_data$image_retake])))
#with retake
df_data$image_number[df_data$image_retake] <- suppressWarnings(
  as.numeric(gsub(pattern = ".+_([0-9]+)_[0-9]\\.czi$",
                  replacement = "\\1",
                  x = df_data$fileName[df_data$image_retake])))

# Add column with information about which group an experiment belongs to
df_data$experiment_group <- NA
df_data$experiment_group[
  grepl(pattern = " Stm", x = df_data$fileName, ignore.case = TRUE)] <- "Stimulation"
df_data$experiment_group[
  grepl(pattern = " poskon", x = df_data$fileName, ignore.case = TRUE)] <- "PositiveControl"
df_data$experiment_group[
  grepl(pattern = " kon", x = df_data$fileName, ignore.case = TRUE)] <- "Control"

# Reorder columns
df_data <- df_data %>% 
  dplyr::relocate(c("date", "experiment_group", "well_number", "image_number", "image_retake", "magnification" ),
                  .after = "fileName")

# Save final data frames
write.csv(x = df_data,
          file = paste(output_dir_data,"image_analysis_summary_complete.csv", sep=""),
          row.names = FALSE)
write.csv2(x = df_data,
           file = paste(output_dir_data,"image_analysis_summary_complete_de.csv", sep=""),
           row.names = FALSE)

write.csv(x = df_metadata,
          file = paste(output_dir_data,"summary_metadata_complete.csv", sep=""),
          row.names = FALSE)
write.csv2(x = df_metadata,
           file = paste(output_dir_data,"summary_metadata_complete_de.csv", sep=""),
           row.names = FALSE)

# Calculations of additional values to be added to the tibble ##############

# Green: corrected total fluorescence inside the entire cell
df_data$green_corrected_total_fluorescence_cell <- NA
df_data$green_corrected_total_fluorescence_cell <-
  df_data$intensity_sum_green_foreground -
  (df_data$intensity_mean_green_background*df_data$number_of_pixels_foreground)

# Green: corrected total fluorescence per number of nuclei inside the entire cell
df_data$green_corrected_total_fluorescence_cell_per_no_of_nuclei <- NA
df_data$green_corrected_total_fluorescence_cell_per_no_of_nuclei <-
  df_data$green_corrected_total_fluorescence_cell /
  df_data$number_of_nuclei

# Green: corrected total fluorescence above the nucleus
df_data$green_corrected_total_fluorescence_nucleus <- NA
df_data$green_corrected_total_fluorescence_nucleus <-
  df_data$intensity_sum_green_nucleus_region -
  (df_data$intensity_mean_green_background*df_data$number_of_pixels_nucleus_region)

# Green: corrected total fluorescence per number of nuclei above the nucleus
df_data$green_corrected_total_fluorescence_nucleus_per_no_of_nuclei <- NA
df_data$green_corrected_total_fluorescence_nucleus_per_no_of_nuclei <-
  df_data$green_corrected_total_fluorescence_nucleus /
  df_data$number_of_nuclei

# Red: corrected total fluorescence inside the entire cell
df_data$red_corrected_total_fluorescence_cell <- NA
df_data$red_corrected_total_fluorescence_cell <-
  df_data$intensity_sum_red_foreground -
  (df_data$intensity_mean_red_background*df_data$number_of_pixels_foreground)

# Red: corrected total fluorescence per number of nuclei inside the entire cell
df_data$red_corrected_total_fluorescence_cell_per_no_of_nuclei <- NA
df_data$red_corrected_total_fluorescence_cell_per_no_of_nuclei <-
  df_data$red_corrected_total_fluorescence_cell /
  df_data$number_of_nuclei

# Red: corrected total fluorescence above the nucleus
df_data$red_corrected_total_fluorescence_nucleus <- NA
df_data$red_corrected_total_fluorescence_nucleus <-
  df_data$intensity_sum_red_nucleus_region -
  (df_data$intensity_mean_red_background*df_data$number_of_pixels_nucleus_region)

# Red: corrected total fluorescence per number of nuclei above the nucleus
df_data$red_corrected_total_fluorescence_nucleus_per_no_of_nuclei <- NA
df_data$red_corrected_total_fluorescence_nucleus_per_no_of_nuclei <-
  df_data$red_corrected_total_fluorescence_nucleus /
  df_data$number_of_nuclei


# Grouping results #########################################################

# ROX red inside entire cells
print("ROX red inside the entire cell per number of nuclei")
# Calculate mean and standard error of mean
df_red_cell <- df_data %>% 
  filter(experiment_group != "PositiveControl") %>% 
  group_by(magnification, experiment_group) %>%
  summarise(mean_intensity = mean(ROX_red_corrected_total_fluorescence_cell_per_no_of_nuclei),
            uncertainty_intensity = qt(p = c(0.975),
                                       df = length(ROX_red_corrected_total_fluorescence_cell_per_no_of_nuclei)-1) *
              sd(ROX_red_corrected_total_fluorescence_cell_per_no_of_nuclei)/
              length(ROX_red_corrected_total_fluorescence_cell_per_no_of_nuclei),
            number_of_measurements= length(ROX_red_corrected_total_fluorescence_cell_per_no_of_nuclei))

df_red_cell_mean <- select(df_red_cell, select=-c(uncertainty_intensity)) %>% tidyr::spread(key = experiment_group, value = mean_intensity)
names(df_red_cell_mean)[names(df_red_cell_mean) == "Control"] <- "mean_control"
names(df_red_cell_mean)[names(df_red_cell_mean) == "Stimulation"] <- "mean_stim"
# names(df_red_cell_mean)[names(df_red_cell_mean) == "PositiveControl"] <- "mean_positiveControl"

df_red_cell_uncertainty <- select(df_red_cell, select=-c(mean_intensity)) %>% tidyr::spread(key = experiment_group, value = uncertainty_intensity)
names(df_red_cell_uncertainty)[names(df_red_cell_uncertainty) == "Control"] <- "uncertainty_control"
names(df_red_cell_uncertainty)[names(df_red_cell_uncertainty) == "Stimulation"] <- "uncertainty_stim"
# names(df_red_cell_uncertainty)[names(df_red_cell_uncertainty) == "PositiveControl"] <- "uncertainty_positiveControl"

df_red_cell_result <- merge(df_red_cell_mean, df_red_cell_uncertainty, by=intersect(names(df_red_cell_mean), names(df_red_cell_uncertainty)))

current_row <- !is.na(df_red_cell_result$mean_control)
df_red_cell_result$stim_over_control[current_row] <- df_red_cell_result$mean_stim[current_row] / df_red_cell_result$mean_control[current_row]
df_red_cell_result$uncertainty_stim_over_control[current_row] <- df_red_cell_result$uncertainty_stim[current_row] / df_red_cell_result$mean_control[current_row] +
  (df_red_cell_result$mean_stim[current_row] * df_red_cell_result$uncertainty_control[current_row])/
  (df_red_cell_result$mean_control[current_row] * df_red_cell_result$mean_control[current_row])

# current_row2 <- is.na(df_red_cell_result$mean_control)
# df_red_cell_result$positiveControl_over_control[current_row2] <- df_red_cell_result$mean_positiveControl[current_row2] / df_red_cell_result$mean_control[current_row]
# df_red_cell_result$uncertainty_positiveControl_over_control[current_row2] <- df_red_cell_result$uncertainty_positiveControl[current_row2] / df_red_cell_result$mean_control[current_row] +
#   (df_red_cell_result$mean_positiveControl[current_row2] * df_red_cell_result$uncertainty_control[current_row])/
#   (df_red_cell_result$mean_control[current_row] * df_red_cell_result$mean_control[current_row])

# ROX green inside nucleus

print("ROX green above/in nuclei per number of nuclei")
# Calculate mean and standard error of mean
df_green_nuc <- df_data %>% 
  filter(experiment_group != "PositiveControl") %>% 
  group_by(magnification, experiment_group) %>%
  summarise(mean_intensity = mean(ROX_green_corrected_total_fluorescence_nucleus_per_no_of_nuclei),
            uncertainty_intensity = qt(p = c(0.975),
                                       df = length(ROX_green_corrected_total_fluorescence_nucleus_per_no_of_nuclei)-1) *
              sd(ROX_green_corrected_total_fluorescence_nucleus_per_no_of_nuclei)/
              length(ROX_green_corrected_total_fluorescence_nucleus_per_no_of_nuclei),
            number_of_measurements= length(ROX_green_corrected_total_fluorescence_nucleus_per_no_of_nuclei))

df_green_nuc_mean <- select(df_green_nuc, select=-c(uncertainty_intensity)) %>% tidyr::spread(key = experiment_group, value = mean_intensity)
names(df_green_nuc_mean)[names(df_green_nuc_mean) == "Control"] <- "mean_control"
names(df_green_nuc_mean)[names(df_green_nuc_mean) == "Stimulation"] <- "mean_stim"

df_green_nuc_uncertainty <- select(df_green_nuc, select=-c(mean_intensity)) %>% tidyr::spread(key = experiment_group, value = uncertainty_intensity)
names(df_green_nuc_uncertainty)[names(df_green_nuc_uncertainty) == "Control"] <- "uncertainty_control"
names(df_green_nuc_uncertainty)[names(df_green_nuc_uncertainty) == "Stimulation"] <- "uncertainty_stim"

df_green_nuc_result <- merge(df_green_nuc_mean, df_green_nuc_uncertainty, by=intersect(names(df_green_nuc_mean), names(df_green_nuc_uncertainty)))

df_green_nuc_result$stim_over_control <- df_green_nuc_result$mean_stim / df_green_nuc_result$mean_control
df_green_nuc_result$uncertainty_stim_over_control <- df_green_nuc_result$uncertainty_stim / df_green_nuc_result$mean_control +
  (df_green_nuc_result$mean_stim*df_green_nuc_result$uncertainty_control)/(df_green_nuc_result$mean_control*df_green_nuc_result$mean_control)

# ROX green inside entire cells
print("ROX green inside the entire cell per number of nuclei")
# Calculate mean and standard error of mean
df_green_cell <- df_data %>% 
  filter(experiment_group != "PositiveControl") %>% 
  group_by(magnification, experiment_group) %>%
  summarise(mean_intensity = mean(ROX_green_corrected_total_fluorescence_cell_per_no_of_nuclei),
            uncertainty_intensity = qt(p = c(0.975),
                                       df = length(ROX_green_corrected_total_fluorescence_cell_per_no_of_nuclei)-1) *
              sd(ROX_green_corrected_total_fluorescence_cell_per_no_of_nuclei)/
              length(ROX_green_corrected_total_fluorescence_cell_per_no_of_nuclei),
            number_of_measurements= length(ROX_green_corrected_total_fluorescence_cell_per_no_of_nuclei))

df_green_cell_mean <- select(df_green_cell, select=-c(uncertainty_intensity)) %>% tidyr::spread(key = experiment_group, value = mean_intensity)
names(df_green_cell_mean)[names(df_green_cell_mean) == "Control"] <- "mean_control"
names(df_green_cell_mean)[names(df_green_cell_mean) == "Stimulation"] <- "mean_stim"

df_green_cell_uncertainty <- select(df_green_cell, select=-c(mean_intensity)) %>% tidyr::spread(key = experiment_group, value = uncertainty_intensity)
names(df_green_cell_uncertainty)[names(df_green_cell_uncertainty) == "Control"] <- "uncertainty_control"
names(df_green_cell_uncertainty)[names(df_green_cell_uncertainty) == "Stimulation"] <- "uncertainty_stim"

df_green_cell_result <- merge(df_green_cell_mean, df_green_cell_uncertainty, by=intersect(names(df_green_cell_mean), names(df_green_cell_uncertainty)))
df_green_cell_result$stim_over_control <- df_green_cell_result$mean_stim / df_green_cell_result$mean_control
df_green_cell_result$uncertainty_stim_over_control <- df_green_cell_result$uncertainty_stim / df_green_cell_result$mean_control +
  (df_green_cell_result$mean_stim*df_green_cell_result$uncertainty_control)/(df_green_cell_result$mean_control*df_green_cell_result$mean_control)

# Beta-Actin
# print("Beta-Actin inside/above nuclei per number of nuclei")
# df_green_nuc <- group_by(df_data, magnification, experiment_group) %>%
#   summarise(mean_intensity = mean(intensity_sum_green_nucleus_region_per_no_of_nuclei))
# print(df_green_nuc)
#
# print("Beta-Actin outside of nuclei per number of nuclei")
# df_green_cyt <- group_by(df_data, magnification, experiment_group) %>%
#   summarise(mean_intensity = mean(intensity_sum_green_outside_of_nucleus_region_per_no_of_nuclei))
# print(df_green_cyt)
#
# print("Beta-Actin in entire image per number of nuclei")
# df_green <- group_by(df_data, magnification, experiment_group) %>%
#   summarise(mean_intensity = mean(intensity_sum_green_per_no_of_nuclei))
# print(df_green)


# Plotting #################################################################
output_dir_plots <- paste0(output_dir, "plots/")
dir.create(output_dir_plots, showWarnings = FALSE)

# Filter experiment group = NA
df_data <- df_data %>% 
  dplyr::filter(!is.na(experiment_group))

# Number (no) of nuclei per image (blue)
df_data$experiment_group <- factor(df_data$experiment_group, levels = c("Control", "Stimulation", "PositiveControl"))

plot_no_nuclei <- ggplot(df_data, aes(y=number_of_nuclei,
                                      fill=experiment_group)) +
  geom_boxplot() +
  #facet_wrap(~end_time, scale="free") +
  theme_bw(base_size = 24) +
  theme(axis.title.y = element_text(size=18),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank()) +
  ylab("Number of nuclei per image") +
  xlab("") +
  ggtitle("Number of nuclei") +
  scale_fill_discrete(name = "Experiment Group",  labels = c("Control", "Stimulation", "Positive Control"))
#print(plot_no_nuclei)

ggsave(filename = paste(output_dir_plots, "number_of_nuclei.pdf", sep="/"),
       width = 297, height = 210, units = "mm")
ggsave(filename = paste(output_dir_plots, "number_of_nuclei.png", sep="/"),
       width = 297, height = 210, units = "mm")



# ROX red

# Total corrected fluorescence per number of nuclei
plot_ROX_red <- ggplot(df_data, aes(y=ROX_red_corrected_total_fluorescence_cell_per_no_of_nuclei,
                                    fill=experiment_group)) +
  geom_boxplot() +
  #facet_wrap(~end_time, scale="free") +
  theme_bw(base_size = 24) +
  theme(axis.title.y = element_text(size=18),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank()) +
  ylab("ROX red corrected total fluorescence in entire cell per number of nuclei") +
  xlab("") +
  ggtitle("ROX red in entire cell") +
  scale_fill_discrete(name = "Experiment Group",  labels = c("Control", "Stimulation", "Positive Control"))
#print(plot_ROX_red)

ggsave(filename = paste(output_dir_plots, "ROXred_cell_per_number_of_nuclei.pdf", sep="/"),
       width = 297, height = 210, units = "mm")
ggsave(filename = paste(output_dir_plots, "ROXred_cell_per_number_of_nuclei.png", sep="/"),
       width = 297, height = 210, units = "mm")


# Ratio of ROX red in cell per number of nuclei for Stim/Control
plot_ROX_red_ratio <- ggplot(df_red_cell_result, aes(x=0, y=stim_over_control)) +
  geom_point(size = 6) +
  geom_errorbar(aes(ymin=stim_over_control-uncertainty_stim_over_control,
                    ymax=stim_over_control+uncertainty_stim_over_control),
                size = 1.5, width=.2) +
  ylim(0,2) +
  xlim(-2,2) +
  theme_bw(base_size = 24) +
  theme(axis.title.y=element_text(size=18),
        #axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank()) +
  ggtitle("ROX red ratio in entire cell") +
  geom_hline(yintercept=1.0, linetype="dashed", size=2) +
  ylab("Ratio of corrected total fluorescence of ROX red\nin entire cell per number of nuclei (Stim/Control)") +
  xlab("")

#print(plot_ROX_red_ratio)

ggsave(filename = paste(output_dir_plots, "ROXred_cell_per_number_of_nuclei_ratio.pdf", sep="/"),
       width = 297, height = 210, units = "mm")
ggsave(filename = paste(output_dir_plots, "ROXred_cell_per_number_of_nuclei_ratio.png", sep="/"),
       width = 297, height = 210, units = "mm")





# ROX green

# Total corrected fluorescence above nuclei per number of nuclei
plot_ROX_green_nuc <- ggplot(df_data, aes(y=ROX_green_corrected_total_fluorescence_nucleus_per_no_of_nuclei,
                                    fill=experiment_group)) +
  geom_boxplot() +
  #facet_wrap(~end_time, scale="free") +
  theme_bw(base_size = 24) +
  theme(axis.title.y = element_text(size=18),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank()) +
  ylab("ROX green corrected total fluorescence\nabove nuclei per number of nuclei") +
  xlab("") +
  ggtitle("ROX green above nucleus") +
  scale_fill_discrete(name = "Experiment Group",  labels = c("Control", "Stimulation", "Positive Control"))
#print(plot_ROX_green_nuc)

ggsave(filename = paste(output_dir_plots, "ROXgreen_nuc_per_number_of_nuclei.pdf", sep="/"),
       width = 297, height = 210, units = "mm")
ggsave(filename = paste(output_dir_plots, "ROXgreen_nuc_per_number_of_nuclei.png", sep="/"),
       width = 297, height = 210, units = "mm")


# Total corrected fluorescence in entire cell per number of nuclei
plot_ROX_green_cell <- ggplot(df_data, aes(y=ROX_green_corrected_total_fluorescence_cell_per_no_of_nuclei,
                                          fill=experiment_group)) +
  geom_boxplot() +
  #facet_wrap(~end_time, scale="free") +
  theme_bw(base_size = 24) +
  theme(axis.title.y = element_text(size=18),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank()) +
  ylab("ROX green corrected total fluorescence\nin entire cell per number of nuclei") +
  xlab("") +
  ggtitle("ROX green in entire cell") +
  scale_fill_discrete(name = "Experiment Group",  labels = c("Control", "Stimulation", "Positive Control"))
#print(plot_ROX_green_cell)

ggsave(filename = paste(output_dir_plots, "ROXgreen_cell_per_number_of_nuclei.pdf", sep="/"),
       width = 297, height = 210, units = "mm")
ggsave(filename = paste(output_dir_plots, "ROXgreen_cell_per_number_of_nuclei.png", sep="/"),
       width = 297, height = 210, units = "mm")


# Total corrected fluorescence above nuclei per number of nuclei
plot_ROX_green_nuc <- ggplot(df_data, aes(y=ROX_green_corrected_total_fluorescence_nucleus_per_no_of_nuclei,
                                           fill=experiment_group)) +
  geom_boxplot() +
  #facet_wrap(~end_time, scale="free") +
  theme_bw(base_size = 24) +
  theme(axis.title.y = element_text(size=18),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank()) +
  ylab("ROX green corrected total fluorescence\nabove nuclei per number of nuclei") +
  xlab("") +
  ggtitle("ROX green above nuclei") +
  scale_fill_discrete(name = "Experiment Group",  labels = c("Control", "Stimulation", "Positive Control"))
#print(plot_ROX_green_nuc)

ggsave(filename = paste(output_dir_plots, "ROXgreen_nuc_per_number_of_nuclei.pdf", sep="/"),
       width = 297, height = 210, units = "mm")
ggsave(filename = paste(output_dir_plots, "ROXgreen_nuc_per_number_of_nuclei.png", sep="/"),
       width = 297, height = 210, units = "mm")




# ROX green inside/above nuclei in entire image per number of nuclei

plot_green_nuc <- ggplot(df_green_nuc_result, aes(x=0, y=stim_over_control)) +
  geom_point(size = 6) +
  geom_errorbar(aes(ymin=stim_over_control-uncertainty_stim_over_control,
                    ymax=stim_over_control+uncertainty_stim_over_control),
                size = 1.5, width=.2) +
  ylim(0,2) +
  xlim(-2,2) +
  theme_bw(base_size = 24) +
  theme(axis.title.y=element_text(size=18),
        #axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank()) +
  ggtitle("ROX green ratio above nucleus") +
  geom_hline(yintercept=1.0, linetype="dashed", size=2) +
  ylab("Ratio of corrected total fluorescence of ROX green\nabove nuclei per number of nuclei (Stim/Control)") +
  xlab("")

#print(plot_green_nuc)

ggsave(filename = paste(output_dir_plots, "ROXgreen_nuc_per_number_of_nuclei_ratio.pdf", sep="/"),
       width = 297, height = 210, units = "mm")
ggsave(filename = paste(output_dir_plots, "ROXgreen_nuc_per_number_of_nuclei_ratio.png", sep="/"),
       width = 297, height = 210, units = "mm")


# ROX green in cell in entire image per number of nuclei
plot_green_cell <- ggplot(df_green_cell_result, aes(x=0, y=stim_over_control)) +
  geom_point(size = 6) +
  geom_errorbar(aes(ymin=stim_over_control-uncertainty_stim_over_control,
                    ymax=stim_over_control+uncertainty_stim_over_control),
                size = 1.5, width=.2) +
  ylim(0,2) +
  xlim(-2,2) +
  theme_bw(base_size = 24) +
  theme(axis.title.y=element_text(size=18),
        #axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank()) +
  ggtitle("ROX green ratio in entire cell") +
  geom_hline(yintercept=1.0, linetype="dashed", size=2) +
  ylab("Ratio of corrected total fluorescence of ROX green\nin entire cell per number of nuclei (Stim/Control)") +
  xlab("")

#print(plot_green_cell)

ggsave(filename = paste(output_dir_plots, "ROXgreen_cell_per_number_of_nuclei_ratio.pdf", sep="/"),
       width = 297, height = 210, units = "mm")
ggsave(filename = paste(output_dir_plots, "ROXgreen_cell_per_number_of_nuclei_ratio.png", sep="/"),
       width = 297, height = 210, units = "mm")



# plot_red_nuc <- ggplot(df_data, aes(x=end_time,
#                                     y=intensity_sum_red_nucleus_region_per_no_of_nuclei,
#                                     fill=experiment_group)) +
#   geom_boxplot() +
#   #facet_wrap(~end_time, scale="free") +
#   theme_bw(base_size = 24) +
#   theme(axis.title.y=element_text(size=18),
#         #axis.text.x = element_blank(),
#         axis.ticks.x = element_blank()) +
#   ylab("Sum of red fluorescence intensity above/\nin nucleus per number of nuclei per image") +
#   xlab("") +
#   scale_fill_discrete(name = "Experiment Group",  labels = c("Control", "Stimulation"))
# print(plot_red_nuc)
#
# ggsave(filename = paste(output_dir, "betacat_nuc.pdf", sep="/"),
#        width = 297, height = 210, units = "mm")
# ggsave(filename = paste(output_dir, "betacat_nuc.png", sep="/"),
#        width = 297, height = 210, units = "mm")

#
# # Beta-Catenin inside/above nuclei in entire image per number of nuclei
# plot_red_cyt <- ggplot(df_data, aes(x=end_time,
#                                     y=intensity_sum_red_outside_of_nucleus_region_per_no_of_nuclei,
#                                     fill=experiment_group)) +
#   geom_boxplot() +
#   #facet_wrap(~end_time, scale="free") +
#   theme_bw(base_size = 24) +
#   theme(axis.title.y=element_text(size=18),
#         #axis.text.x = element_blank(),
#         axis.ticks.x = element_blank()) +
#   ylab("Sum of red fluorescence intensity outside of\nnuclei regions per number of nuclei per image") +
#   xlab("") +
#   scale_fill_discrete(name = "Experiment Group",  labels = c("Control", "Stimulation"))
# print(plot_red_cyt)
#
# ggsave(filename = paste(output_dir, "betacat_cyt.pdf", sep="/"),
#        width = 297, height = 210, units = "mm")
# ggsave(filename = paste(output_dir, "betacat_cyt.png", sep="/"),
#        width = 297, height = 210, units = "mm")
#
# # Beta-Actin (green)
#
# # Beta-Actin (green) outside of nuclei regions in entire image per number of nuclei
# plot_green_cyt <- ggplot(df_data, aes(x=end_time,
#                                       y=intensity_sum_green_outside_of_nucleus_region_per_no_of_nuclei ,
#                                       fill=experiment_group)) +
#   geom_boxplot() +
#   #facet_wrap(~magnification, scale="free") +
#   theme_bw(base_size = 24) +
#   theme(axis.title.y=element_text(size=18),
#         #axis.text.x = element_blank(),
#         axis.ticks.x = element_blank()) +
#   ylab("Sum of green fluorescence intensity outside of\nnuclei regions per number of nuclei per image") +
#   xlab("") +
#   scale_fill_discrete(name = "Experiment Group",  labels = c("Control", "Stimulation"))
# print(plot_green_cyt)
#
# ggsave(filename = paste(output_dir, "bactin_cyt.pdf", sep="/"),
#        width = 297, height = 210, units = "mm")
# ggsave(filename = paste(output_dir, "bactin_cyt.png", sep="/"),
#        width = 297, height = 210, units = "mm")
#
#
# # Beta-Actin (green) insinde of/above nuclei regions in entire image per number of nuclei
# plot_green_nuc <- ggplot(df_data, aes(x=end_time,
#                                       y=intensity_sum_green_nucleus_region_per_no_of_nuclei ,
#                                       fill=experiment_group)) +
#   geom_boxplot() +
#   #facet_wrap(~magnification, scale="free") +
#   theme_bw(base_size = 24) +
#   theme(axis.title.y=element_text(size=18),
#         #axis.text.x = element_blank(),
#         axis.ticks.x = element_blank()) +
#   ylab("Sum of green fluorescence intensity above/\nin nucleus per number of nuclei per image") +
#   xlab("") +
#   scale_fill_discrete(name = "Experiment Group",  labels = c("Control", "Stimulation"))
# print(plot_green_nuc)
#
# ggsave(filename = paste(output_dir, "bactin_nuc.pdf", sep="/"),
#        width = 297, height = 210, units = "mm")
# ggsave(filename = paste(output_dir, "bactin_nuc.png", sep="/"),
#        width = 297, height = 210, units = "mm")
#
#
# # Beta-Actin (green) in entire image per number of nuclei
# plot_green <- ggplot(df_data, aes(x=end_time,
#                                       y=intensity_sum_green_per_no_of_nuclei ,
#                                       fill=experiment_group)) +
#   geom_boxplot() +
#   #facet_wrap(~magnification, scale="free") +
#   theme_bw(base_size = 24) +
#   theme(axis.title.y=element_text(size=18),
#         #axis.text.x = element_blank(),
#         axis.ticks.x = element_blank()) +
#   ylab("Sum of green fluorescence intensity\nin entire image per number of nuclei per image") +
#   xlab("") +
#   scale_fill_discrete(name = "Experiment Group",  labels = c("Control", "Stimulation"))
# print(plot_green)
#
# ggsave(filename = paste(output_dir, "bactin_nuc.pdf", sep="/"),
#        width = 297, height = 210, units = "mm")
# ggsave(filename = paste(output_dir, "bactin_nuc.png", sep="/"),
#        width = 297, height = 210, units = "mm")

