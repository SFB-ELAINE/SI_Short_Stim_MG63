# Script for analyzing czi images                                  +++++++++
# Author: Kai Budde
# Created: 2022/06/30
# Last changed: 2022/07/05

# Delete everything in the environment
rm(list = ls())
# close all open plots in RStudio
graphics.off()

# Please adapt the following parameters ####################################
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

# Directory of the images with CellROX measurements
input_directories_LipidPeroxidation <- c(
  "E:/PhD/Daten/ShortStim_ZellBio/LipidPeroxidation/01.04.2022 24h/",
  "E:/PhD/Daten/ShortStim_ZellBio/LipidPeroxidation/05.04.2022 2h/")

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
  input_directories <- c(input_directories_LipidPeroxidation)
  
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
  input_directories <- c(input_directories_LipidPeroxidation)
  
  for(i in 1:length(input_directories)){
    
    current_input_dir <- input_directories[i]
    # Data frame with certain parameter values
    file_names <- list.files(path = current_input_dir)
    file_names_czi <- file_names[grepl("czi", file_names)]
    file_names_czi <- paste(current_input_dir, file_names_czi, sep="/")
    number_of_czi_files <- length(file_names_czi)
    
    
    for(j in 1:number_of_czi_files){
      if(j==1){
        df_metadata <- readCzi::readCziMetadata(input_file = file_names_czi[j], save_metadata = FALSE)
      }else{
        df_dummy <- readCzi::readCziMetadata(input_file = file_names_czi[j], save_metadata = FALSE)
        df_metadata <- rbind(df_metadata, df_dummy)
        rm(df_dummy)
      }
    }
    rm(j)
    
    output_dir <- paste0(current_input_dir, "output/")
    dir.create(output_dir, showWarnings = FALSE)
    
    write.csv(x = df_metadata,
              file = paste(output_dir,"summary_metadata.csv", sep=""),
              row.names = FALSE)
    write.csv2(x = df_metadata,
               file = paste(output_dir,"summary_metadata_de.csv", sep=""),
               row.names = FALSE)
    
    
  }
  rm(i)
  
}

## Data anylsis ------------------------------------------------------------

for(i in 1:length(input_directories_LipidPeroxidation)){
  # Import data
  df_data <- readr::read_csv(file = paste0(input_directories_LipidPeroxidation[i], "output/image_analysis_summary_en.csv"), name_repair = "universal")
  df_metadata <- readr::read_csv(file = paste0(input_directories_LipidPeroxidation[i], "output/summary_metadata.csv"), name_repair = "universal")
  
  output_dir <- paste0(input_directories_LipidPeroxidation[i], "output/plots/")
  dir.create(output_dir, showWarnings = FALSE)
  
  # Delete all rows that do not contain "NA" in manual_quality_check column anymore
  df_data <- df_data[is.na(df_data$manual_quality_check),]
  
  # Add column with information of magnification of the objective of the microscope
  df_data$magnification <- NA
  df_data$magnification <- df_metadata$objective_magnification[match(df_metadata$fileName, df_data$fileName)]
  
  # Add column with information about which group an experiment belongs to
  df_data$experiment_number <- NA
  df_data$experiment_number <- suppressWarnings(
    as.numeric(gsub(pattern = ".+z1_([0-9]{0,2})_.+$",
                    replacement = "\\1",
                    x = df_data$fileName)))
  # df_data$experiment_number[is.na(df_data$experiment_number)] <- 1
  
  # Add column with information about which group an experiment belongs to
  df_data$experiment_group <- NA
  df_data$experiment_group[
    grepl(pattern = "^Stim", x = df_data$fileName, ignore.case = TRUE)] <- "Stimulation"
  df_data$experiment_group[
    grepl(pattern = "^PosKon", x = df_data$fileName, ignore.case = TRUE)] <- "PositiveControl"
  df_data$experiment_group[
    grepl(pattern = "^Kon", x = df_data$fileName, ignore.case = TRUE)] <- "Control"
  
  # Calculations of additional values to be added to the df ################
  
  # Red (unoxidized): corrected total fluorescence inside the entire cell
  df_data$unoxidized_red_corrected_total_fluorescence_cell <- NA
  df_data$unoxidized_red_corrected_total_fluorescence_cell <-
    df_data$intensity_sum_red_foreground -
    (df_data$intensity_mean_red_background*df_data$number_of_pixels_foreground)
  
  df_data$unoxidized_red_corrected_total_fluorescence_cell_per_no_of_nuclei <- NA
  df_data$unoxidized_red_corrected_total_fluorescence_cell_per_no_of_nuclei <-
    df_data$unoxidized_red_corrected_total_fluorescence_cell /
    df_data$number_of_nuclei
  
  
  # Green (oxidized): corrected total fluorescence inside the nucleus and inside
  # the entire cell
  df_data$oxidized_green_corrected_total_fluorescence_nucleus <- NA
  df_data$oxidized_green_corrected_total_fluorescence_nucleus <-
    df_data$intensity_sum_green_nucleus_region -
    (df_data$intensity_mean_green_background*df_data$number_of_pixels_nucleus_region)
  
  df_data$oxidized_green_corrected_total_fluorescence_nucleus_per_no_of_nuclei <- NA
  df_data$oxidized_green_corrected_total_fluorescence_nucleus_per_no_of_nuclei <-
    df_data$oxidized_green_corrected_total_fluorescence_nucleus /
    df_data$number_of_nuclei
  
  
  df_data$oxidized_green_corrected_total_fluorescence_cell <- NA
  df_data$oxidized_green_corrected_total_fluorescence_cell <-
    df_data$intensity_sum_green_foreground -
    (df_data$intensity_mean_green_background*df_data$number_of_pixels_foreground)
  
  df_data$oxidized_green_corrected_total_fluorescence_cell_per_no_of_nuclei <- NA
  df_data$oxidized_green_corrected_total_fluorescence_cell_per_no_of_nuclei <-
    df_data$oxidized_green_corrected_total_fluorescence_cell /
    df_data$number_of_nuclei
  
  # Grouping results #########################################################
  
  # Unoxidized red inside entire cells
  
  print("Unoxidized red inside the entire cell per number of nuclei")
  
  # Calculate mean and standard error of mean
  df_red_cell <- df_data %>% 
    filter(experiment_group != "PositiveControl") %>% 
    group_by(magnification, experiment_group) %>%
    summarise(mean_intensity = mean(unoxidized_red_corrected_total_fluorescence_cell_per_no_of_nuclei),
              uncertainty_intensity = qt(p = c(0.975),
                                         df = length(unoxidized_red_corrected_total_fluorescence_cell_per_no_of_nuclei)-1) *
                sd(unoxidized_red_corrected_total_fluorescence_cell_per_no_of_nuclei)/
                length(unoxidized_red_corrected_total_fluorescence_cell_per_no_of_nuclei),
              number_of_measurements= length(unoxidized_red_corrected_total_fluorescence_cell_per_no_of_nuclei))
  
  df_red_cell <- df_red_cell %>% select(-"number_of_measurements")
  
  df_red_cell_mean <- select(df_red_cell, select=-c(uncertainty_intensity)) %>% tidyr::spread(key = experiment_group, value = mean_intensity)
  names(df_red_cell_mean)[names(df_red_cell_mean) == "Control"] <- "mean_control"
  names(df_red_cell_mean)[names(df_red_cell_mean) == "Stimulation"] <- "mean_stim"
  # names(df_red_cell_mean)[names(df_red_cell_mean) == "PositiveControl"] <- "mean_positiveControl"
  
  df_red_cell_uncertainty <- select(df_red_cell, select=-c(mean_intensity)) %>% tidyr::spread(key = experiment_group, value = uncertainty_intensity)
  names(df_red_cell_uncertainty)[names(df_red_cell_uncertainty) == "Control"] <- "uncertainty_control"
  names(df_red_cell_uncertainty)[names(df_red_cell_uncertainty) == "Stimulation"] <- "uncertainty_stim"
  # names(df_red_cell_uncertainty)[names(df_red_cell_uncertainty) == "PositiveControl"] <- "uncertainty_positiveControl"
  
  df_red_cell_result <- merge(df_red_cell_mean, df_red_cell_uncertainty, by=intersect(names(df_red_cell_mean), names(df_red_cell_uncertainty)))
  
  df_red_cell_result$stim_over_control <- df_red_cell_result$mean_stim / df_red_cell_result$mean_control
  df_red_cell_result$uncertainty_stim_over_control <- df_red_cell_result$uncertainty_stim / df_red_cell_result$mean_control +
    (df_red_cell_result$mean_stim * df_red_cell_result$uncertainty_control)/
    (df_red_cell_result$mean_control * df_red_cell_result$mean_control)
  
  # current_row2 <- is.na(df_red_cell_result$mean_control)
  # df_red_cell_result$positiveControl_over_control[current_row2] <- df_red_cell_result$mean_positiveControl[current_row2] / df_red_cell_result$mean_control[current_row]
  # df_red_cell_result$uncertainty_positiveControl_over_control[current_row2] <- df_red_cell_result$uncertainty_positiveControl[current_row2] / df_red_cell_result$mean_control[current_row] +
  #   (df_red_cell_result$mean_positiveControl[current_row2] * df_red_cell_result$uncertainty_control[current_row])/
  #   (df_red_cell_result$mean_control[current_row] * df_red_cell_result$mean_control[current_row])
  
  # Oxidized green inside nucleus
  
  print("Oxidized green above/in nuclei per number of nuclei")
  # Calculate mean and standard error of mean
  df_green_nuc <- df_data %>% 
    filter(experiment_group != "PositiveControl") %>% 
    group_by(magnification, experiment_group) %>%
    summarise(mean_intensity = mean(oxidized_green_corrected_total_fluorescence_nucleus_per_no_of_nuclei),
              uncertainty_intensity = qt(p = c(0.975),
                                         df = length(oxidized_green_corrected_total_fluorescence_nucleus_per_no_of_nuclei)-1) *
                sd(oxidized_green_corrected_total_fluorescence_nucleus_per_no_of_nuclei)/
                length(oxidized_green_corrected_total_fluorescence_nucleus_per_no_of_nuclei),
              number_of_measurements= length(oxidized_green_corrected_total_fluorescence_nucleus_per_no_of_nuclei))
  
  df_green_nuc <- df_green_nuc %>% select(-"number_of_measurements")
  
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
  
  # Oxidized green inside entire cells
  print("Oxidized green inside the entire cell per number of nuclei")
  # Calculate mean and standard error of mean
  df_green_cell <- df_data %>% 
    filter(experiment_group != "PositiveControl") %>% 
    group_by(magnification, experiment_group) %>%
    summarise(mean_intensity = mean(oxidized_green_corrected_total_fluorescence_cell_per_no_of_nuclei),
              uncertainty_intensity = qt(p = c(0.975),
                                         df = length(oxidized_green_corrected_total_fluorescence_cell_per_no_of_nuclei)-1) *
                sd(oxidized_green_corrected_total_fluorescence_cell_per_no_of_nuclei)/
                length(oxidized_green_corrected_total_fluorescence_cell_per_no_of_nuclei),
              number_of_measurements= length(oxidized_green_corrected_total_fluorescence_cell_per_no_of_nuclei))
  
  df_green_cell <- df_green_cell %>% select(-"number_of_measurements")
  
  
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
  
  # Plotting #################################################################
  
  # Number (no) of nuclei per image (blue)
  df_data$experiment_group <- factor(df_data$experiment_group,
                                     levels = c("Control", "Stimulation", "PositiveControl"))
  
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
  
  ggsave(filename = paste(output_dir, "lipidPeroxidatoin_number_of_nuclei.pdf", sep="/"),
         width = 297, height = 210, units = "mm")
  ggsave(filename = paste(output_dir, "lipidPeroxidatoin_number_of_nuclei.png", sep="/"),
         width = 297, height = 210, units = "mm")
  
  
  
  # Unoxidized red
  
  # Total corrected fluorescence per number of nuclei
  plot_unoxidized_red <- ggplot(df_data, aes(y=unoxidized_red_corrected_total_fluorescence_cell_per_no_of_nuclei,
                                             fill=experiment_group)) +
    geom_boxplot() +
    #facet_wrap(~end_time, scale="free") +
    theme_bw(base_size = 24) +
    theme(axis.title.y = element_text(size=18),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank()) +
    ylab("Corrected total red fluorescence\nin entire cell per number of nuclei") +
    xlab("") +
    ggtitle("Unoxidized BODIPY(TM) 581/591 C11 in entire cell") +
    scale_fill_discrete(name = "Experiment Group",  labels = c("Control", "Stimulation", "Positive Control"))
  #print(plot_unoxidized_red)
  
  ggsave(filename = paste(output_dir, "UnoxidizedRed_cell_per_number_of_nuclei.pdf", sep="/"),
         width = 297, height = 210, units = "mm")
  ggsave(filename = paste(output_dir, "UnoxidizedRed_cell_per_number_of_nuclei.png", sep="/"),
         width = 297, height = 210, units = "mm")
  
  
  # Ratio of Unoxidized red in cell per number of nuclei for Stim/Control
  plot_unoxidized_red_ratio <- ggplot(df_red_cell_result, aes(x=0, y=stim_over_control)) +
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
    ggtitle("Unoxidized BODIPY(TM) 581/591 C11 ratio in entire cell") +
    geom_hline(yintercept=1.0, linetype="dashed", size=2) +
    ylab("Ratio of corrected total red fluorescence (unoxidized reagent)\nin entire cell per number of nuclei (Stim/Control)") +
    xlab("")
  
  #print(plot_unoxidized_red_ratio)
  
  ggsave(filename = paste(output_dir, "UnoxidizedRed_cell_per_number_of_nuclei_ratio.pdf", sep="/"),
         width = 297, height = 210, units = "mm")
  ggsave(filename = paste(output_dir, "UnoxidizedRed_cell_per_number_of_nuclei_ratio.png", sep="/"),
         width = 297, height = 210, units = "mm")
  
  
  
  
  
  # Oxidized green
  
  # Total corrected fluorescence above nuclei per number of nuclei
  plot_oxidized_green_nuc <- ggplot(df_data, aes(y=oxidized_green_corrected_total_fluorescence_nucleus_per_no_of_nuclei,
                                                 fill=experiment_group)) +
    geom_boxplot() +
    #facet_wrap(~end_time, scale="free") +
    theme_bw(base_size = 24) +
    theme(axis.title.y = element_text(size=18),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank()) +
    ylab("Corrected total green fluorescence\nabove nuclei per number of nuclei") +
    xlab("") +
    ggtitle("Oxidized BODIPY(TM) 581/591 C11 above nucleus") +
    scale_fill_discrete(name = "Experiment Group",  labels = c("Control", "Stimulation", "Positive Control"))
  #print(plot_oxidized_green_nuc)
  
  ggsave(filename = paste(output_dir, "OxidizedGreen_nuc_per_number_of_nuclei.pdf", sep="/"),
         width = 297, height = 210, units = "mm")
  ggsave(filename = paste(output_dir, "OxidizedGreen_nuc_per_number_of_nuclei.png", sep="/"),
         width = 297, height = 210, units = "mm")
  
  
  # Total corrected fluorescence in entire cell per number of nuclei
  plot_oxidized_green_cell <- ggplot(df_data, aes(y=oxidized_green_corrected_total_fluorescence_cell_per_no_of_nuclei,
                                                  fill=experiment_group)) +
    geom_boxplot() +
    #facet_wrap(~end_time, scale="free") +
    theme_bw(base_size = 24) +
    theme(axis.title.y = element_text(size=18),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank()) +
    ylab("Corrected total green fluorescence\nin entire cell per number of nuclei") +
    xlab("") +
    ggtitle("Oxidized BODIPY(TM) 581/591 C11 in entire cell") +
    scale_fill_discrete(name = "Experiment Group",  labels = c("Control", "Stimulation", "Positive Control"))
  #print(plot_oxidized_green_cell)
  
  ggsave(filename = paste(output_dir, "OxidizedGreen_cell_per_number_of_nuclei.pdf", sep="/"),
         width = 297, height = 210, units = "mm")
  ggsave(filename = paste(output_dir, "OxidizedGreen_cell_per_number_of_nuclei.png", sep="/"),
         width = 297, height = 210, units = "mm")
  
  
  # Total corrected fluorescence above nuclei per number of nuclei
  plot_oxidized_green_nuc <- ggplot(df_data, aes(y=oxidized_green_corrected_total_fluorescence_nucleus_per_no_of_nuclei,
                                                 fill=experiment_group)) +
    geom_boxplot() +
    #facet_wrap(~end_time, scale="free") +
    theme_bw(base_size = 24) +
    theme(axis.title.y = element_text(size=18),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank()) +
    ylab("Corrected total green fluorescence\nabove nuclei per number of nuclei") +
    xlab("") +
    ggtitle("Oxidized BODIPY(TM) 581/591 C11 above nuclei") +
    scale_fill_discrete(name = "Experiment Group",  labels = c("Control", "Stimulation", "Positive Control"))
  #print(plot_oxidized_green_nuc)
  
  ggsave(filename = paste(output_dir, "OxidizedGreen_nuc_per_number_of_nuclei.pdf", sep="/"),
         width = 297, height = 210, units = "mm")
  ggsave(filename = paste(output_dir, "OxidizedGreen_nuc_per_number_of_nuclei.png", sep="/"),
         width = 297, height = 210, units = "mm")
  
  
  
  
  # Oxidized green inside/above nuclei in entire image per number of nuclei
  
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
    ggtitle("Oxidized BODIPY(TM) 581/591 C11 ratio above nucleus") +
    geom_hline(yintercept=1.0, linetype="dashed", size=2) +
    ylab("Ratio of corrected total green fluorescence (oxidized reagent)\nabove nuclei per number of nuclei (Stim/Control)") +
    xlab("")
  
  #print(plot_green_nuc)
  
  ggsave(filename = paste(output_dir, "OxidizedGreen_nuc_per_number_of_nuclei_ratio.pdf", sep="/"),
         width = 297, height = 210, units = "mm")
  ggsave(filename = paste(output_dir, "OxidizedGreen_nuc_per_number_of_nuclei_ratio.png", sep="/"),
         width = 297, height = 210, units = "mm")
  
  
  # Oxidized green in cell in entire image per number of nuclei
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
    ggtitle("Oxidized BODIPY(TM) 581/591 C11 ratio in entire cell") +
    geom_hline(yintercept=1.0, linetype="dashed", size=2) +
    ylab("Ratio of corrected total green fluorescence of (oxidized reagent)\nin entire cell per number of nuclei (Stim/Control)") +
    xlab("")
  
  #print(plot_green_cell)
  
  ggsave(filename = paste(output_dir, "OxidizedGreen_cell_per_number_of_nuclei_ratio.pdf", sep="/"),
         width = 297, height = 210, units = "mm")
  ggsave(filename = paste(output_dir, "OxidizedGreen_cell_per_number_of_nuclei_ratio.png", sep="/"),
         width = 297, height = 210, units = "mm")
  
  
  
}
rm(i)


