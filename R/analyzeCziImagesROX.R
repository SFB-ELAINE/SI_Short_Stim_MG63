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
dir.create(output_dir_data, showWarnings = FALSE)

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






# Filter data frame and then calculate addition values and plots them ######

# Filter experiment group = NA
df_data <- df_data %>% 
  dplyr::filter(!is.na(experiment_group))

# Add "x" to magnification column
df_data$magnification <- paste0(df_data$magnification, "x")


filter_well_number <- c(1, 2)

for(i in 0:length(filter_well_number)){
  
  if(i == 0){
    # Do not filter
    df_data_dummy <- df_data
    append_text <- ""
  }else{
    # Filter
    df_data_dummy <- df_data %>% 
      dplyr::filter(well_number == i)
    append_text <- paste0("(Well ", i, ")")
  }
  
  
  # Calculations of additional values to be added to the tibble ##############
  
  # Green: corrected total fluorescence inside the entire cell
  df_data_dummy$green_corrected_total_fluorescence_cell <- NA
  df_data_dummy$green_corrected_total_fluorescence_cell <-
    df_data_dummy$intensity_sum_green_foreground -
    (df_data_dummy$intensity_mean_green_background*df_data_dummy$number_of_pixels_foreground)
  
  # Green: corrected total fluorescence per number of nuclei inside the entire cell
  df_data_dummy$green_corrected_total_fluorescence_cell_per_no_of_nuclei <- NA
  df_data_dummy$green_corrected_total_fluorescence_cell_per_no_of_nuclei <-
    df_data_dummy$green_corrected_total_fluorescence_cell /
    df_data_dummy$number_of_nuclei
  
  # Green: corrected total fluorescence above the nucleus
  df_data_dummy$green_corrected_total_fluorescence_nucleus <- NA
  df_data_dummy$green_corrected_total_fluorescence_nucleus <-
    df_data_dummy$intensity_sum_green_nucleus_region -
    (df_data_dummy$intensity_mean_green_background*df_data_dummy$number_of_pixels_nucleus_region)
  
  # Green: corrected total fluorescence per number of nuclei above the nucleus
  df_data_dummy$green_corrected_total_fluorescence_nucleus_per_no_of_nuclei <- NA
  df_data_dummy$green_corrected_total_fluorescence_nucleus_per_no_of_nuclei <-
    df_data_dummy$green_corrected_total_fluorescence_nucleus /
    df_data_dummy$number_of_nuclei
  
  # Red: corrected total fluorescence inside the entire cell
  df_data_dummy$red_corrected_total_fluorescence_cell <- NA
  df_data_dummy$red_corrected_total_fluorescence_cell <-
    df_data_dummy$intensity_sum_red_foreground -
    (df_data_dummy$intensity_mean_red_background*df_data_dummy$number_of_pixels_foreground)
  
  # Red: corrected total fluorescence per number of nuclei inside the entire cell
  df_data_dummy$red_corrected_total_fluorescence_cell_per_no_of_nuclei <- NA
  df_data_dummy$red_corrected_total_fluorescence_cell_per_no_of_nuclei <-
    df_data_dummy$red_corrected_total_fluorescence_cell /
    df_data_dummy$number_of_nuclei
  
  # Red: corrected total fluorescence above the nucleus
  df_data_dummy$red_corrected_total_fluorescence_nucleus <- NA
  df_data_dummy$red_corrected_total_fluorescence_nucleus <-
    df_data_dummy$intensity_sum_red_nucleus_region -
    (df_data_dummy$intensity_mean_red_background*df_data_dummy$number_of_pixels_nucleus_region)
  
  # Red: corrected total fluorescence per number of nuclei above the nucleus
  df_data_dummy$red_corrected_total_fluorescence_nucleus_per_no_of_nuclei <- NA
  df_data_dummy$red_corrected_total_fluorescence_nucleus_per_no_of_nuclei <-
    df_data_dummy$red_corrected_total_fluorescence_nucleus /
    df_data_dummy$number_of_nuclei
  
  
  # Grouping results #########################################################
  
  # RED ----------------------------------------------------------------------
  # ROX red inside entire cells
  print("ROX red inside the entire cell per number of nuclei")
  # Calculate mean and standard error of mean
  df_red_cell <- df_data_dummy %>% 
    filter(experiment_group != "PositiveControl") %>% 
    group_by(magnification, experiment_group) %>%
    summarise(mean_intensity = mean(red_corrected_total_fluorescence_cell_per_no_of_nuclei),
              uncertainty_intensity = qt(p = c(0.975),
                                         df = length(red_corrected_total_fluorescence_cell_per_no_of_nuclei)-1) *
                sd(red_corrected_total_fluorescence_cell_per_no_of_nuclei)/
                length(red_corrected_total_fluorescence_cell_per_no_of_nuclei),
              number_of_measurements= length(red_corrected_total_fluorescence_cell_per_no_of_nuclei))
  
  df_red_cell_mean <- select(df_red_cell, select=-c(uncertainty_intensity, number_of_measurements)) %>% tidyr::spread(key = experiment_group, value = mean_intensity)
  names(df_red_cell_mean)[names(df_red_cell_mean) == "Control"] <- "mean_control"
  names(df_red_cell_mean)[names(df_red_cell_mean) == "Stimulation"] <- "mean_stim"
  # names(df_red_cell_mean)[names(df_red_cell_mean) == "PositiveControl"] <- "mean_positiveControl"
  
  df_red_cell_uncertainty <- select(df_red_cell, select=-c(mean_intensity, number_of_measurements)) %>% tidyr::spread(key = experiment_group, value = uncertainty_intensity)
  names(df_red_cell_uncertainty)[names(df_red_cell_uncertainty) == "Control"] <- "uncertainty_control"
  names(df_red_cell_uncertainty)[names(df_red_cell_uncertainty) == "Stimulation"] <- "uncertainty_stim"
  # names(df_red_cell_uncertainty)[names(df_red_cell_uncertainty) == "PositiveControl"] <- "uncertainty_positiveControl"
  
  df_red_cell_result <- merge(df_red_cell_mean, df_red_cell_uncertainty, by=intersect(names(df_red_cell_mean), names(df_red_cell_uncertainty)))
  
  current_row <- !is.na(df_red_cell_result$mean_control)
  df_red_cell_result$stim_over_control[current_row] <- df_red_cell_result$mean_stim[current_row] / df_red_cell_result$mean_control[current_row]
  df_red_cell_result$uncertainty_stim_over_control[current_row] <- df_red_cell_result$uncertainty_stim[current_row] / df_red_cell_result$mean_control[current_row] +
    (df_red_cell_result$mean_stim[current_row] * df_red_cell_result$uncertainty_control[current_row])/
    (df_red_cell_result$mean_control[current_row] * df_red_cell_result$mean_control[current_row])
  
  # ROX red inside nucleus
  print("ROX red above/in nuclei per number of nuclei")
  # Calculate mean and standard error of mean
  df_red_nuc <- df_data_dummy %>% 
    filter(experiment_group != "PositiveControl") %>% 
    group_by(magnification, experiment_group) %>%
    summarise(mean_intensity = mean(red_corrected_total_fluorescence_nucleus_per_no_of_nuclei),
              uncertainty_intensity = qt(p = c(0.975),
                                         df = length(red_corrected_total_fluorescence_nucleus_per_no_of_nuclei)-1) *
                sd(red_corrected_total_fluorescence_nucleus_per_no_of_nuclei)/
                length(red_corrected_total_fluorescence_nucleus_per_no_of_nuclei),
              number_of_measurements= length(red_corrected_total_fluorescence_nucleus_per_no_of_nuclei))
  
  df_red_nuc_mean <- select(df_red_nuc, select=-c(uncertainty_intensity, number_of_measurements)) %>% tidyr::spread(key = experiment_group, value = mean_intensity)
  names(df_red_nuc_mean)[names(df_red_nuc_mean) == "Control"] <- "mean_control"
  names(df_red_nuc_mean)[names(df_red_nuc_mean) == "Stimulation"] <- "mean_stim"
  
  df_red_nuc_uncertainty <- select(df_red_nuc, select=-c(mean_intensity, number_of_measurements)) %>% tidyr::spread(key = experiment_group, value = uncertainty_intensity)
  names(df_red_nuc_uncertainty)[names(df_red_nuc_uncertainty) == "Control"] <- "uncertainty_control"
  names(df_red_nuc_uncertainty)[names(df_red_nuc_uncertainty) == "Stimulation"] <- "uncertainty_stim"
  
  df_red_nuc_result <- merge(df_red_nuc_mean, df_red_nuc_uncertainty, by=intersect(names(df_red_nuc_mean), names(df_red_nuc_uncertainty)))
  
  df_red_nuc_result$stim_over_control <- df_red_nuc_result$mean_stim / df_red_nuc_result$mean_control
  df_red_nuc_result$uncertainty_stim_over_control <- df_red_nuc_result$uncertainty_stim / df_red_nuc_result$mean_control +
    (df_red_nuc_result$mean_stim*df_red_nuc_result$uncertainty_control)/(df_red_nuc_result$mean_control*df_red_nuc_result$mean_control)
  
  
  # GREEN --------------------------------------------------------------------
  
  # ROX green inside entire cells
  print("ROX green inside the entire cell per number of nuclei")
  # Calculate mean and standard error of mean
  df_green_cell <- df_data_dummy %>% 
    filter(experiment_group != "PositiveControl") %>% 
    group_by(magnification, experiment_group) %>%
    summarise(mean_intensity = mean(green_corrected_total_fluorescence_cell_per_no_of_nuclei),
              uncertainty_intensity = qt(p = c(0.975),
                                         df = length(green_corrected_total_fluorescence_cell_per_no_of_nuclei)-1) *
                sd(green_corrected_total_fluorescence_cell_per_no_of_nuclei)/
                length(green_corrected_total_fluorescence_cell_per_no_of_nuclei),
              number_of_measurements= length(green_corrected_total_fluorescence_cell_per_no_of_nuclei))
  
  df_green_cell_mean <- select(df_green_cell, select=-c(uncertainty_intensity, number_of_measurements)) %>% tidyr::spread(key = experiment_group, value = mean_intensity)
  names(df_green_cell_mean)[names(df_green_cell_mean) == "Control"] <- "mean_control"
  names(df_green_cell_mean)[names(df_green_cell_mean) == "Stimulation"] <- "mean_stim"
  
  df_green_cell_uncertainty <- select(df_green_cell, select=-c(mean_intensity, number_of_measurements)) %>% tidyr::spread(key = experiment_group, value = uncertainty_intensity)
  names(df_green_cell_uncertainty)[names(df_green_cell_uncertainty) == "Control"] <- "uncertainty_control"
  names(df_green_cell_uncertainty)[names(df_green_cell_uncertainty) == "Stimulation"] <- "uncertainty_stim"
  
  df_green_cell_result <- merge(df_green_cell_mean, df_green_cell_uncertainty, by=intersect(names(df_green_cell_mean), names(df_green_cell_uncertainty)))
  df_green_cell_result$stim_over_control <- df_green_cell_result$mean_stim / df_green_cell_result$mean_control
  df_green_cell_result$uncertainty_stim_over_control <- df_green_cell_result$uncertainty_stim / df_green_cell_result$mean_control +
    (df_green_cell_result$mean_stim*df_green_cell_result$uncertainty_control)/(df_green_cell_result$mean_control*df_green_cell_result$mean_control)
  
  # ROX green inside nucleus
  
  print("ROX green above/in nuclei per number of nuclei")
  # Calculate mean and standard error of mean
  df_green_nuc <- df_data_dummy %>% 
    filter(experiment_group != "PositiveControl") %>% 
    group_by(magnification, experiment_group) %>%
    summarise(mean_intensity = mean(green_corrected_total_fluorescence_nucleus_per_no_of_nuclei),
              uncertainty_intensity = qt(p = c(0.975),
                                         df = length(green_corrected_total_fluorescence_nucleus_per_no_of_nuclei)-1) *
                sd(green_corrected_total_fluorescence_nucleus_per_no_of_nuclei)/
                length(green_corrected_total_fluorescence_nucleus_per_no_of_nuclei),
              number_of_measurements= length(green_corrected_total_fluorescence_nucleus_per_no_of_nuclei))
  
  df_green_nuc_mean <- select(df_green_nuc, select=-c(uncertainty_intensity, number_of_measurements)) %>% tidyr::spread(key = experiment_group, value = mean_intensity)
  names(df_green_nuc_mean)[names(df_green_nuc_mean) == "Control"] <- "mean_control"
  names(df_green_nuc_mean)[names(df_green_nuc_mean) == "Stimulation"] <- "mean_stim"
  
  df_green_nuc_uncertainty <- select(df_green_nuc, select=-c(mean_intensity, number_of_measurements)) %>% tidyr::spread(key = experiment_group, value = uncertainty_intensity)
  names(df_green_nuc_uncertainty)[names(df_green_nuc_uncertainty) == "Control"] <- "uncertainty_control"
  names(df_green_nuc_uncertainty)[names(df_green_nuc_uncertainty) == "Stimulation"] <- "uncertainty_stim"
  
  df_green_nuc_result <- merge(df_green_nuc_mean, df_green_nuc_uncertainty, by=intersect(names(df_green_nuc_mean), names(df_green_nuc_uncertainty)))
  
  df_green_nuc_result$stim_over_control <- df_green_nuc_result$mean_stim / df_green_nuc_result$mean_control
  df_green_nuc_result$uncertainty_stim_over_control <- df_green_nuc_result$uncertainty_stim / df_green_nuc_result$mean_control +
    (df_green_nuc_result$mean_stim*df_green_nuc_result$uncertainty_control)/(df_green_nuc_result$mean_control*df_green_nuc_result$mean_control)
  
  
  # Plotting ###############################################################
  output_dir_plots <- paste0(output_dir, "plots/")
  dir.create(output_dir_plots, showWarnings = FALSE)
  
  plot_title_nuc <- "ROX"
  plot_title_green <- "CellROX Green"
  plot_title_red <- "CellROX Red"
  
  # Number (no) of nuclei per image (blue) ---------------------------------
  
  df_data_dummy$experiment_group <- factor(df_data_dummy$experiment_group, levels = c("Control", "Stimulation", "PositiveControl"))
  
  plot_no_nuclei <- ggplot(df_data_dummy, aes(x=experiment_group, y=number_of_nuclei)) +
    geom_violin() +
    # stat_boxplot(geom ='errorbar', width = 0.3) +
    geom_boxplot(width=0.05) +
    facet_wrap(.~ magnification, scale="free") +
    # geom_jitter(color="black", size=0.5, alpha=0.9) +
    #ylim(0,20) +
    theme_bw(base_size = 24) +
    theme(axis.title.y = element_text(size=18)) +
    xlab("") +
    ylab("Number of nuclei per image") +
    ggtitle(paste0("Number of nuclei (", plot_title_nuc, ") ", append_text)) +
    scale_fill_discrete(name = "Experiment Group",  labels = c("Control", "Stimulation", "Positive Control"))
  
  # print(plot_no_nuclei)
  
  ggsave(filename = paste(output_dir_plots, plot_title_nuc, "_number_of_nuclei_violin", append_text, ".pdf", sep=""),
         width = 297, height = 210, units = "mm")
  ggsave(filename = paste(output_dir_plots, plot_title_nuc, "_number_of_nuclei_violin", append_text, ".png", sep=""),
         width = 297, height = 210, units = "mm")
  
  plot_no_nuclei <- ggplot(df_data_dummy, aes(y=number_of_nuclei,
                                              fill=experiment_group)) +
    geom_boxplot() +
    facet_wrap(.~ magnification, scale="free") +
    theme_bw(base_size = 24) +
    theme(axis.title.y = element_text(size=18),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank()) +
    ylab("Number of nuclei per image") +
    xlab("") +
    ggtitle(paste0("Number of nuclei (", plot_title_nuc, ") ", append_text)) +
    scale_fill_discrete(name = "Experiment Group",  labels = c("Control", "Stimulation", "Positive Control"))
  
  #print(plot_no_nuclei)
  
  ggsave(filename = paste(output_dir_plots, plot_title_nuc, "_number_of_nuclei", append_text, ".pdf", sep=""),
         width = 297, height = 210, units = "mm")
  ggsave(filename = paste(output_dir_plots, plot_title_nuc, "_number_of_nuclei", append_text, ".png", sep=""),
         width = 297, height = 210, units = "mm")
  
  
  # Total corrected red fluorescence ---------------------------------------
  
  plot_red_nuc <- ggplot(df_data_dummy, aes(y=red_corrected_total_fluorescence_nucleus_per_no_of_nuclei,
                                            x=experiment_group)) +
    geom_violin() +
    geom_boxplot(width=0.05) +
    facet_wrap(.~ magnification, scale="free") +
    scale_y_continuous(trans='log10') +
    theme_bw(base_size = 24) +
    theme(axis.title.y = element_text(size=18)) +
    ylab("Total corrected red fluorescence\nabove nuclei per number of nuclei") +
    xlab("") +
    ggtitle(paste0(plot_title_red, " above nucleus ", append_text)) +
    scale_fill_discrete(name = "Experiment Group",  labels = c("Control", "Stimulation"))
  #print(plot_red_nuc)
  
  ggsave(filename = paste(output_dir_plots, plot_title_nuc, "_red_nuc_per_number_of_nuclei_violin", append_text, ".pdf", sep=""),
         width = 297, height = 210, units = "mm")
  ggsave(filename = paste(output_dir_plots, plot_title_nuc, "_red_nuc_per_number_of_nuclei_violin", append_text, ".png", sep=""),
         width = 297, height = 210, units = "mm")
  
  
  plot_red_nuc <- ggplot(df_data_dummy, aes(y=red_corrected_total_fluorescence_nucleus_per_no_of_nuclei,
                                            fill=experiment_group)) +
    geom_boxplot() +
    #facet_wrap(~end_time, scale="free") +
    facet_wrap(.~ magnification, scale="free") +
    scale_y_continuous(trans='log10') +
    theme_bw(base_size = 24) +
    theme(axis.title.y = element_text(size=18),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank()) +
    ylab("Total corrected red fluorescence\nabove nuclei per number of nuclei") +
    xlab("") +
    ggtitle(paste0(plot_title_red, " above nucleus ", append_text)) +
    scale_fill_discrete(name = "Experiment Group",  labels = c("Control", "Stimulation"))
  #print(plot_red_nuc)
  
  ggsave(filename = paste(output_dir_plots, plot_title_nuc, "_red_nuc_per_number_of_nuclei", append_text, ".pdf", sep=""),
         width = 297, height = 210, units = "mm")
  ggsave(filename = paste(output_dir_plots, plot_title_nuc, "_red_nuc_per_number_of_nuclei", append_text, ".png", sep=""),
         width = 297, height = 210, units = "mm")
  
  
  
  plot_red_cell <- ggplot(df_data_dummy, aes(y=red_corrected_total_fluorescence_cell_per_no_of_nuclei,
                                             x=experiment_group)) +
    geom_violin() +
    facet_wrap(.~ magnification, scale="free") +
    scale_y_continuous(trans='log10') +
    geom_boxplot(width=0.05) +
    theme_bw(base_size = 24) +
    theme(axis.title.y = element_text(size=18)) +
    ylab("Total corrected red fluorescence\nin entire cell per number of nuclei") +
    xlab("") +
    ggtitle(paste0(plot_title_red, " in entire cell ", append_text)) +
    scale_fill_discrete(name = "Experiment Group",  labels = c("Control", "Stimulation"))
  #print(plot_red_cell)
  
  ggsave(filename = paste(output_dir_plots, plot_title_nuc, "_red_cell_per_number_of_nuclei_violin", append_text, ".pdf", sep=""),
         width = 297, height = 210, units = "mm")
  ggsave(filename = paste(output_dir_plots, plot_title_nuc, "_red_cell_per_number_of_nuclei_violin", append_text, ".png", sep=""),
         width = 297, height = 210, units = "mm")
  
  plot_red_cell <- ggplot(df_data_dummy, aes(y=red_corrected_total_fluorescence_cell_per_no_of_nuclei,
                                             fill=experiment_group)) +
    geom_boxplot() +
    facet_wrap(.~ magnification, scale="free") +
    scale_y_continuous(trans='log10') +
    #facet_wrap(~end_time, scale="free") +
    theme_bw(base_size = 24) +
    theme(axis.title.y = element_text(size=18),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank()) +
    ylab("Total corrected red fluorescence\nin entire cell per number of nuclei") +
    xlab("") +
    ggtitle(paste0(plot_title_red, " in entire cell ", append_text)) +
    scale_fill_discrete(name = "Experiment Group",  labels = c("Control", "Stimulation"))
  #print(plot_red_cell)
  
  ggsave(filename = paste(output_dir_plots, plot_title_nuc, "_red_cell_per_number_of_nuclei", append_text, ".pdf", sep=""),
         width = 297, height = 210, units = "mm")
  ggsave(filename = paste(output_dir_plots, plot_title_nuc, "_red_cell_per_number_of_nuclei", append_text, ".png", sep=""),
         width = 297, height = 210, units = "mm")
  
  
  # Ratio of red fluorescence in cell per number of nuclei for Stim/Control
  plot_red_ratio_cell <- ggplot(df_red_cell_result, aes(x=0, y=stim_over_control)) +
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
    ggtitle("ROX red ratio in entire cell ", append_text) +
    geom_hline(yintercept=1.0, linetype="dashed", size=2) +
    ylab("Ratio of total corrected red fluorescence\nin entire cell per number of nuclei (Stim/Control)") +
    xlab("")
  
  #print(plot_red_ratio_cell)
  
  ggsave(filename = paste(output_dir_plots, plot_title_nuc, "_red_cell_per_number_of_nuclei_ratio", append_text, ".pdf", sep=""),
         width = 297, height = 210, units = "mm")
  ggsave(filename = paste(output_dir_plots, plot_title_nuc, "_red_cell_per_number_of_nuclei_ratio", append_text, ".png", sep=""),
         width = 297, height = 210, units = "mm")
  
  
  # Ratio of red fluorescence in cell per number of nuclei for Stim/Control
  plot_red_ratio_nuc <- ggplot(df_red_nuc_result, aes(x=0, y=stim_over_control)) +
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
    ggtitle("ROX red ratio above nuclei") +
    geom_hline(yintercept=1.0, linetype="dashed", size=2) +
    ylab("Ratio of total corrected red fluorescence\nabove nuclei per number of nuclei (Stim/Control)") +
    xlab("")
  
  #print(plot_red_ratio_nuc)
  
  ggsave(filename = paste(output_dir_plots, plot_title_nuc, "_red_nuc_per_number_of_nuclei_ratio", append_text, ".pdf", sep=""),
         width = 297, height = 210, units = "mm")
  ggsave(filename = paste(output_dir_plots, plot_title_nuc, "_red_nuc_per_number_of_nuclei_ratio", append_text, ".png", sep=""),
         width = 297, height = 210, units = "mm")
  
  
  
  # Total corrected green fluorescence -------------------------------------
  
  
  plot_green_nuc <- ggplot(df_data_dummy, aes(y=green_corrected_total_fluorescence_nucleus_per_no_of_nuclei,
                                              x=experiment_group)) +
    geom_violin() +
    geom_boxplot(width=0.05) +
    facet_wrap(.~ magnification, scale="free") +
    scale_y_continuous(trans='log10') +
    theme_bw(base_size = 24) +
    theme(axis.title.y = element_text(size=18)) +
    ylab("Total corrected green fluorescence\nabove nuclei per number of nuclei") +
    xlab("") +
    ggtitle(paste0(plot_title_green, " above nucleus ", append_text)) +
    scale_fill_discrete(name = "Experiment Group",  labels = c("Control", "Stimulation"))
  #print(plot_green_nuc)
  
  ggsave(filename = paste(output_dir_plots, plot_title_nuc, "_green_nuc_per_number_of_nuclei_violin", append_text, ".pdf", sep=""),
         width = 297, height = 210, units = "mm")
  ggsave(filename = paste(output_dir_plots, plot_title_nuc, "_green_nuc_per_number_of_nuclei_violin", append_text, ".png", sep=""),
         width = 297, height = 210, units = "mm")
  
  
  plot_green_nuc <- ggplot(df_data_dummy, aes(y=green_corrected_total_fluorescence_nucleus_per_no_of_nuclei,
                                              fill=experiment_group)) +
    geom_boxplot() +
    #facet_wrap(~end_time, scale="free") +
    facet_wrap(.~ magnification, scale="free") +
    scale_y_continuous(trans='log10') +
    theme_bw(base_size = 24) +
    theme(axis.title.y = element_text(size=18),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank()) +
    ylab("Total corrected green fluorescence\nabove nuclei per number of nuclei") +
    xlab("") +
    ggtitle(paste0(plot_title_green, " above nucleus ", append_text)) +
    scale_fill_discrete(name = "Experiment Group",  labels = c("Control", "Stimulation"))
  #print(plot_green_nuc)
  
  ggsave(filename = paste(output_dir_plots, plot_title_nuc, "_green_nuc_per_number_of_nuclei", append_text, ".pdf", sep=""),
         width = 297, height = 210, units = "mm")
  ggsave(filename = paste(output_dir_plots, plot_title_nuc, "_green_nuc_per_number_of_nuclei", append_text, ".png", sep=""),
         width = 297, height = 210, units = "mm")
  
  
  plot_green_cell <- ggplot(df_data_dummy, aes(y=green_corrected_total_fluorescence_cell_per_no_of_nuclei,
                                               x=experiment_group)) +
    geom_violin() +
    facet_wrap(.~ magnification, scale="free") +
    scale_y_continuous(trans='log10') +
    geom_boxplot(width=0.05) +
    theme_bw(base_size = 24) +
    theme(axis.title.y = element_text(size=18)) +
    ylab("Total corrected green fluorescence\nin entire cell per number of nuclei") +
    xlab("") +
    ggtitle(paste0(plot_title_green, " in entire cell ", append_text)) +
    scale_fill_discrete(name = "Experiment Group",  labels = c("Control", "Stimulation"))
  #print(plot_green_cell)
  
  ggsave(filename = paste(output_dir_plots, plot_title_nuc, "_green_cell_per_number_of_nuclei_violin", append_text, ".pdf", sep=""),
         width = 297, height = 210, units = "mm")
  ggsave(filename = paste(output_dir_plots, plot_title_nuc, "_green_cell_per_number_of_nuclei_violin", append_text, ".png", sep=""),
         width = 297, height = 210, units = "mm")
  
  plot_green_cell <- ggplot(df_data_dummy, aes(y=green_corrected_total_fluorescence_cell_per_no_of_nuclei,
                                               fill=experiment_group)) +
    geom_boxplot() +
    facet_wrap(.~ magnification, scale="free") +
    scale_y_continuous(trans='log10') +
    #facet_wrap(~end_time, scale="free") +
    theme_bw(base_size = 24) +
    theme(axis.title.y = element_text(size=18),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank()) +
    ylab("Total corrected green fluorescence\nin entire cell per number of nuclei") +
    xlab("") +
    ggtitle(paste0(plot_title_green, " in entire cell ", append_text)) +
    scale_fill_discrete(name = "Experiment Group",  labels = c("Control", "Stimulation"))
  #print(plot_green_cell)
  
  ggsave(filename = paste(output_dir_plots, plot_title_nuc, "_green_cell_per_number_of_nuclei", append_text, ".pdf", sep=""),
         width = 297, height = 210, units = "mm")
  ggsave(filename = paste(output_dir_plots, plot_title_nuc, "_green_cell_per_number_of_nuclei", append_text, ".png", sep=""),
         width = 297, height = 210, units = "mm")
  
  
  # Ratio of green fluorescence in cell per number of nuclei for Stim/Control
  plot_green_ratio_cell <- ggplot(df_green_cell_result, aes(x=0, y=stim_over_control)) +
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
    ggtitle("ROX green ratio in entire cell ", append_text) +
    geom_hline(yintercept=1.0, linetype="dashed", size=2) +
    ylab("Ratio of total corrected green fluorescence\nin entire cell per number of nuclei (Stim/Control)") +
    xlab("")
  
  #print(plot_green_ratio_cell)
  
  ggsave(filename = paste(output_dir_plots, plot_title_nuc, "_green_cell_per_number_of_nuclei_ratio", append_text, ".pdf", sep=""),
         width = 297, height = 210, units = "mm")
  ggsave(filename = paste(output_dir_plots, plot_title_nuc, "_green_cell_per_number_of_nuclei_ratio", append_text, ".png", sep=""),
         width = 297, height = 210, units = "mm")
  
  
  # Ratio of green fluorescence in cell per number of nuclei for Stim/Control
  plot_green_ratio_nuc <- ggplot(df_green_nuc_result, aes(x=0, y=stim_over_control)) +
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
    ggtitle("ROX green ratio above nuclei") +
    geom_hline(yintercept=1.0, linetype="dashed", size=2) +
    ylab("Ratio of total corrected green fluorescence\nabove nuclei per number of nuclei (Stim/Control)") +
    xlab("")
  
  #print(plot_green_ratio_nuc)
  
  ggsave(filename = paste(output_dir_plots, plot_title_nuc, "_green_nuc_per_number_of_nuclei_ratio", append_text, ".pdf", sep=""),
         width = 297, height = 210, units = "mm")
  ggsave(filename = paste(output_dir_plots, plot_title_nuc, "_green_nuc_per_number_of_nuclei_ratio", append_text, ".png", sep=""),
         width = 297, height = 210, units = "mm")
  
  
}





