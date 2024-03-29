# Script for plotting waveforms and measurements from ++++++++++++++++++++++
# oscilloscope recordings (zipped files)              ++++++++++++++++++++++
# Author: Kai Budde
# Created: 2022/05/18
# Last changed: 2023/02/21

# Delete everything in the environment
rm(list = ls())
# close all open plots in RStudio
graphics.off()

# Load packages ############################################################warnings

# Set groundhog day for reproducibility (see https://groundhogr.com)
groundhog.day <- "2022-03-01"

if(!any(grepl(pattern = "groundhog", x = installed.packages(), ignore.case = TRUE))){
  install.packages("groundhog")
}

# Load packages
library(groundhog)
pkgs <- c("tidyverse", "yaml", "scattermore")
groundhog.library(pkgs, groundhog.day)


# Load functions from package ##############################################
# TODO: This should be changed to loading the package once its done.

source("C:/Users/Kai/Documents/git/gitHub/oscilloscopeR/R/convertMeasurementListToTibble.R")
source("C:/Users/Kai/Documents/git/gitHub/oscilloscopeR/R/getMeasurements.R")
source("C:/Users/Kai/Documents/git/gitHub/oscilloscopeR/R/getWaveforms.R")
source("C:/Users/Kai/Documents/git/gitHub/oscilloscopeR/R/plotMeasurements.R")
source("C:/Users/Kai/Documents/git/gitHub/oscilloscopeR/R/plotWaveforms.R")

# Data input ###############################################################
# Directory (containing directories with recording (yaml files))

# input_directory <- "data/stimulation"
# output_dir <- "plots/stimulation"
input_directory <- "E:/PhD/Daten/ShortStim_ZellBio/Voltage3"
stim_data_name_pattern <- "mA" #ShortStimCellBio
output_dir <- file.path(input_directory, "plots/")

ES_channel_name <- "ElectricalStim"
function_generator_channel_name <- "FunGen"
resistor_channel_name <- "Resistor"

# Go through every recording ###############################################
dir.create(output_dir, showWarnings = FALSE)

recording_files <- list.files(path = input_directory, recursive = FALSE, full.names = TRUE)
# Keep those with stim_data_name_pattern in name

# Check if files are folders or zips
zip_files <- any(grepl(pattern = stim_data_name_pattern, x = recording_files, ignore.case = TRUE) &
  grepl(pattern = "zip", x = recording_files, ignore.case = TRUE))

if(zip_files){
  recording_files <- recording_files[
    grepl(pattern = stim_data_name_pattern, x = recording_files, ignore.case = TRUE) &
      grepl(pattern = "zip", x = recording_files, ignore.case = TRUE)]
  
}else{
  recording_files <- recording_files[
    grepl(pattern = stim_data_name_pattern, x = recording_files, ignore.case = TRUE)]
  
}

# Create empty tibble for statistics on p2p data
df_p2p_statistics <- tibble(
  "directoryName" = basename(recording_files),
  "date" = as.Date(NA),
  "ES_p2p_voltage_mean" = as.numeric(NA),
  "ES_p2p_voltage_median" = as.numeric(NA),
  "ES_p2p_voltage_sd" = as.numeric(NA),
  "ES_avg_voltage_mean" = as.numeric(NA),
  "ES_avg_voltage_median" = as.numeric(NA),
  "ES_avg_voltage_sd" = as.numeric(NA),
  "shunt_p2p_voltage_mean" = as.numeric(NA),
  "shunt_p2p_voltage_median" = as.numeric(NA),
  "shunt_p2p_voltage_sd" = as.numeric(NA),
  "shunt_avg_voltage_mean" = as.numeric(NA),
  "shunt_avg_voltage_median" = as.numeric(NA),
  "shunt_avg_voltage_sd" = as.numeric(NA),
  "frequency" = as.numeric(NA))

for(current_file in recording_files){
  
  # current_output_dir <- file.path(output_dir, basename(current_file))
  # current_output_dir <- gsub(pattern = "\\.zip", replacement = "", x = current_output_dir, ignore.case = TRUE)
  # dir.create(current_output_dir, showWarnings = FALSE)
  # plot_title <- gsub(pattern = "(.+)_.+", replacement = "\\1", x = basename(current_file))
  plot_title <- gsub(pattern = "(.+)\\..*$", replacement = "\\1", x = basename(current_file))
  
  plot_title1 <- gsub(pattern = "([0-9]+)_.+", replacement = "\\1", x = plot_title)
  plot_title1 <- lubridate::as_date(plot_title1)
  plot_title2 <- gsub(pattern = "[0-9]+_(.+)$", replacement = "\\1", x = plot_title)
  plot_title2 <- gsub(pattern = "_", replacement = " ", x = plot_title2)
  
  plot_title <- paste(plot_title1, plot_title2, sep=" ")
    
  # Plot measurement data ##################################################
  if(zip_files){
    df_data <- getMeasurements(input_file = current_file)
  }else{
    df_data <- getMeasurements(input_directory = current_file)
  }
  
  
  # Get information about the channels (which one is the function generator and which one the stimulator)
  # Use the channel (out of two) with least "noise" as the function generator channel
  channels <- unique(df_data$Channel)
  
  if(length(channels) == 1){
    df_data$Channel <- "StimPulse"
    # channel_stimulation_pulse <- paste("CHAN", channels, sep="")
  }else if(length(channels) == 2){
    
    sd_ch1 <- df_data %>% filter(Channel == channels[1] & VPP <1e3 & VPP > 1e-3) %>% summarise(sd(VPP))
    sd_ch2 <- df_data %>% filter(Channel == channels[2] & VPP <1e3 & VPP > 1e-3) %>% summarise(sd(VPP))
    
    # lowest sd -> function generator
    # highest sd -> stimulator output
    
    ES_channel <- which(c(sd_ch1, sd_ch2) == max(sd_ch1, sd_ch2))
    FunGen_channel <- which(c(sd_ch1, sd_ch2) == min(sd_ch1, sd_ch2))
    
    df_data$Channel[df_data$Channel == channels[ES_channel]]  <- ES_channel_name
    df_data$Channel[df_data$Channel == channels[FunGen_channel]]  <- function_generator_channel_name
    
  }else if(length(channels) == 3){
    if(length(df_data) > 10){
      sd_ch1 <- df_data %>% filter(Channel == channels[1] & VPP <1e3 & VPP > 1e-3) %>% summarise(sd(VPP))
      sd_ch2 <- df_data %>% filter(Channel == channels[2] & VPP <1e3 & VPP > 1e-3) %>% summarise(sd(VPP))
      sd_ch3 <- df_data %>% filter(Channel == channels[3] & VPP <1e3 & VPP > 1e-3) %>% summarise(sd(VPP))
      
      # lowest sd -> function generator
      # middle sd -> measurement at resistor
      # highest sd -> stimulator output
      
      ES_channel <- which(c(sd_ch1, sd_ch2, sd_ch3) == max(sd_ch1, sd_ch2, sd_ch3))
      FunGen_channel <- which(c(sd_ch1, sd_ch2, sd_ch3) == min(sd_ch1, sd_ch2, sd_ch3))
      resistor_channel <- which(c(sd_ch1, sd_ch2, sd_ch3) != max(sd_ch1, sd_ch2, sd_ch3) & c(sd_ch1, sd_ch2, sd_ch3) != min(sd_ch1, sd_ch2, sd_ch3))
      
      ES_channel <- channels[ES_channel]
      FunGen_channel <- channels[FunGen_channel]
      resistor_channel <- channels[resistor_channel]
      
    }else{
      # Not enough points
      
      # lowest Vpp -> function generator
      # middle Vpp -> measurement at resistor
      # highest Vpp -> stimulator output
      
      FunGen_channel   <- df_data$Channel[which(df_data$VPP == min(df_data$VPP[df_data$VPP < 1e3]))]
      ES_channel       <- df_data$Channel[which(df_data$VPP == max(df_data$VPP[df_data$VPP < 1e3]))]
      resistor_channel <- channels[!channels %in% c(FunGen_channel, ES_channel[df_data$VPP < 1e3])][1]
      
    }
    
    df_data$Channel[df_data$Channel == ES_channel]  <- ES_channel_name
    df_data$Channel[df_data$Channel == FunGen_channel]  <- function_generator_channel_name
    df_data$Channel[df_data$Channel == resistor_channel]  <- resistor_channel_name
    
    
  }else{
    print("More than 2 channels used. Please check.")
  }
  
  max_voltage <- max(df_data$VPP[df_data$Channel == ES_channel_name &
                                   df_data$VPP < 1e3])
  if(is.na(max_voltage)){
    print("Max voltage is na. Please check")
  }else if(max_voltage < 10){
    voltage_limit <- 10
  }else if(max_voltage < 20){
    voltage_limit <- 20
  }else if(max_voltage < 30){
    voltage_limit <- 30
  }else if(max_voltage < 40){
    voltage_limit <- 40
  }else{
    voltage_limit <- 1.01*max_voltage
  }
  
  plotMeasurements(input_data = df_data,
                   output_dir = output_dir,
                   vpp_min = 0,
                   vpp_max = voltage_limit,
                   vavg_min = -2,
                   vavg_max = 2,
                   plot_title = plot_title)
  
  # Save results in separate tibble
  df_p2p_statistics$date[df_p2p_statistics$directoryName == basename(current_file)] <- unique(as.Date(df_data$date_time))[1]
  
  df_p2p_statistics$ES_p2p_voltage_mean[df_p2p_statistics$directoryName == basename(current_file)] <- 
    mean(df_data$VPP[df_data$Channel == ES_channel_name])
  df_p2p_statistics$ES_p2p_voltage_median[df_p2p_statistics$directoryName == basename(current_file)] <- 
    median(df_data$VPP[df_data$Channel == ES_channel_name])
  df_p2p_statistics$ES_p2p_voltage_sd[df_p2p_statistics$directoryName == basename(current_file)] <- 
    sd(df_data$VPP[df_data$Channel == ES_channel_name])
  df_p2p_statistics$ES_avg_voltage_mean[df_p2p_statistics$directoryName == basename(current_file)] <- 
    mean(df_data$VAVG[df_data$Channel == ES_channel_name])
  df_p2p_statistics$ES_avg_voltage_median[df_p2p_statistics$directoryName == basename(current_file)] <- 
    median(df_data$VAVG[df_data$Channel == ES_channel_name])
  df_p2p_statistics$ES_avg_voltage_sd[df_p2p_statistics$directoryName == basename(current_file)] <- 
    sd(df_data$VAVG[df_data$Channel == ES_channel_name])
  
  if(length(channels) == 3){
      df_p2p_statistics$shunt_p2p_voltage_mean[df_p2p_statistics$directoryName == basename(current_file)] <- 
      mean(df_data$VPP[df_data$Channel == resistor_channel_name])
    df_p2p_statistics$shunt_p2p_voltage_median[df_p2p_statistics$directoryName == basename(current_file)] <- 
      median(df_data$VPP[df_data$Channel == resistor_channel_name])
    df_p2p_statistics$shunt_p2p_voltage_sd[df_p2p_statistics$directoryName == basename(current_file)] <- 
      sd(df_data$VPP[df_data$Channel == resistor_channel_name])
    df_p2p_statistics$shunt_avg_voltage_mean[df_p2p_statistics$directoryName == basename(current_file)] <- 
      mean(df_data$VAVG[df_data$Channel == resistor_channel_name])
    df_p2p_statistics$shunt_avg_voltage_median[df_p2p_statistics$directoryName == basename(current_file)] <- 
      median(df_data$VAVG[df_data$Channel == resistor_channel_name])
    df_p2p_statistics$shunt_avg_voltage_sd[df_p2p_statistics$directoryName == basename(current_file)] <- 
      sd(df_data$VAVG[df_data$Channel == resistor_channel_name])
  }
  
  if("FREQ" %in% names(df_data)){
    df_p2p_statistics$frequency[df_p2p_statistics$directoryName == basename(current_file)] <- 
      as.numeric(names(sort(-table(df_data$FREQ[df_data$Channel == function_generator_channel_name])))[1])
  }
  
  
  # Plot waveform data #####################################################
  if(zip_files){
    df_data <- getWaveforms(input_file = current_file)
  }else{
    df_data <- getWaveforms(input_directory = current_file)
  }
  
  
  if(exists("ES_channel")){
    df_data$Channel[df_data$Channel == ES_channel]  <- ES_channel_name
  }
  if(exists("FunGen_channel")){
    df_data$Channel[df_data$Channel == FunGen_channel]  <- function_generator_channel_name
  }
  if(exists("resistor_channel")){
    df_data$Channel[df_data$Channel == resistor_channel]  <- resistor_channel_name
  }
  
  max_voltage <- max(df_data$U[df_data$Channel == ES_channel_name &
                                 df_data$U < 1e3])
  if(is.na(max_voltage)){
    print("Max voltage is na. Please check")
  }else if(max_voltage < 10){
    voltage_limit <- 10
  }else if(max_voltage < 15){
    voltage_limit <- 15
  }else if(max_voltage < 20){
    voltage_limit <- 20
  }else if(max_voltage < 30){
    voltage_limit <- 30
  }else{
    voltage_limit <- 1.01*max_voltage
  }

  plotWaveforms(input_data = df_data,
                output_dir = output_dir,
                show_time_in_us = FALSE,
                channel_function_generator = function_generator_channel_name,
                channel_stimulation_pulse = ES_channel_name,
                plot_waveforms = "one",
                voltage_limits_of_plot = voltage_limit,
                filter_stim_off = TRUE,
                epsilon_for_filtering = 1,
                plot_title = plot_title)
  
}

write_csv(x = df_p2p_statistics, file = file.path(output_dir, "p2p_statistics.csv"))
