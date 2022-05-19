# Script for plotting waveforms and measurements from ++++++++++++++++++++++
# oscilloscope recordings (zipped files)              ++++++++++++++++++++++
# Author: Kai Budde
# Created: 2022/05/18
# Last changed: 2022/05/19

# Delete everything in the environment
rm(list = ls())
# close all open plots in RStudio
graphics.off()

# Load packages ############################################################

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
input_directory <- "data/stimulation"
output_dir <- "plots/stimulation"

ES_channel_name <- "ElectricalStim"
function_generator_channel_name <- "FunGen"


# Go through every recording ###############################################
dir.create(output_dir, showWarnings = FALSE)

recording_files <- list.files(path = input_directory, recursive = FALSE, full.names = TRUE)
# Keep those with "ShortStimCellBio" in name
recording_files <- recording_files[
  grepl(pattern = "ShortStimCellBio", x = recording_files, ignore.case = TRUE) &
    grepl(pattern = "zip", x = recording_files, ignore.case = TRUE)]

for(current_zip_file in recording_files){
  
  current_output_dir <- file.path(output_dir, basename(current_zip_file))
  current_output_dir <- gsub(pattern = "\\.zip", replacement = "", x = current_output_dir, ignore.case = TRUE)
  dir.create(current_output_dir, showWarnings = FALSE)
  plot_title <- gsub(pattern = "(.+)_.+", replacement = "\\1", x = basename(current_zip_file))
  plot_title <- lubridate::as_date(plot_title)
  
  # Plot measurement data ####################################################
  df_data <- getMeasurements(input_file = current_zip_file)
  
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
    
    # if(sd_ch1 > sd_ch2){
    #   # ch1 was the stimulus
    #   ES_channel <- 1
    #   FunGen_channel <- 2
    #   df_data$Channel[df_data$Channel == channels[ES_channel]]  <- ES_channel_name
    #   df_data$Channel[df_data$Channel == channels[FunGen_channel]]  <- function_generator_channel_name
    #   
    #   channel_stimulation_pulse <- channels[1]
    #   channel_function_generator <- channels[2]
    #   
    # }else{
    #   # ch2 was the stimulus
    #   ES_channel <- 2
    #   FunGen_channel <- 1
    #   df_data$Channel[df_data$Channel == channels[ES_channel]]  <- ES_channel_name
    #   df_data$Channel[df_data$Channel == channels[FunGen_channel]]  <- function_generator_channel_name
    #   
    #   channel_stimulation_pulse <- paste("CHAN", channels[2], sep="")
    #   channel_function_generator <- paste("CHAN", channels[1], sep="")
    # }
    # 
    # channel_stimulation_pulse <- ES_channel_name
    # channel_function_generator <- function_generator_channel_name
    
  }else if(length(channels) == 3){
    
    sd_ch1 <- df_data %>% filter(Channel == channels[1] & VPP <1e3 & VPP > 1e-3) %>% summarise(sd(VPP))
    sd_ch2 <- df_data %>% filter(Channel == channels[2] & VPP <1e3 & VPP > 1e-3) %>% summarise(sd(VPP))
    sd_ch3 <- df_data %>% filter(Channel == channels[3] & VPP <1e3 & VPP > 1e-3) %>% summarise(sd(VPP))
    
    # lowest sd -> function generator
    # middle sd -> measurement at resistor
    # highest sd -> stimulator output
    
    ES_channel <- which(c(sd_ch1, sd_ch2, sd_ch3) == max(sd_ch1, sd_ch2, sd_ch3))
    FunGen_channel <- which(c(sd_ch1, sd_ch2, sd_ch3) == min(sd_ch1, sd_ch2, sd_ch3))
    resistor_channel <- which(c(sd_ch1, sd_ch2, sd_ch3) != max(sd_ch1, sd_ch2, sd_ch3) & c(sd_ch1, sd_ch2, sd_ch3) != min(sd_ch1, sd_ch2, sd_ch3))
    
    df_data$Channel[df_data$Channel == channels[ES_channel]]  <- ES_channel_name
    df_data$Channel[df_data$Channel == channels[FunGen_channel]]  <- function_generator_channel_name
    df_data$Channel[df_data$Channel == channels[resistor_channel]]  <- resistor_channel_name
    
  }else{
    print("More than 2 channels used. Please check.")
  }
  
  plotMeasurements(input_data = df_data, output_dir = current_output_dir, vpp_min = 0, vpp_max = 20, vavg_min = -2, vavg_max = 2, plot_title = plot_title)
  
  # Plot waveform data #######################################################
  df_data <- getWaveforms(input_file = current_zip_file)
  
  if(exists("ES_channel_name")){
    df_data$Channel[df_data$Channel == channels[ES_channel]]  <- ES_channel_name
  }
  if(exists("function_generator_channel_name")){
    df_data$Channel[df_data$Channel == channels[FunGen_channel]]  <- function_generator_channel_name
  }
  if(exists("resistor_channel_name")){
    df_data$Channel[df_data$Channel == channels[resistor_channel]]  <- resistor_channel_name
  }

  plotWaveforms(input_data = df_data, output_dir = current_output_dir,
                show_time_in_us = FALSE,
                channel_function_generator = function_generator_channel_name,
                channel_stimulation_pulse = ES_channel_name,
                voltage_limits_of_plot = 10, filter_stim_off = TRUE,
                epsilon_for_filtering = 1)
  
}
