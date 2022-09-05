# Script for analyzing czi images of Live/Dead staining            +++++++++
# Author: Kai Budde
# Created: 2022/07/21
# Last changed: 2022/07/22

# Delete everything in the environment
rm(list = ls())
# close all open plots in RStudio
graphics.off()

# Please adapt the following parameters ####################################
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

# Directory of the images with CellROX measurements
input_directories_livedead <- c(
  "E:/PhD/Daten/ShortStim_ZellBio/LiveDead/Versuch 1/",
  "E:/PhD/Daten/ShortStim_ZellBio/LiveDead/Versuch 2/",
  "E:/PhD/Daten/ShortStim_ZellBio/LiveDead/Versuch 3/",
  "E:/PhD/Daten/ShortStim_ZellBio/LiveDead/Versuch 4/",
  "E:/PhD/Daten/ShortStim_ZellBio/LiveDead/Versuch 5/")

output_dir <- "E:/PhD/Daten/ShortStim_ZellBio/LiveDead/"

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
  input_directories <- c(input_directories_livedead)
  
  for(i in 1:length(input_directories)){
    if(i == 1){
      df_results <- cellPixels::cellPixels(input_dir = input_directories[i],
                                           nucleus_color = "blue",
                                           protein_in_nuc_color = "red",
                                           protein_in_cytosol_color = "green",
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
  input_directories <- c(input_directories_livedead)
  
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
        print(j)
        df_dummy <- readCzi::readCziMetadata(input_file = file_names_czi[j], save_metadata = FALSE)
        df_metadata <- dplyr::bind_rows(df_metadata, df_dummy)
        rm(df_dummy)
      }
    }
    rm(j)
    
    current_output_dir <- paste0(current_input_dir, "output/")
    dir.create(current_output_dir, showWarnings = FALSE)
    
    write.csv(x = df_metadata,
              file = paste(current_output_dir,"summary_metadata.csv", sep=""),
              row.names = FALSE)
    write.csv2(x = df_metadata,
               file = paste(current_output_dir,"summary_metadata_de.csv", sep=""),
               row.names = FALSE)
    
    
  }
  rm(i)
  
}

## Data analysis -----------------------------------------------------------

# Import data
input_directories <- c(input_directories_livedead)
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
# Still unknown dates
df_data$date <- NA
df_data$date <- c(rep("2022-06-20", 60),
                  rep("2022-07-05", 62),
                  rep("2022-07-06", 60),
                  rep("2022-07-07", 60),
                  rep("2022-07-08", 58))
df_data$date <- lubridate::ymd(df_data$date)

# Add column with information of magnification of the objective of the microscope
df_data$magnification <- NA
df_data$magnification <- df_metadata$objective_magnification[match(df_data$fileName, df_metadata$fileName)]

# Add column with information about whether an image was retaken
df_data$image_retake <- FALSE
df_data$image_retake <- grepl(pattern = ".+_[0-9]+_[0-9]+\\.czi$", x = df_data$fileName, ignore.case = TRUE)

# Add column with information about which well the cells were from
df_data$well_number <- NA
df_data$well_number <- suppressWarnings(
  as.numeric(gsub(pattern = "^.+([0-9])+_z.+",
                  replacement = "\\1",
                  x = df_data$fileName)))
# #with retake
# df_data$well_number[df_data$image_retake] <- suppressWarnings(
#   as.numeric(gsub(pattern = ".+([0-9]+)_[0-9]+_[0-9]\\.czi$",
#                   replacement = "\\1",
#                   x = df_data$fileName[df_data$image_retake])))


# Add column with information about which image number the cells were from
df_data$image_number <- NA

#without retake
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
  grepl(pattern = "stm", x = df_data$fileName, ignore.case = TRUE)] <- "Stimulation"
# df_data$experiment_group[
#   grepl(pattern = " poskon", x = df_data$fileName, ignore.case = TRUE)] <- "PositiveControl"
df_data$experiment_group[
  grepl(pattern = "kon", x = df_data$fileName, ignore.case = TRUE)] <- "Control"


# Add column with information about which experiment_number an experiment belongs to (n)
df_data$experiment_number <- NA
number_of_dates <- length(unique(df_data$date))
number_of_wells <- length(unique(df_data$well_number[!is.na(df_data$well_number)]))
number_of_independent_experiments <- number_of_dates * number_of_wells

df_data$experiment_number <-  (match(df_data$date, unique(df_data$date))-1) * number_of_wells +
  df_data$well_number

# Calculate fractions of dead/live cells
df_data$fraction_dead_live <- df_data$number_of_nuclei_with_second_protein /
  df_data$number_of_nuclei


# Reorder columns
df_data <- df_data %>% 
  dplyr::relocate(c("experiment_group", "experiment_number", "well_number", "image_number", "image_retake", "magnification",
                    "fraction_dead_live", "number_of_nuclei", "number_of_nuclei_with_second_protein"),
                  .after = "fileName")
  # dplyr::relocate(c("date", "experiment_group", "experiment_number", "well_number", "image_number", "image_retake", "magnification" ),
  #                 .after = "fileName")




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


filter_experiment_number <- unique(df_data$experiment_number)

for(i in 0:length(filter_experiment_number)){
  
  if(i == 0){
    # Do not filter
    df_data_dummy <- df_data
    append_text <- ""
  }else{
    # Filter
    df_data_dummy <- df_data %>% 
      dplyr::filter(experiment_number == i)
    append_text <- paste0("(Experiment ", i, ")")
  }
  
  # Plotting ###############################################################
  output_dir_plots <- paste0(output_dir, "plots/")
  dir.create(output_dir_plots, showWarnings = FALSE)
  
  plot_title_nuc <- "Live-Dead"
  
  # Number (no) of nuclei per image (blue) ---------------------------------
  
  df_data_dummy$experiment_group <- factor(df_data_dummy$experiment_group, levels = c("Control", "Stimulation"))
  
  plot_no_nuclei <- ggplot(df_data_dummy, aes(x=experiment_group, y=number_of_nuclei)) +
    geom_violin() +
    # stat_boxplot(geom ='errorbar', width = 0.3) +
    geom_boxplot(width=0.05) +
    stat_summary(fun=mean, geom="point", size = 3, shape=23, color="black", fill="black") +
    facet_wrap(.~ magnification, scale="free") +
    # geom_jitter(color="black", size=0.5, alpha=0.9) +
    #ylim(0,20) +
    theme_bw(base_size = 24) +
    theme(axis.title.y = element_text(size=18)) +
    xlab("") +
    ylab("Number of nuclei per image") +
    ggtitle(paste0("Number of nuclei (", plot_title_nuc, ") ", append_text)) +
    scale_fill_discrete(name = "Experiment Group",  labels = c("Control", "Stimulation"))
  
  # print(plot_no_nuclei)
  
  ggsave(filename = paste(output_dir_plots, plot_title_nuc, "_number_of_nuclei_violin", append_text, ".pdf", sep=""),
         width = 297, height = 210, units = "mm")
  ggsave(filename = paste(output_dir_plots, plot_title_nuc, "_number_of_nuclei_violin", append_text, ".png", sep=""),
         width = 297, height = 210, units = "mm")
  
  plot_no_nuclei <- ggplot(df_data_dummy, aes(x=experiment_group,
                                              y=number_of_nuclei,
                                              fill=experiment_group)) +
    geom_boxplot() +
    stat_summary(fun = mean, geom="point", size = 3, shape=23, color="black", fill="black") +
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
    scale_fill_discrete(name = "Experiment Group",  labels = c("Control", "Stimulation"))
  
  #print(plot_no_nuclei)
  
  ggsave(filename = paste(output_dir_plots, plot_title_nuc, "_number_of_nuclei", append_text, ".pdf", sep=""),
         width = 297, height = 210, units = "mm")
  ggsave(filename = paste(output_dir_plots, plot_title_nuc, "_number_of_nuclei", append_text, ".png", sep=""),
         width = 297, height = 210, units = "mm")
  
  # Ratio of dead/live cell (no) of nuclei per image (blue) ----------------
  plot_deadlive_fraction <- ggplot(df_data_dummy, aes(x=experiment_group, y=fraction_dead_live)) +
    geom_violin() +
    # stat_boxplot(geom ='errorbar', width = 0.3) +
    geom_boxplot(width=0.05) +
    stat_summary(fun=mean, geom="point", size = 3, shape=23, color="black", fill="black") +
    facet_wrap(.~ magnification, scale="free") +
    # geom_jitter(color="black", size=0.5, alpha=0.9) +
    ylim(0,0.25) +
    theme_bw(base_size = 24) +
    theme(axis.title.y = element_text(size=18)) +
    xlab("") +
    ylab("Fraction of dead cells") +
    ggtitle(paste0("Fraction of dead cells (", plot_title_nuc, ") ", append_text)) +
    scale_fill_discrete(name = "Experiment Group",  labels = c("Control", "Stimulation"))
  
  # print(plot_deadlive_fraction)
  
  ggsave(filename = paste(output_dir_plots, plot_title_nuc, "_fraction_violin", append_text, ".pdf", sep=""),
         width = 297, height = 210, units = "mm")
  ggsave(filename = paste(output_dir_plots, plot_title_nuc, "_fraction_violin", append_text, ".png", sep=""),
         width = 297, height = 210, units = "mm")
  
  plot_deadlive_fraction <- ggplot(df_data_dummy, aes(x=experiment_group,
                                              y=fraction_dead_live,
                                              fill=experiment_group)) +
    geom_boxplot() +
    stat_summary(fun = mean, geom="point", size = 3, shape=23, color="black", fill="black") +
    facet_wrap(.~ magnification, scale="free") +
    theme_bw(base_size = 24) +
    theme(axis.title.y = element_text(size=18),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank()) +
    ylim(0,0.25) +
    ylab("Fraction of dead cells") +
    xlab("") +
    ggtitle(paste0("Fraction of dead cells (", plot_title_nuc, ") ", append_text)) +
    scale_fill_discrete(name = "Experiment Group",  labels = c("Control", "Stimulation"))
  
  #print(plot_deadlive_fraction)
  
  ggsave(filename = paste(output_dir_plots, plot_title_nuc, "_fraction", append_text, ".pdf", sep=""),
         width = 297, height = 210, units = "mm")
  ggsave(filename = paste(output_dir_plots, plot_title_nuc, "_fraction", append_text, ".png", sep=""),
         width = 297, height = 210, units = "mm")
  
  
}




# Get fraction of dead/live for stim/control looking at every experiment ---
# number separately                                                      ---

df_result <- df_data %>% 
  dplyr::group_by(experiment_number, experiment_group) %>% 
  dplyr::summarise(fraction_mean = mean(fraction_dead_live))

write.csv(x = df_result,
          file = paste(output_dir_data,"fraction_dead_live_eaggregated.csv", sep=""),
          row.names = FALSE)
write.csv2(x = df_result,
           file = paste(output_dir_data,"fraction_dead_live_aggregated_de.csv", sep=""),
           row.names = FALSE)



# Calculate ratios
df_result_ratio <- df_result %>% 
  dplyr::group_by(experiment_number) %>% 
  dplyr::mutate(ratio = 1/((sum(fraction_mean)/fraction_mean)-1))


df_result_ratio <- df_result_ratio %>% 
  dplyr::filter(experiment_group == "Stimulation") %>% 
  dplyr::select(experiment_number, ratio)

mean_ratio <- mean(df_result_ratio$ratio)
sd_ratio <- sd(df_result_ratio$ratio)
se_ratio <- sd_ratio/sqrt(number_of_independent_experiments)

# Ratio (Stim/Control) of all fluorescence
plot_ratios <- ggplot() +
  geom_point(aes(x = 0, y = mean_ratio), size = 6) +
  geom_errorbar(aes(x = 0, ymin=mean_ratio-se_ratio,
                    ymax=mean_ratio+se_ratio),
                size = 1.5, width=.2) +
  ylim(0,2) +
  # xlim(-2,2) +
  theme_bw(base_size = 24) +
  theme(axis.title.y=element_text(size=18),
        #axis.text.x = element_blank(),
        # axis.ticks.x = element_blank(),
        # axis.text.x = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.title.x=element_blank()) +
  scale_x_discrete(labels=c("green_in_cell_per_nuc" = "Green\nin cell",
                            "green_in_nuc_per_nuc" = "Green\nin nucleus",
                            "red_in_cell_per_nuc" = "Red\nin cell",
                            "red_in_nuc_per_nuc" = "Red\nin nucleus")) +
  ggtitle("Ratio of dead/alive fractions (Stimulation/Control)") +
  geom_hline(yintercept=1.0, linetype="dashed", size=2) +
  ylab("Ratio of dead/alive fractions (mean+-SE)")

#print(plot_ratios)

ggsave(filename = paste(output_dir_plots, "deadalive_ratios_StimOverControl_combined.pdf", sep=""),
       width = 297, height = 210, units = "mm")
ggsave(filename = paste(output_dir_plots,"deadalive_ratios_StimOverControl_combined.png", sep=""),
       width = 297, height = 210, units = "mm")


