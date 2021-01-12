rm(list = ls()) 

#Plotting libs
library(ggplot2)
library(ggthemes)
#Used to transform dataframes befor plotting
library(reshape)

col_HRV4T_blue <-rgb(20/256, 46/256, 81/256, 1)
col_HRV4T_lightblue <-rgb(52/256, 176/256, 227/256, 1)
col_HRV4T_yellow <-rgb(251/256, 207/256, 48/256, 1)
col_HRV4T_gray <-rgb(151/256, 151/256, 151/256, 1)
col_HRV4T_gray_alpha <-rgb(151/256, 151/256, 151/256, 0.3)

#function to compute rMSSD
hrv_features_rmssd <- function(values)
{
  valuesDiff <- diff(values)
  return (sqrt(mean(valuesDiff^2, na.rm = TRUE)))
}

output_to_pdf = TRUE

#EDIT based on your machine settings
files_path_root <- paste("~/Dropbox/R workspace/github/hera leto two/", sep = "")
files_path_data <- paste(files_path_root, "data/", sep = "")
files_path <- paste(files_path_root, "figures/", sep = "")
setwd(files_path_root)
source(paste(files_path_root, "multiplot.R", sep = ""))

subjects <- c("001")
df_rr <- data.frame()
for(index_subject in 1:length(subjects))
{
  curr_subject <- subjects[index_subject]
  
  #Load reference (Polar H7 data)
  rr_h7 = read.csv(paste(files_path_data, curr_subject, "/polar.csv", sep = ""), header=TRUE)
  rr_h7 <- rr_h7[, c(2:3)] 
  names(rr_h7) <- c("rr", "since_start")
  rr_h7$since_start <- rr_h7$since_start  / 1000 / 1.024 # convert to seconds, this needs to be done only for Polar sensors 
  rr_h7[, "sensor"] <- "Polar H7"
  
  #Load Hera Leto Two
  rr_hera = read.csv(paste(files_path_data, curr_subject, "/hera.csv", sep = ""), header=TRUE)
  rr_hera <- rr_hera[, c(2:3)] 
  names(rr_hera) <- c("rr", "since_start")
  rr_hera$since_start <- rr_hera$since_start  / 1000 #convert to seconds
  rr_hera[, "sensor"] <- "Hera Leto Two"
  
  # sync
  if(curr_subject == '001')
  {
    rr_hera <- rr_hera[15:nrow(rr_hera), ]
    rr_hera$since_start <- rr_hera$since_start - rr_hera$since_start[1]
  }

  #Create data frame to plot using ggplot
  df_rr_curr_subj <- rbind(rr_h7[, c('rr', 'since_start', 'sensor')], rr_hera[, c('rr', 'since_start', 'sensor')])
  df_rr_curr_subj[, "subject_ID"] <- curr_subject
  df_rr <- rbind(df_rr, df_rr_curr_subj)
}

#Segment windows for HRV computation and plotting (1 minute)
min_window <- 1
max_window <- round(max(df_rr$since_start)/60)
window_size <- 60
max_window <- max_window/(window_size/60)
df_rr[, "window_min"] <- NA
for(index_window_min in min_window:max_window)
{
  df_rr[df_rr$since_start >= (window_size*(index_window_min-1)) & 
          df_rr$since_start < (window_size*index_window_min), "window_min"] <- index_window_min
}
df_rr <- df_rr[!is.na(df_rr$window_min), ]

#Compute features over segmented windows
df_features <- data.frame()
for(index_subject in 1:length(subjects))
{
  curr_subject <- subjects[index_subject]
  curr_subject_data <- df_rr[df_rr$subject_ID == curr_subject, ]
  
  for(index_window_min in 1:max_window)
  {
    #Reference feature
    curr_window_h7 <- curr_subject_data[!is.na(curr_subject_data$window_min) & 
                                             curr_subject_data$window_min == index_window_min &
                                             curr_subject_data$sensor == "Polar H7", "rr"]
    rMSSD_h7 <- round(hrv_features_rmssd(curr_window_h7), 1)
    
    #Hera
    curr_window_hera <- curr_subject_data[!is.na(curr_subject_data$window_min) & 
                                          curr_subject_data$window_min == index_window_min &
                                          curr_subject_data$sensor == "Hera Leto Two", "rr"]
    rMSSD_hera <- round(hrv_features_rmssd(curr_window_hera), 1)
    
    curr_features <- data.frame(rMSSD_h7, rMSSD_hera)
    names(curr_features) <- c("Polar H7", "Hera Leto Two")
    curr_features[, "window"] <- index_window_min 
    curr_features[, "subject_ID"] <- curr_subject 
    df_features <- rbind(df_features, curr_features)
  }
}

#Plot data, RR intervals first (synch is not perfect but signals overlap decently, won't be shifting or aligning them any further)
for(index_subject in 1:length(subjects))
{
  curr_subject <- subjects[index_subject]
  curr_subject_data <- df_rr[df_rr$subject_ID == curr_subject, ]
  
  #rr intervals 
  p1 <- ggplot(curr_subject_data, aes(since_start, rr, col = sensor)) +
    geom_line(size = 1) +
    facet_wrap(sensor~window_min, scale = "free_x", ncol = max_window-min_window+1) +
    ggtitle(paste("Comparison: Polar H7 vs Hera Leto Two (RR intervals) - Subject", curr_subject)) +
    theme_minimal() + 
    scale_colour_manual(values = c(col_HRV4T_blue, col_HRV4T_yellow, col_HRV4T_lightblue)) +
    scale_fill_manual(values = c(col_HRV4T_blue, col_HRV4T_yellow, col_HRV4T_lightblue)) +
    theme(panel.background = element_rect(fill = 'white', colour = col_HRV4T_gray_alpha)) + 
    labs(col = 'Sensor') +
    theme(legend.position="none") +
    scale_x_continuous("Time (seconds)")+
    scale_y_continuous("RR (ms)", limits = c(650, 1350))
  if(output_to_pdf)
  {
    pdf(paste(files_path,"fig_rr_", curr_subject, ".pdf", sep=""), width=16, height=10)
  }
  multiplot(p1)
  if(output_to_pdf)
  {
    dev.off()
  }
}

#Plot features (rMSSD) for all sensors, subjects and windows
df_rmssd <- melt(df_features[, c("Polar H7", "Hera Leto Two", "window", "subject_ID")], id = c("window", "subject_ID"))
names(df_rmssd)[3:4] <- c("Sensor", "rMSSD")

p1 <- ggplot(df_rmssd[!is.na(df_rmssd$rMSSD), ], aes(window, rMSSD, fill = Sensor)) +
  geom_bar(stat = 'identity', position = 'dodge') +
  facet_wrap(~subject_ID, scales = 'free') +
  ggtitle("Comparison: Polar H7 vs Hera Leto Two(rMSSD in ms)") +
  xlab("Time window") +
  scale_y_continuous('ms', limits = c(0, 180)) +
  theme_minimal() + 
  scale_colour_manual(values = c(col_HRV4T_blue, col_HRV4T_yellow, col_HRV4T_lightblue)) +
  scale_fill_manual(values = c(col_HRV4T_blue, col_HRV4T_yellow, col_HRV4T_lightblue)) +
  theme(panel.background = element_rect(fill = 'white', colour = col_HRV4T_gray_alpha)) + 
  labs(fill = 'Sensor')
if(output_to_pdf)
{
  pdf(paste(files_path,"fig_rmssd_all.pdf", sep=""), width=10, height=5)
}
multiplot(p1)
if(output_to_pdf)
{
  dev.off()
}

