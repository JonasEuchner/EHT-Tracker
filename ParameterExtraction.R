# ParameterExtraction_v2.0.1

# Check if all packages are available, and install it if not
is_readr_available <- require("readr")
if (!is_readr_available) {
  install.packages("readr")
}
is_tidyverse_available <- require("tidyverse")
if (!is_tidyverse_available) {
  install.packages("tidyverse")
}
is_rstudioapi_available <- require("rstudioapi")
if (!is_rstudioapi_available) {
  install.packages("rstudioapi")
}
is_svDialogs_available <- require("svDialogs")
if (!is_svDialogs_available) {
  install.packages("svDialogs")
}
is_signal_available <- require("signal")
if (!is_signal_available) {
  install.packages("signal")
}
library(signal)
library(readr)
library(tidyverse)
library(rstudioapi)
library(svDialogs)

if(!is.null(dev.list())) dev.off()

dir<-selectDirectory(
  caption = "Select Directory of ImageJ EHT-tracking results (*.txt)",
  label = "Select",
  path = getActiveProject()
)

setwd(dir)
FileList<-list.files(pattern="*.txt")

FrameRate <- as.numeric(dlgInput("Framerate [Frames per second]", "53.28")$res)
pixelsize <- as.numeric(dlgInput("Pixelsiz [mm/px]", "0.0036205")$res)
DisplacementConversionFactor<- as.numeric(dlgInput("Conversion factor [mN/mm]: see https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0266834", "0.26")$res)
pillarCount <- as.numeric(dlgInput("Number of pillars per EHT [1 or 2]", "1")$res)

# Create empty data frames for Trace extraction
Trace<-data.frame(dist=NA)
output<-data.frame(rep(NA,length(FileList)),rep(NA,length(FileList)),rep(NA,length(FileList)),rep(NA,length(FileList)),rep(NA,length(FileList)),rep(NA,length(FileList)),rep(NA,length(FileList)),rep(NA,length(FileList)))
colnames(output)<-c("Peak to peak [s]","Contraction duration [ms]","Time to peak [ms]","Time to peak (5% to peak) [ms]","Relaxation time [ms]","Relaxation time (peak to 5%) [ms]","Contraction amplitude [mm]","Contraction force [µN]")

# Iterate over each file in the 'FileList'
for (i in 1:length(FileList)) {
  Pillar <- read_delim(paste(dir,FileList[i],sep="/"),";", escape_double = FALSE, trim_ws = TRUE)
  
  # Find resting x and y position
  medianX<-median(Pillar$xPos)
  posXQ<-quantile(Pillar$xPos, c(0.01, 0.99))
  minX<-as.numeric(posXQ[1])
  maxX<-as.numeric(posXQ[2])
  medianY<-median(Pillar$yPos)
  posYQ<-quantile(Pillar$yPos, c(0.01, 0.99))
  minY<-as.numeric(posYQ[1])
  maxY<-as.numeric(posYQ[2])
  
  if (abs(maxX-medianX)<abs(minX-medianX)) {
    xResting <- maxX
  } else if (abs(maxX-medianX)>abs(minX-medianX)){
    xResting <- minX
  } else {
    xResting = mean(c(maxX, minX))
  }
  
  if (abs(maxY-medianY)<abs(minY-medianY)) {
    yResting <- maxY
  } else if (abs(maxY-medianY)>abs(minY-medianY)){
    yResting <- minY
  } else {
    yResting = mean(c(maxY, minY))
  }
  
  xDistance<-xResting-Pillar$xPos
  yDistance<-yResting-Pillar$yPos
  Pillar$dist<-sqrt(xDistance^2+yDistance^2)

  # Apply shape-preserving smoothing filter
  smoothed_savgol <- sgolayfilt(Pillar$dist, p = 3, n = floor(floor(FrameRate)/3)*2+1)
  window_size=5
  split_size <- floor(length(smoothed_savgol) / window_size)
  subsets <- split(smoothed_savgol, ceiling(seq_along(1:length(smoothed_savgol)) / split_size))
  
  # Correct baseline
    baseline_values <- numeric(window_size)
  for (j in 1:window_size) {
    baseline_values[j] <- quantile(subsets[[j]], probs = 0.005)
  }
  
  center_indices <- sapply(1:window_size, function(j) {
    start_idx <- (j - 1) * split_size + 1
    end_idx <- min(j * split_size, length(smoothed_savgol))
    return(floor((start_idx + end_idx) / 2)) 
  })
  
  poly_fit <- lm(baseline_values ~ stats::poly(center_indices, degree=window_size-1, raw=TRUE))
  predicted_baseline <- predict(poly_fit, newdata = data.frame(center_indices = 1:length(smoothed_savgol)))
  corrected_dist <- smoothed_savgol-predicted_baseline
  
  Trace_corrected=data.frame(corrected_dist)
  names(Trace_corrected)[1] <- "Displacement"
  occurences<-data.frame(table(unlist(round(Trace_corrected$Displacement,digits=1))))
  corFactor<-match(max(occurences$Freq),occurences$Freq)
  Trace_corrected<-Trace_corrected-as.numeric(levels(occurences$Var1))[occurences$Var1][corFactor]
  
  # Peak detection
  find_all_rise <- function(data) {
    max_value <- max(data)
    threshold <- 0.3 * max_value
    
    rise_indices <- numeric(0)
    
    for (j in 6:(length(data) - 5)) {
      if (all(data[(j - 5):(j - 1)] < threshold) && sum(data[j:(j + 2)] > threshold) >= 3) {
        rise_indices <- c(rise_indices, j)  
      }
    }
    
    if (length(rise_indices) > 0) {
      return(rise_indices)
    } else {
      return(NULL)
    }
  }
  
  find_all_fall <- function(data) {
    max_value <- max(data)
    threshold <- 0.3 * max_value  
    
    fall_indices <- numeric(0)  
    
    for (j in 6:(length(data) - 5)) {
      if (sum(data[(j-2):j] > threshold) >= 3 && all(data[(j +1):(j +5)] < threshold)) {
        fall_indices <- c(fall_indices, j) 
      }
    }
    
    if (length(fall_indices) > 0) {
      return(fall_indices)
    } else {
      return(NULL)
    }
  }
  
  rise_indices <- find_all_rise(Trace_corrected$Displacement)
  fall_indices <- find_all_fall(Trace_corrected$Displacement)

  
  fall_indices_to_keep <- fall_indices[sapply(fall_indices, function(fall_value) {
    any(rise_indices < fall_value)
  })]
  
  rise_indices_to_keep <- rise_indices[sapply(rise_indices, function(rise_value) {
    any(fall_indices > rise_value)
  })]
  
  detected_peaks <- data.frame(rise_indices = rise_indices_to_keep, fall_indices = fall_indices_to_keep)
  
  max_values <- numeric(nrow(detected_peaks))
  max_values_idx <- numeric(nrow(detected_peaks))
  min_values_pre <- numeric(nrow(detected_peaks))
  min_values_pre_idx<- numeric(nrow(detected_peaks))
  min_values_post <- numeric(nrow(detected_peaks))
  min_values_post_idx<- numeric(nrow(detected_peaks))
  
  # find maxima
  for (j in 1:nrow(detected_peaks)) {
    max_values[j] <- NA
    rise_idx <- detected_peaks$rise_indices[j]
    fall_idx <- detected_peaks$fall_indices[j]
    peak_segment <- Trace_corrected$Displacement[rise_idx:fall_idx]
    max_values[j] <- max(peak_segment)
    max_values_idx[j] <- which.max(peak_segment)+rise_idx-1
  }
  
  #find pre contractionk minimum
  for (k in 1:(nrow(detected_peaks))) {
    rise_idx <- detected_peaks$rise_indices[k]
    min_values_pre_idx[k] <- NA
    min_values_pre[k] <- NA
    
    if (k == 1) {
      if (fall_indices[1]<rise_idx) {
        bl_segment <- Trace_corrected$Displacement[fall_indices[1]:rise_idx]
        peak_height <- max_values[k]
        threshold <- 0.02 * peak_height
        
        min_val <- min(bl_segment)
        max_val <- max(bl_segment)
        bins <- cut(bl_segment, breaks=100, include.lowest=TRUE)
        bin_counts <- table(bins)
        most_frequent_bin <- names(bin_counts)[which.max(bin_counts)]
        upper_value_str <- sub(".*\\,(.*\\])", "\\1", most_frequent_bin)
        upper_value <- as.numeric(substr(upper_value_str, 1, nchar(upper_value_str) - 1))
        
        for (j in length(bl_segment):1) {
          if (bl_segment[j] < threshold && bl_segment[j]< upper_value) {
            min_values_pre_idx[k] <- rise_idx-(length(bl_segment)-j)
            min_values_pre[k]<-bl_segment[j]
            break 
          }
        }
      }
    } else {
      fall_idx <- detected_peaks$fall_indices[k-1]
      bl_segment <- Trace_corrected$Displacement[fall_idx:rise_idx]
      peak_height <- max_values[k]
      threshold <- 0.02 * peak_height
      
      min_val <- min(bl_segment)
      max_val <- max(bl_segment)
      bins <- cut(bl_segment, breaks=100, include.lowest=TRUE)
      bin_counts <- table(bins)
      most_frequent_bin <- names(bin_counts)[which.max(bin_counts)]
      upper_value_str <- sub(".*\\,(.*\\])", "\\1", most_frequent_bin)
      upper_value <- as.numeric(substr(upper_value_str, 1, nchar(upper_value_str) - 1))
      
      for (j in length(bl_segment):1) {
        if (bl_segment[j] < threshold && bl_segment[j]< upper_value) {
          min_values_pre_idx[k] <- rise_idx-(length(bl_segment)-j)
          min_values_pre[k]<-bl_segment[j]
          break 
        }
      }
    }
  }
  
  #find post contractionk minimum
  for (k in 1:(nrow(detected_peaks))) {
    min_values_post_idx[k] <- NA
    min_values_post[k] <- NA
    fall_idx <- detected_peaks$fall_indices[k]
    if (k == nrow(detected_peaks)) {
      if (fall_idx<rise_indices[length(rise_indices)]) {
        bl_segment <- Trace_corrected$Displacement[fall_idx:rise_indices[length(rise_indices)]]
        peak_height <- max_values[k]
        threshold <- 0.02 * peak_height
        
        min_val <- min(bl_segment)
        max_val <- max(bl_segment)
        bins <- cut(bl_segment, breaks=100, include.lowest=TRUE)
        bin_counts <- table(bins)
        most_frequent_bin <- names(bin_counts)[which.max(bin_counts)]
        upper_value_str <- sub(".*\\,(.*\\])", "\\1", most_frequent_bin)
        upper_value <- as.numeric(substr(upper_value_str, 1, nchar(upper_value_str) - 1))
        
        for (j in 1:length(bl_segment)) {
          if (bl_segment[j] < threshold && bl_segment[j]< upper_value){
            min_values_post_idx[k] <- fall_idx+j-1
            min_values_post[k]<-bl_segment[j]
            break 
          }
        }
      }
    } else {
      rise_idx <- detected_peaks$rise_indices[k+1]
      bl_segment <- Trace_corrected$Displacement[fall_idx:rise_idx]
      peak_height <- max_values[k]
      threshold <- 0.02 * peak_height
      
      min_val <- min(bl_segment)
      max_val <- max(bl_segment)
      bins <- cut(bl_segment, breaks=100, include.lowest=TRUE)
      bin_counts <- table(bins)
      most_frequent_bin <- names(bin_counts)[which.max(bin_counts)]
      upper_value_str <- sub(".*\\,(.*\\])", "\\1", most_frequent_bin)
      upper_value <- as.numeric(substr(upper_value_str, 1, nchar(upper_value_str) - 1))
      
      for (j in 1:length(bl_segment)) {
        if (bl_segment[j] < threshold && bl_segment[j]< upper_value){
          min_values_post_idx[k] <- fall_idx+j-1
          min_values_post[k]<-bl_segment[j]
          break 
        }
      }
    }
  }

  #store all peak information
  peaks <- data.frame(max_value = (max_values-0.5*min_values_pre-0.5*min_values_post),
                      max_idx = max_values_idx,
                      min_pre_idx = min_values_pre_idx,
                      min_post_idx = min_values_post_idx
                      )
  peaks <- peaks[complete.cases(peaks), ]

  Frame<-c(1:length(Trace_corrected$Displacement))
  threshold_level <- 0.05
  threshold <- threshold_level * Trace_corrected$Displacement[peaks[1,2]]
  peaks <- cbind(peaks,rep(NA, nrow(peaks)),rep(NA, nrow(peaks)))
  for (m in 1:(length(peaks[,2]))) {
    for (n in seq(peaks[m,2] - 1, 1)) {
      if (Trace_corrected$Displacement[n] <= threshold) {
        peaks[m,5] <- n
        break
      }
    }
    for (n in seq(peaks[m,2] + 1, length(Trace_corrected$Displacement))) {
      if (Trace_corrected$Displacement[n] <= threshold) {
        peaks[m,6] <- n
        break
      }
    }
  }
  if (is.na(peaks[1,5])) {
    peaks[m,5]<-peaks[m,3]
  }
  if (is.na(peaks[m,6])) {
    peaks[m,6]<-peaks[m,4]
  }
 
  
  #Calculate all important parameters per contraction
  intermediate<-data.frame("temp"=rep(NA,length(peaks[,1])))
  
  intermediate$`PeakPos [Frame]` <-sort(peaks[,2])
  intermediate$`Time to peak [ms]`<-(intermediate$`PeakPos [Frame]`-sort(peaks[,3]))/FrameRate*1000
  intermediate$`Time to peak (5% to peak) [ms]`<-(intermediate$`PeakPos [Frame]`-sort(peaks[,5]))/FrameRate*1000
  intermediate$`Relaxation time [ms]`<-(sort(peaks[,4])-intermediate$`PeakPos [Frame]`)/FrameRate*1000
  intermediate$`Relaxation time (peak to 5%) [ms]`<-(sort(peaks[,6])-intermediate$`PeakPos [Frame]`)/FrameRate*1000  
  intermediate$`Contraction duration [ms]`<-(sort(peaks[,4])-sort(peaks[,3]))/FrameRate*1000
  intermediate$`Contraction duration (5% to 5%) [ms]`<-(sort(peaks[,6])-sort(peaks[,5]))/FrameRate*1000  
  intermediate$`Peak to peak [s]`<-(intermediate$`PeakPos [Frame]` - lag(intermediate$`PeakPos [Frame]`))/FrameRate
  intermediate$`Contraction amplitude [mm]`<-peaks[,1]*pixelsize
  intermediate$`Contraction force [µN]`<-intermediate$`Contraction amplitude [mm]`*DisplacementConversionFactor*1000
  
  intermediate$temp<- NULL
  
  subDir=paste0(str_split(FileList[i],"[.]")[[1]][1])
  if (!file.exists(subDir)){
    dir.create(file.path(dir, subDir))
  }
  
  # Save all peak-related information per trace
  name<-paste0(subDir,"/",strsplit(FileList[i], ".txt")[[1]][1],"_allPeaks.csv")
  write.csv(intermediate, paste(dir,name,sep="/"), row.names=TRUE)
  
  if(length(intermediate$`PeakPos [Frame]`)>3){
    intermediate_adj <- intermediate[2:(nrow(intermediate)- 1),]
  }
  
  Q_RT <- quantile(intermediate_adj$`Relaxation time [ms]`, probs=c(.25, 0.5,.75), na.rm = FALSE)
  iqr_RT <- IQR(intermediate_adj$`Relaxation time [ms]`)
  Q_TP <- quantile(intermediate_adj$`Time to peak [ms]`, probs=c(.25, 0.5,.75), na.rm = FALSE)
  iqr_TP <- IQR(intermediate_adj$`Time to peak [ms]`)
  Q_PTP <- quantile(intermediate_adj$`Peak to peak [s]`, probs=c(.25, 0.5,.75), na.rm = FALSE)
  iqr_PTP <- IQR(intermediate_adj$`Peak to peak [s]`)
  Q_CA <- quantile(intermediate_adj$`Contraction amplitude [mm]`, probs=c(.25, 0.5,.75), na.rm = FALSE)
  iqr_CA <- IQR(intermediate_adj$`Contraction amplitude [mm]`)
  
  intermediate_adj_temp<- subset(intermediate_adj, intermediate_adj$`Relaxation time [ms]` > (Q_RT[1]*0.8 - 1.5*iqr_RT) & intermediate_adj$`Relaxation time [ms]`< (Q_RT[3]*1.25+1.5*iqr_RT) & +
                                   intermediate_adj$`Time to peak [ms]` > (Q_TP[1]*0.8 - 1.5*iqr_TP) & intermediate_adj$`Time to peak [ms]`< (Q_TP[3]*1.25+1.5*iqr_TP) & +
                                   intermediate_adj$`Peak to peak [s]` > (Q_PTP[1]*0.8 - 1.5*iqr_PTP) & intermediate_adj$`Peak to peak [s]`< (Q_PTP[3]*1.25+1.5*iqr_PTP) & +
                                   intermediate_adj$`Contraction amplitude [mm]` > (Q_CA[1]*0.8 - 1.5*iqr_CA) & intermediate_adj$`Contraction amplitude [mm]`< (Q_CA[3]*1.25+1.5*iqr_CA))
  
  if(length(intermediate_adj_temp$`PeakPos [Frame]`)>3){
    intermediate_adj <- intermediate_adj_temp
    name<-paste0(subDir,"/",strsplit(FileList[i], ".txt")[[1]][1],"_cleanedPeaks.csv")
    write.csv(intermediate_adj, paste(dir,name,sep="/"), row.names=TRUE)
  }
  
  
  setwd(dir)
  png(file=paste0(subDir,"/",str_split(FileList[i],"[.]")[[1]][1],"_ForceOverTime_raw.png"), width=1200, height=600)
  plot(Frame/FrameRate,Pillar$Distance*pixelsize*DisplacementConversionFactor*1000,type="l",xlab="Time [s]", ylab="Contraction Force [µN]",main=strsplit(FileList[i], ".txt")[[1]][1])+abline(h = 0, col="red")+abline(v = peaks[1:(length(peaks[,1])),2]/FrameRate, col="blue")+abline(v = peaks[1:(length(peaks[,1])),3]/FrameRate, col="cyan")+abline(v = peaks[1:(length(peaks[,1])),4]/FrameRate, col="green")+abline(v = intermediate_adj$PeakPos..Frame./FrameRate, col="magenta",lwd =2)
  dev.off()
  
  png(file=paste0(subDir,"/",str_split(FileList[i],"[.]")[[1]][1],"_ForceOverTime_corrected.png"), width=1200, height=600)
  plot(Frame/FrameRate,Trace_corrected$Displacement*pixelsize*DisplacementConversionFactor*1000,type="l",xlab="Time [s]", ylab="Contraction Force [µN]",main=strsplit(FileList[i], ".txt")[[1]][1])+abline(h = 0, col="red")+abline(v = peaks[1:(length(peaks[,1])),2]/FrameRate, col="blue")+abline(v = peaks[1:(length(peaks[,1])),3]/FrameRate, col="cyan")+abline(v = peaks[1:(length(peaks[,1])),4]/FrameRate, col="green")+abline(v = intermediate_adj$PeakPos..Frame./FrameRate, col="magenta",lwd =2)
  dev.off()
  
  png(file=paste0(subDir,"/",str_split(FileList[i],"[.]")[[1]][1],"_ForceOverTime_corrected_5topeakto5.png"), width=1200, height=600)
  plot(Frame/FrameRate,Trace_corrected$Displacement*pixelsize*DisplacementConversionFactor*1000,type="l",xlab="Time [s]", ylab="Contraction Force [µN]",main=strsplit(FileList[i], ".txt")[[1]][1])+abline(h = 0, col="red")+abline(v = peaks[1:(length(peaks[,1])),2]/FrameRate, col="blue")+abline(v = peaks[1:(length(peaks[,1])),5]/FrameRate, col="cyan")+abline(v = peaks[1:(length(peaks[,1])),6]/FrameRate, col="green")+abline(v = intermediate_adj$PeakPos..Frame./FrameRate, col="magenta",lwd =2)
  dev.off()
  
  svg(file=paste0(subDir,"/",str_split(FileList[i],"[.]")[[1]][1],"_ForceOverTime_corrected_5topeakto5.svg"), width=8, height=4)
  plot(Frame/FrameRate,Trace_corrected$Displacement*pixelsize*DisplacementConversionFactor*1000,type="l",xlab="Time [s]", ylab="Contraction Force [µN]",main=strsplit(FileList[i], ".txt")[[1]][1])+abline(h = 0, col="red")+abline(v = peaks[1:(length(peaks[,1])),2]/FrameRate, col="blue")+abline(v = peaks[1:(length(peaks[,1])),5]/FrameRate, col="cyan")+abline(v = peaks[1:(length(peaks[,1])),6]/FrameRate, col="green")+abline(v = intermediate_adj$PeakPos..Frame./FrameRate, col="magenta",lwd =2)
  dev.off()
  
  # Calculate and store summary statistics in the 'output' data frame for each pillar (remove first and last)
  output$`Peak to peak [s]`[i]<-signif(mean(intermediate_adj$`Peak to peak [s]`, na.rm=TRUE),3)
  output$`Peak to peak - standard dev [s]`[i]<-signif(sd(intermediate_adj$`Peak to peak [s]`, na.rm=TRUE),3)
  
  output$`Contraction duration [ms]`[i]<-round(mean(intermediate_adj$`Contraction duration [ms]`, na.rm=TRUE),0)
  output$`Contraction duration - standard dev [ms]`[i]<-round(sd(intermediate_adj$`Contraction duration [ms]`, na.rm=TRUE),0)
  
  output$`Time to peak [ms]`[i]<-round(mean(intermediate_adj$`Time to peak [ms]`, na.rm=TRUE),0)
  output$`Time to peak - standard dev [ms]`[i]<-round(sd(intermediate_adj$`Time to peak [ms]`, na.rm=TRUE),0)
  
  output$`Time to peak (5% to peak) [ms]`[i]<-round(mean(intermediate_adj$`Time to peak (5% to peak) [ms]`, na.rm=TRUE),0)
  output$`Time to peak (5% to peak) [ms] - standard dev [ms]`[i]<-round(sd(intermediate_adj$`Time to peak (5% to peak) [ms]`, na.rm=TRUE),0)
  
  output$`Relaxation time [ms]`[i]<-round(mean(intermediate_adj$`Relaxation time [ms]`, na.rm=TRUE),0)
  output$`Relaxation time - standard dev [ms]`[i]<-round(sd(intermediate_adj$`Relaxation time [ms]`, na.rm=TRUE),0)
  
  output$`Relaxation time (peak to 5%) [ms]`[i]<-round(mean(intermediate_adj$`Relaxation time (peak to 5%) [ms]`, na.rm=TRUE),0)
  output$`Relaxation time (peak to 5%) [ms] - standard dev [ms]`[i]<-round(sd(intermediate_adj$`Relaxation time (peak to 5%) [ms]`, na.rm=TRUE),0)
  
  output$`Contraction amplitude [mm]`[i]<-signif(mean(intermediate_adj$`Contraction amplitude [mm]`, na.rm=TRUE),3)
  output$`Contraction amplitude - standard dev [mm]`[i]<-signif(sd(intermediate_adj$`Contraction amplitude [mm]`, na.rm=TRUE),3)
  
  output$`Contraction force [µN]`[i]<-signif(mean(intermediate_adj$`Contraction force [µN]`, na.rm=TRUE),3)
  output$`Contraction force - standard dev [µN]`[i]<-signif(sd(intermediate_adj$`Contraction force [µN]`, na.rm=TRUE),3)
  
  plot(Frame/FrameRate,Trace_corrected$Displacement*pixelsize*DisplacementConversionFactor*1000,type="l",xlab="Time [s]", ylab="Contraction Force [µN]",main=strsplit(FileList[i], ".txt")[[1]][1])+abline(h = 0, col="red")+abline(v = peaks[1:(length(peaks[,1])),2]/FrameRate, col="blue")+abline(v = peaks[1:(length(peaks[,1])),3]/FrameRate, col="cyan")+abline(v = peaks[1:(length(peaks[,1])),4]/FrameRate, col="green")+abline(v = intermediate_adj$PeakPos..Frame./FrameRate, col="magenta",lwd =2)
}
# Save summary of peak-related information per Pillar
rownames(output)<-FileList
write.csv(output, paste(dir,"SummaryPerPillar.csv",sep="/"), row.names=TRUE)

# Save summary of peak-related information per EHT
if (pillarCount == 2) {
  outputSummary<-do.call(rbind,
                         lapply(seq(1, nrow(output), 2), function(i){
                           x <- output[ i:(i + 1), , drop = FALSE]
                           res <- rbind(colSums(x)/2)
                           rownames(res)[ nrow(res) ] <- strsplit(rownames(x), "_")[[1]][1]
                           res
                         }))
  outputSummary<-data.frame(outputSummary)
  colnames(outputSummary)<-c("Peak to Peak [s]","Contraction duration [ms]","Time to peak [ms]","Relaxation time [ms]","Contraction amplitude [mm]","Total contraction force [µN]")
  outputSummary$`Total contraction force [µN]`<-2*outputSummary$`Total contraction force [µN]`
  output$`Estimated total contraction force - standard dev [µN]`<-2*output$`Contraction force - standard dev [µN]`
  write.csv(outputSummary, paste(dir,"SummaryPerEHT.csv",sep="/"), row.names=TRUE)
} else if (pillarCount == 1){
  output$`Estimated total contraction force [µN]`<-2*output$`Contraction force [µN]`
  output$`Estimated total contraction force - standard dev [µN]`<-2*output$`Contraction force - standard dev [µN]`
  write.csv(output, paste(dir,"SummaryPerEHT(estimate).csv",sep="/"), row.names=TRUE)
}
rm(list=ls())