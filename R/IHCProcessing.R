#' @title IHCRProcessing: a series of functions for or ImageJ CD61 Megakaryocytes data processing.
#'
#' @description This packages is specifically developed for Linda Resar Lab to fully processing the CD61 IHC data or simliar
#' kind of cell recognition raw data output from ImageJ. You can input your .csv raw data, and use functions in
#' this package to reorganize, categorize, label and plot the data.
#'
#' @details
#' This package contains 5 functions:
#' \code{readraw} read the raw data and organize them into a list
#' \code{labind} bind and label each group of data
#' \code{numplot} plot number proportion
#' \code{areaplot} plot area proportion
#' \code{test} return p-value of multiple statistical test methoods
#' 1. load the raw data into R
#' 2. label and bind each groups and get the summary and percentage of mega and nuclei, respectively.
#' 3. plot (which can be processed in R or Prism 8)
#'
#' @author Ruitao Zhang
#'
#' @section Maintainer:
#' Ruitao Zhang <rzhangh49@jhu.edu>
#'
#' @docType package
#' @name IHCProcessing
NULL


#' Read the Megakaryocytes and nuclei raw data from the document
#'
#' This function is used to read the raw data read from the document and make them a list
#'
#' @param wd the working directory of the raw data
#' @return a list with all groups of raw data
#' @author Ruitao Zhang
#' @details
#' This function is used to read the raw data read from the document and make them a list. Each element of the list is a plot.
#' @export
#' @import Rmisc


readraw <- function(wd){
  setwd(wd)
  temp = list.files(pattern="*.csv")
  raw = lapply(temp, read.csv)
  return(raw)
}


#' Organizing the Megakaryocytes and nuclei together and label
#'
#' This function is used to transform the raw data read from last step into well organized and labeled data.
#' And also calculating the percentage of megakaryocytes number and area propotion for further analyze.
#'
#' @param nucleiraw the nuclei raw data
#' @param megaraw the megakaryocytes raw data
#' @return A list contains: 1. labeled and binded raw data; 2. labeled data with the megakaryocytes coount percentage of each group;
#' 3. labeled data with the megakaryocytes area percentage of each group;
#' @author Ruitao Zhang
#' @details
#' This function is used to transform the raw data read from last step into well organized and labeled data.
#' And also calculating the percentage of megakaryocytes number and area propotion for further analyze.
#' This function will return us a list contains 1. labeled and binded raw data; 2. labeled data with the megakaryocytes count percentage of each group;
#' 3. labeled data with the megakaryocytes area percentage of each group;
#' @export
#' @import forcats
#' @import dplyr
#' @importFrom  magrittr %>%
#' @import tidyr

labind <- function(nucleiraw, megaraw, Type = NULL){
  area <- data.frame()
  for(i in 1:length(nucleiraw)){
    n <- as.data.frame(nucleiraw[[i]][,2]) # get area data
    n$Label <- rep("nuclei", nrow(nucleiraw[[i]]))
    colnames(n) <- c('Area', 'Label')
    mega <- as.data.frame(megaraw[[i]][,2])
    mega$Label <- rep("Mega", nrow(megaraw[[i]]))
    colnames(mega) <- c('Area', 'Label')
    dat <- rbind(n, mega)
    dat$type <- rep(Type, nrow(dat))
    dat$group <- rep(i, nrow(dat))
    # area calculation
    Percent <- sum(dat$Area[dat$Label == "Mega"])/sum(dat$Area)
    area <- rbind(area, dat)
    # counting
    if(i == 1){
      count <- aggregate(data.frame(count = dat$Label), list(value = dat$Label), length)
      names(count)[2] <- i
    }
    else{
      c <- aggregate(data.frame(count = dat$Label), list(value = dat$Label), length)
      count <- cbind(count, c[,2])
      colnames(count)[i+1] <- i
    }
  }
  # get count clean data
  datl <- count[,-1]
  rownames(datl) <- count[,1]
  dat <- as.data.frame(t(datl))
  percentage <- dat$Mega/(dat$Mega+dat$nuclei)

  dat$percentage <- percentage
  dat$type <- rep(Type, nrow(dat))

  # get mean for wt
  mwt <- vector()
  for(k in 1: length(megaraw)){
    m <- mean(megaraw[[k]]$Area)
    mwt <- c(mwt, m)
  }

  # get sd for wt
  sdwt <- vector()
  for(k in 1: length(megaraw)){
    sd <- sd(megaraw[[k]]$Area)
    sdwt <- c(sdwt, sd)
  }

  # get area clean data
  area <- area %>% group_by(group) %>% mutate(percentage = (Area/sum(Area)*100))
  a <- aggregate(subset(area, Label == 'Mega')$percentage, list(value = subset(area, Label == 'Mega')$group), sum)
  a$type <- rep(Type, nrow(a))
  a$mean <- mwt
  a$sd <- sdwt
  colnames(a)[2] <- "percentage"

  # return
  return(list(area,dat,a))
}


#' Extracting and output labeled general data for downward analysis
#'
#' This function is used to generate the processed general labeled data read from last step,
#' which can also output the data file for further Prism analyze.
#'
#' @param data a \code{list} data output from the function \code{labind}
#' @param output if you need to output this data for analysis in other software
#' @return A data frame with labeled and binded general data.
#' @author Ruitao Zhang
#' @details
#' This function is used to extract and output labeled general data for downward analysis, generate
#' the processed general labeled data read from last step, which can also output the data file for
#' further Prism analyze. It will return a data frame with labeled and binded general data, and it's
#' your choice if output it or not.
#' @export

gendata <- function(data, output = FALSE, name = NULL){
  return(as.data.frame(data[[1]]))
  if(output == TRUE){write.csv(as.data.frame(data[[1]]), file = data)}
  else{NULL}
}


#' Extracting and output labeled count data for downward analysis
#'
#' This function is used to generate the processed labeled count data read from last step,
#' which can also output the data file for further Prism analyze.
#'
#' @param data a \code{list} data output from the function \code{labind}
#' @param output if you need to output this data for analysis in other software
#' @param name the name of your ooutput file
#' @return A data frame with labeled and binded general data.
#' @author Ruitao Zhang
#' @details
#' This function is used to extract and output labeled count data for downward analysis, generate
#' the processed general labeled data read from last step, which can also output the data file for
#' further Prism analyze. It will return a data frame with labeled and binded general data, and it's
#' your choice if output it or not.
#' @export

ctdata <- function(data, output = FALSE, name = NULL){
  return(as.data.frame(data[[2]]))
  if(output == TRUE){write.csv(as.data.frame(data[[2]]), file = name)}
  else{NULL}
}


#' Extracting and output labeled area data for downward analysis
#'
#' This function is used to generate the processed labeled area data read from last step,
#' which can also output the data file for further Prism analyze.
#'
#' @param data a \code{list} data output from the function \code{labind}
#' @param output if you need to output this data for analysis in other software
#' @return A data frame with labeled and binded area data.
#' @author Ruitao Zhang
#' @details
#' This function is used to extract and output labeled area data for downward analysis, generate
#' the processed general labeled data read from last step, which can also output the data file for
#' further Prism analyze. It will return a data frame with labeled and binded general data, and it's
#' your choice if output it or not.
#' @export

areadata <- function(data, output = FALSE){
  return(as.data.frame(data[[3]]))
  if(output == TRUE){write.csv(as.data.frame(data[[3]]), file = data)}
  else{NULL}
}


#' plotting the count percentage
#'
#' This function is used organizing and plotting the number prportion: how many megakaryocytes in
#' all nuclei in this plot
#'
#' @param ... the organized count data of each type of phenotype
#' @param type which kind of plot you want to output
#' @return a boxplot/violin plot of mega number proprtion
#' @author Ruitao Zhang
#' @details
#' This function is used organizing and plotting the number prportion: how many megakaryocytes in
#' all nuclei in this plot. This function will return a boxplot/violin plot of mega number proprtion
#' @export
#' @import ggplot2

numplot <- function(..., type = c("box", "violin")){
  # count percentage boxplot
  ct <- rbind(...)
  ctmeans <- aggregate(percentage ~  type, ct, mean)
  ctmeans$percentage <- round(ctmeans$percentage, 3)
  box <- ggplot(ct, aes(x = type, y = percentage*100, fill = type)) +
    theme_bw() +
    geom_boxplot(width = 0.5) + geom_jitter(width=0.1, alpha=0.5) +
    geom_text(data = ctmeans, aes(label = percentage), position = position_nudge(y = 0.3)) +
    labs(x = "", y = "Mega Percent(%)") +
    ggtitle("Megakaryocytes Proportion") +
    theme(plot.title = element_text(hjust = 0.5, size = 35, face = "bold"),
          axis.text.x = element_text(size = 20), axis.title.y = element_text(size = 20),
          legend.title=element_text(size=rel(2)), legend.text=element_text(size=rel(1.5)), legend.key.size = unit(1.5, "cm"))

  # count percentage violin plot
  violin <- ggplot(ct, aes(x = type, y = percentage*100, fill = type)) +
    theme_bw() +
    geom_violin() + geom_boxplot(width=0.03, fill = 'white') +
    geom_jitter(width=0.1,alpha=0.5) +
    geom_text(data = ctmeans, aes(label = percentage), position = position_nudge(y = 0.3)) +
    labs(x = "", y = "Mega Percent(%)") +
    ggtitle("Megakaryocytes Proportion") +
    theme(plot.title = element_text(hjust = 0.5, size = 35, face = "bold"),
          axis.text.x = element_text(size = 20), axis.title.y = element_text(size = 20),
          legend.title=element_text(size=rel(2)), legend.text=element_text(size=rel(1.5)), legend.key.size = unit(1.5, "cm"))

  if(type == "violin"){return(violin)}
  else(return(box))
}

#' plotting the area percentage
#'
#' This function is used organizing and plotting the area prportion: the area that megakaryocytes
#' occupied in all area occupied by cells
#'
#' @param ... the organized count data of each type of phenotype
#' @param type which kind of plot you want to output
#' @return a boxplot/violin plot/barplot of mega area proprtion
#' @author Ruitao Zhang
#' @details
#' This function is used organizing and plotting the area prportion: the area that megakaryocytes
#' occupied in all area occupied by cells. This function will return a boxplot/violin plot/barplot of mega area proprtion
#' @export
#' @import ggplot2

areaplot <- function(..., type = c("box", "violin", "bar")){
  Area <- rbind(...)
  Area$type <- factor(Area$type)


  ameans <- aggregate(percentage ~  type, Area, mean)
  ameans$x <- round(ameans$percentage, 3)
  box <- ggplot(Area, aes(x = type, y = percentage, fill = type)) +
    theme_bw() +
    geom_boxplot(width = 0.5) + geom_jitter(width=0.1,alpha=0.2) +
    geom_text(data = ameans, aes(label = percentage), position = position_nudge(y = 8))+
    labs(x = "", y = "Mega Percent(%)") +
    ggtitle("Megakaryocytes Area Occupied Proportion") +
    theme(plot.title = element_text(hjust = 0.5, size = 35, face = "bold"),
          axis.text.x = element_text(size = 20), axis.title.y = element_text(size = 20),
          legend.title=element_text(size=rel(2)), legend.text=element_text(size=rel(1.5)), legend.key.size = unit(1.5, "cm"))

  # area percentage violin plot
  violin <- ggplot(Area, aes(x = type, y = percentage, fill = type)) +
    theme_bw() +
    theme(panel.grid.major.x = element_blank())+
    #opts(axis.title.x = theme_blank(), axis.title.y = theme_blank())+
    geom_violin() + geom_boxplot(width=0.05, fill="white") +
    geom_jitter(width=0.1,alpha=0.5) +
    geom_text(data = ameans, aes(label = percentage), position = position_nudge(y = 8))+
    labs(x = "", y = "Mega Percent(%)") +
    ggtitle("Megakaryocytes Area Occupied Proportion") +
    theme(plot.title = element_text(hjust = 0.5, size = 35, face = "bold"),
          axis.text.x = element_text(size = 20), axis.title.y = element_text(size = 20),
          legend.title=element_text(size=rel(2)), legend.text=element_text(size=rel(1.5)), legend.key.size = unit(1.5, "cm"))

  # area percentage barplot
  ggplot(Area) +
    geom_bar(aes(x = type, y = percentage, fill = type), position="dodge", stat = "summary") +
    labs(x = "", y = "Mega Percent(%)") +
    ggtitle("Megakaryocytes Area Occupied Proportion") +
    theme(plot.title = element_text(hjust = 0.5, size = 35, face = "bold"),
          axis.text.x = element_text(size = 20), axis.title.y = element_text(size = 20),
          legend.title=element_text(size=rel(2)), legend.text=element_text(size=rel(1.5)), legend.key.size = unit(1.5, "cm"))

  if(type == "violin"){return(violin)}
  else if(type == "bar"){return(bar)}
  else(return(box))

}


#' Get the p-value of multiple statistical test methoods
#'
#' This function provides multiple statistical test methods for the unpaired-two sample testing.
#'
#' @param a the first data list you want to test
#' @param b the second data list you want to test
#' @param statistics the methods you use for statistical test
#' @return the p-value of the test result
#' @author Ruitao Zhang
#' @details
#' This function provides multiple statistical test methods for the unpaired-two sample testing. And it will give
#' us the p-value of the test result as an output.
#' @export

test <- function(a, b, statistics = c("Student","Mann","Chisq","Kolmogorov","Fisher")){
  if(statistics == "Chisq"){
  p <- chisq.test(a$percentage, b$percentage)$p.value
  cat("Chi Squared test: p-value =",p,sep = " ")
  }
  else if(statistics == "Mann"){
    p <- wilcox.test(a$percentage, b$percentage)$p.value
    cat("Mann-Whiethy U test: p-value =",p,sep = " ")
  }
  else if(statistics == "Kolmogorov"){
    p <- ks.test(a$percentage, b$percentage)$p.value
    cat("Kolmogorov & Smirnov test: p-value =",p,sep = " ")
  }
  else if(statistics == "Fisher"){
    p <- var.test(a$percentage, b$percentage)$p.value
    cat("Fisher’s F test: p-value =",p,sep = " ")
  }
  else{
    p <- t.test(a$percentage, b$percentage)$p.value
    cat("Student's-t test: p-value =",p,sep = " ")
  }
  }
