# Script as an example and test for ImageJ CD61 Megakaryocytes data processing and IHCResarLab package
# 1. load the raw data into R
# 2. label and bind each groups and get the summary and percentage of mega and nuclei, respectively.
# 3. plot (which can be processed in R or Prism 8)

install.packages(c("ggplot2",
                 "Rmisc",
                 "forcats",
                 "dplyr",
                 "magrittr",
                 "tidyr")) # install all packages


library(ggplot2) # package for plotting
library(Rmisc) # package for read all sheets in one list data
library(forcats)
library(dplyr)
library(tidyr)
library(IHCResarLab) #Our main package for data processing

#1. Raw data loading

# load and list WT nuclei raw data in order
wildtype_nuclei_raw_data <- readraw("~/Desktop/Resar Lab/imageJ processing/IHC CD61 image processing kit/testdata/wt/nucleus/")
# load and list WT megakaryocytes raw data in order
wildtype_mega_raw_data <- readraw("~/Desktop/Resar Lab/imageJ processing/IHC CD61 image processing kit/testdata/wt/mega/")

# load and list VF nuclei raw data in order
vf_nuclei_raw_data <- readraw("~/Desktop/Resar Lab/imageJ processing/IHC CD61 image processing kit/testdata/vf/nucleus/")
# load and list VF megakaryocytes raw data in order
vf_mega_raw_data <- readraw("~/Desktop/Resar Lab/imageJ processing/IHC CD61 image processing kit/testdata/vf/mega/")

# load and list TM1B nuclei raw data in order
tm1b_nuclei_raw_data <- readraw("~/Desktop/Resar Lab/imageJ processing/IHC CD61 image processing kit/testdata/tm1b/nucleus/")
# load and list TM1B megakaryocytes raw data in order
tm1b_mega_raw_data <- readraw("~/Desktop/Resar Lab/imageJ processing/IHC CD61 image processing kit/testdata/tm1b/mega/")

# load and list TPOKO nuclei raw data in order
tpoko_nuclei_raw_data <- readraw("~/Desktop/Resar Lab/imageJ processing/IHC CD61 image processing kit/testdata/nucleus/")
# load and list TPOKO megakaryocytes raw data in order
tpoko_mega_raw_data <- readraw("~/Desktop/Resar Lab/imageJ processing/IHC CD61 image processing kit/testdata/mega/")


#2. labelling and binding the raw data
wildtype <- labind(wildtype_nuclei_raw_data, wildtype_mega_raw_data, Type = "WT")
jak2vf <- labind(vf_nuclei_raw_data, vf_mega_raw_data, Type = "VF")
vftm1b <- labind(tm1b_nuclei_raw_data, tm1b_mega_raw_data,Type = "VF/TM1B")
vftpoko <- labind(tpoko_nuclei_raw_data, tpoko_mega_raw_data,Type = "VF/TPO KO")


#3. extracting the data we need for plotting and testing
wt <- gendata(wildtype)
wt_count <- ctdata(wildtype)
wt_area <- areadata(wildtype)

vf <- gendata(jak2vf)
vf_count <- ctdata(jak2vf)
vf_area <- areadata(jak2vf)

tm1b <- gendata(vftm1b)
tm1b_count <- ctdata(vftm1b)
tm1b_area <- areadata(vftm1b)

tpoko <- gendata(vftpoko)
tpoko_count <- ctdata(vftpoko)
tpoko_area <- areadata(vftpoko)

#4. plotting

#number proportion
numplot(wt_count, vf_count, tm1b_count, tpoko_count, type = "box") # boxplot
numplot(wt_count, vf_count, tm1b_count, tpoko_count, type = "violin") # violin plot

#area proportion
areaplot(wt_area, vf_area, tm1b_area, tpoko_area, type = "box") # boxplot
areaplot(wt_area, vf_area, tm1b_area, tpoko_area, type = "bar") # barplot
areaplot(wt_area, vf_area, tm1b_area, tpoko_area, type = "violin") # violin plot

#5. statistical
test(wt_area, vf_area, statistics = "Student") # get p-value of student's t-test
test(vf_area, tm1b_area, statistics = "Mann") # get p-value of Mann-Whiethy U test
test(vf_area, tm1b_area, statistics = "Chisq") # get p-value of Chi Squared test
test(vf_area, tm1b_area, statistics = "Kolmogorov") # get p-value of Kolmogorov and Smirnov test
test(vf_area, tm1b_area, statistics = "Fisher") # get p-value of Fisher’s F test

