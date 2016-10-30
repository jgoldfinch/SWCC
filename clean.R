# Cleaning SWCC data 
# Jessie Golding
# 10/30/2016

###############################################
# Load Packages
install.packages("readr")
install.packages("tidyr")
install.packages("dplyr")
require(readr)
require(tidyr)
require(dplyr)



#######################################
# Source functions from GitHub
script_url <-"https://raw.githubusercontent.com/Huh/Rexamples/master/DataValidation/Check_char_num.R"
downloader::source_url(script_url, prompt=F, quiet=T)

# Create directories
data_dir <-"C:/Users/jgolding/Documents/USFS_R1_Carnivore_Monitoring/data/SWCC_do_not_distribute/data"
result_dir <-"C:/Users/jgolding/Documents/USFS_R1_Carnivore_Monitoring/data/SWCC_do_not_distribute/results"

#######################################
# Read in data
rawtd <-read_csv(#from readr function)
  file.path(
    data_dir,
    list.files(data_dir, pattern = "csv$")
  )
)

##############################################
# Subset columns to the ones you might use
use_col <-c("Year","Grid Cell","Survey Route ID","Date","Overall Tracking Conditions","Suspected Species","Confidence",
            "Species results","Carnivore tracks present?", "Genetics collected?")

keep <-rawtd %>%
  select(which(colnames(.)%in% use_col))
  
keep <-filter(keep, keep$`Genetics collected?` == "YES")

#remove the rows with no genetic ID
keep2 <-keep[!(is.na(keep$`Species results`) | keep$`Species results`== "POOR DNA"| keep$`Species results` == "no sample"),]

keep2$correct <-ifelse(keep2$`Suspected Species`==keep2$`Species results`,1,0)
keep2$correct[is.na(keep2$correct),] <- 0

pc<-sum(keep2$correct)/nrow(keep2)
