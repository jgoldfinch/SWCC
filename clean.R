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
# Read in track data
rawtd <-read_csv(#from readr function)
  file.path(
    data_dir,
    list.files(data_dir, pattern = "csv$")
  )
)

##############################################
# Subset columns to the ones you might use
use_col <-c("Year","Grid Cell","Survey Route ID","Date","Overall Tracking Conditions","Suspected Species","Confidence",
            "Species results","Carnivore tracks present?", "Genetics collected?", "Individual ID")

keep <-rawtd %>%
  select(which(colnames(.)%in% use_col))
  
keep <-filter(keep, keep$`Genetics collected?` == "YES")

# remove the rows with no genetic ID
keep2 <-keep[!(is.na(keep$`Species results`) | keep$`Species results`== "POOR DNA"| keep$`Species results` == "no sample"),]

# Assign a 1 if species was correctly identified, 0 if it was not
keep2$correct <-ifelse(keep2$`Suspected Species`==keep2$`Species results`,1,0)
# A single row has an NA for "Suspeceted species" column - assign a 0 for correct id 
keep2$correct[is.na(keep2$correct)] <- 0

# Calculate percent correct identified
pc<-sum(keep2$correct)/nrow(keep2)

# Rough calculation of successful effort for track surveys
nrow(keep2)/nrow(rawtd)
# Approximately 11% were successful


# Subset to look at individuals that were misidentified
keep3 <-keep2[keep2$correct==0,]

# Group by suspected species 
keep3 %>% group_by(`Suspected Species`) %>% summarise(length(unique(`Species results`)))



# Read in bait station data
rawbd <-read.csv("C:/Users/jgolding/Documents/USFS_R1_Carnivore_Monitoring/data/SWCC_do_not_distribute/data/SWCC_bait.csv")

# Subset columns to the ones you might use
use_col2 <-c("Year","Grid.Cell","Station.ID","Station.Easting","Station.Northing",
            "Set.up.Date", "Revisit.Date","Genetics.collected.","Sample.ID", 
            "Sample.Source", "Suspected.Species","Species","Individual.ID")

# Subset the raw data to keep only the columns you're interested in
keepb <-rawbd %>%
  select(which(colnames(.)%in% use_col2))

# Create a list of species to keep
keepsp<-c("BOBCAT","WOLVERINE","MARTEN","LION","LYNX","RED SQUIRREL","SKUNK","ERMINE","N FLYING SQUIRREL",
          "LONG-TAILED WEASEL","RED FOX", "BEAVER", "mink")

keepb2 <-filter(keepb, keepb$Species == "BOBCAT"|keepb$Species =="WOLVERINE"|)
