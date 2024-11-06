install.packages(c("camtrapR", "tidyr","raster","secr"))
library(camtrapR)
library(raster)
library(secr)
library(tidyr)

#exiftool is required for reading data from camera trap images 
# for more info read https://jniedballa.github.io/camtrapR/ 

addToPath( "C:/Program Files/exiftool-12.96_64/exiftool-12.96_64")
Sys.which("exiftool.exe")



########density estimation for 2022

## Read and plot the shapefile 


#Reading the shapefile for the wildlife sanctuary
NWLS<-shapefile('nwls_utm44N.shp')
plot(NWLS)



## Reading the Camera operation file to create a matrix
NWLS_camop_2022<-read.csv("2022/gps_2022_utm44.csv")
head(NWLS_camop_2022)
cameraOperation

#creating on-off matrix/trap matrix for camera trap stations to understand effort 
NWLS_camop_table_2022<-cameraOperation(NWLS_camop_2022,
                                  stationCol = "station",
                                  setupCol = "Deployment",
                                  retrievalCol = "Retrieving",
                                  dateFormat = "%d-%m-%Y",
                                  writecsv= TRUE,
                                  outDir= "2022/Output_2022")

head(NWLS_camop_table_2022)

# this command is used for creating a capture matrix for the tigers
# It is not used here as directory 'experiment_2022' which has tiger images cannot 
# be shared without permission from Wildlife Institute of India
tiger_record_table_2022<-recordTableIndividual(inDir = "2022/experiment_2022",
                                           hasStationFolders = TRUE,           
                                           IDfrom = "directory",
                                           camerasIndependent = TRUE,
                                           minDeltaTime = 01,
                                           deltaTimeComparedTo = "lastIndependentRecord",
                                           timeZone = "Asia/Calcutta",
                                           stationCol = "station",
                                           writecsv = TRUE,
                                           outDir = "2022/Output_2022")

# output for the capture matrix is provided  to allow analysis using secr package
tiger_record_table_2022<-read.csv("record_table_individuals1min_deltaT_2022-11-06.csv")
#

# creating a capthist object required for secr
# need trap matrix and capture matrix
tiger_capthist_2022 <-spatialDetectionHistory(tiger_record_table_2022,
                                          species = "experiment_2022",
                                          output = "binary",
                                          camOp = NWLS_camop_table_2022,
                                          CTtable = NWLS_camop_2022,
                                          stationCol = "station",
                                          Xcol = "X",
                                          Ycol= "Y",
                                          individualCol = "Individual",
                                          recordDateTimeCol     = "DateTimeOriginal",
                                          recordDateTimeFormat  = "%Y-%m-%d %H:%M",
                                          occasionLength        = 1,
                                          day1 = "station",
                                          includeEffort         = TRUE,
                                          timeZone              = "Asia/Calcutta"
)
write.csv(tiger_capthist_2022, "2022/Output_2022/tiger_capthist_2022.csv")

summary(tiger_capthist_2022)
plot(tiger_capthist_2022, tracks = TRUE)

# command to understand how much buffer is required 
suggest.buffer(tiger_capthist_2022)
hist(unlist(moves(tiger_capthist_2022)))

#creating a habitat Mask
nwls_mask_2022<-make.mask(traps(tiger_capthist_2022), buffer = 9000, spacing = 800, type = "trapbuffer", poly = NWLS, poly.habitat = TRUE)
plot(nwls_mask_2022)

####using different models for density estimation using secr package
##Null Model-
tiger_secr_2022<-secr.fit(tiger_capthist_2022, mask = nwls_mask_2022)
##detection probablity + heterogenity-
tiger_secr_go_het_2022<-secr.fit(tiger_capthist_2022, model = g0~h2, mask = nwls_mask_2022)
##movement parameters + heterogenity-
tiger_secr_sig_het_2022<-secr.fit(tiger_capthist_2022, model = sigma~h2, mask = nwls_mask_2022)
##Combination
tiger_secr_go_sig_het_2022<-secr.fit(tiger_capthist_2022, model = list(g0~h2, sigma~h2), mask = nwls_mask_2022)


#checking the AIC values for the all the models used
AIC(tiger_secr_2022, tiger_secr_go_het_2022, tiger_secr_sig_het_2022, tiger_secr_go_sig_het_2022)



#summary for best fit model
summary(tiger_secr_go_sig_het_2022)

# realised and calculated estimates of tiger number 
region.N(tiger_secr_sig_het_2022)








#################################density estimation for 2023 


## Reading the Camera operation file to create a matrix
NWLS_camop_2023 <- read.csv("2023/gps_2023_utm44.csv")
head(NWLS_camop_2023)
cameraOperation

# creating on-off matrix/trap matrix for camera trap stations to understand effort
NWLS_camop_table_2023 <- cameraOperation(NWLS_camop_2023,
                                         stationCol = "station",
                                         setupCol = "Deployment",
                                         retrievalCol = "Retrieving",
                                         dateFormat = "%d-%m-%Y",
                                         writecsv = TRUE,
                                         outDir = "2023/Output_2023")

head(NWLS_camop_table_2023)

# this command is used for creating a capture matrix for the tigers
# It is not used here as directory 'experiment_2023' which has tiger images cannot
# be shared without permission from Wildlife Institute of India
tiger_record_table_2023 <- recordTableIndividual(inDir = "2023/experiment_2023",
                                                 hasStationFolders = TRUE,
                                                 IDfrom = "directory",
                                                 camerasIndependent = TRUE,
                                                 minDeltaTime = 01,
                                                 deltaTimeComparedTo = "lastIndependentRecord",
                                                 timeZone = "Asia/Calcutta",
                                                 stationCol = "station",
                                                 writecsv = TRUE,
                                                 outDir = "2023/Output_2023")

# output for the capture matrix is provided to allow analysis using secr package
tiger_record_table_2023 <- read.csv("record_table_individuals1min_deltaT_2023-11-06.csv")

# creating a capthist object required for secr
# need trap matrix and capture matrix
tiger_capthist_2023 <- spatialDetectionHistory(tiger_record_table_2023,
                                               species = "experiment_2023",
                                               output = "binary",
                                               camOp = NWLS_camop_table_2023,
                                               CTtable = NWLS_camop_2023,
                                               stationCol = "station",
                                               Xcol = "X",
                                               Ycol = "Y",
                                               individualCol = "Individual",
                                               recordDateTimeCol = "DateTimeOriginal",
                                               recordDateTimeFormat = "%Y-%m-%d %H:%M",
                                               occasionLength = 1,
                                               day1 = "station",
                                               includeEffort = TRUE,
                                               timeZone = "Asia/Calcutta"
)
write.csv(tiger_capthist_2023, "2023/Output_2023/tiger_capthist_2023.csv")

summary(tiger_capthist_2023)
plot(tiger_capthist_2023, tracks = TRUE)

# command to understand how much buffer is required
suggest.buffer(tiger_capthist_2023)
hist(unlist(moves(tiger_capthist_2023)))

# creating a habitat Mask
nwls_mask_2023 <- make.mask(traps(tiger_capthist_2023), buffer = 10000, spacing = 950, type = "trapbuffer", poly = NWLS, poly.habitat = TRUE)
plot(nwls_mask_2023)

#### using different models for density estimation using secr package
## Null Model-
tiger_secr_2023 <- secr.fit(tiger_capthist_2023, mask = nwls_mask_2023)
## detection probability + heterogeneity-
tiger_secr_go_het_2023 <- secr.fit(tiger_capthist_2023, model = g0 ~ h2, mask = nwls_mask_2023)
## movement parameters + heterogeneity-
tiger_secr_sig_het_2023 <- secr.fit(tiger_capthist_2023, model = sigma ~ h2, mask = nwls_mask_2023)
## Combination
tiger_secr_go_sig_het_2023 <- secr.fit(tiger_capthist_2023, model = list(g0 ~ h2, sigma ~ h2), mask = nwls_mask_2023)

# checking the AIC values for the all the models used
AIC(tiger_secr_2023, tiger_secr_go_het_2023, tiger_secr_sig_het_2023, tiger_secr_go_sig_het_2023)

# summary for best fit model
summary(tiger_secr_go_sig_het_2023)

# realised and calculated estimates of tiger number
region.N(tiger_secr_go_sig_het_2023)














########################################density estimation for 2024 

## Reading the Camera operation file to create a matrix
NWLS_camop_2024 <- read.csv("2024/gps_2024_utm44.csv")
head(NWLS_camop_2024)
cameraOperation

# creating on-off matrix/trap matrix for camera trap stations to understand effort
NWLS_camop_table_2024 <- cameraOperation(NWLS_camop_2024,
                                         stationCol = "station",
                                         setupCol = "Deployment",
                                         retrievalCol = "Retrieving",
                                         dateFormat = "%d-%m-%Y",
                                         writecsv = TRUE,
                                         outDir = "2024/Output_2024")

head(NWLS_camop_table_2024)

# this command is used for creating a capture matrix for the tigers
# It is not used here as directory 'experiment_2024' which has tiger images cannot
# be shared without permission from Wildlife Institute of India
tiger_record_table_2024 <- recordTableIndividual(inDir = "2024/experiment_2024",
                                                 hasStationFolders = TRUE,
                                                 IDfrom = "directory",
                                                 camerasIndependent = TRUE,
                                                 minDeltaTime = 01,
                                                 deltaTimeComparedTo = "lastIndependentRecord",
                                                 timeZone = "Asia/Calcutta",
                                                 stationCol = "station",
                                                 writecsv = TRUE,
                                                 outDir = "2024/Output_2024")

# output for the capture matrix is provided to allow analysis using secr package
tiger_record_table_2024 <- read.csv("record_table_individuals1min_deltaT_2024-11-06.csv")

# creating a capthist object required for secr
# need trap matrix and capture matrix
tiger_capthist_2024 <- spatialDetectionHistory(tiger_record_table_2024,
                                               species = "experiment_2024",
                                               output = "binary",
                                               camOp = NWLS_camop_table_2024,
                                               CTtable = NWLS_camop_2024,
                                               stationCol = "station",
                                               Xcol = "X",
                                               Ycol = "Y",
                                               individualCol = "Individual",
                                               recordDateTimeCol = "DateTimeOriginal",
                                               recordDateTimeFormat = "%Y-%m-%d %H:%M",
                                               occasionLength = 1,
                                               day1 = "station",
                                               includeEffort = TRUE,
                                               timeZone = "Asia/Calcutta"
)
write.csv(tiger_capthist_2024, "2024/Output_2024/tiger_capthist_2024.csv")

summary(tiger_capthist_2024)
plot(tiger_capthist_2024, tracks = TRUE)

# command to understand how much buffer is required
suggest.buffer(tiger_capthist_2024)
hist(unlist(moves(tiger_capthist_2024)))

# creating a habitat Mask
nwls_mask_2024 <- make.mask(traps(tiger_capthist_2024), buffer = 10000, spacing = 1000, type = "trapbuffer", poly = NWLS, poly.habitat = TRUE)
plot(nwls_mask_2024)

#### using different models for density estimation using secr package
## Null Model-
tiger_secr_2024 <- secr.fit(tiger_capthist_2024, mask = nwls_mask_2024)
## detection probability + heterogeneity-
tiger_secr_go_het_2024 <- secr.fit(tiger_capthist_2024, model = g0 ~ h2, mask = nwls_mask_2024)
## movement parameters + heterogeneity-
tiger_secr_sig_het_2024 <- secr.fit(tiger_capthist_2024, model = sigma ~ h2, mask = nwls_mask_2024)
## Combination
tiger_secr_go_sig_het_2024 <- secr.fit(tiger_capthist_2024, model = list(g0 ~ h2, sigma ~ h2), mask = nwls_mask_2024)

# checking the AIC values for the all the models used
AIC(tiger_secr_2024, tiger_secr_go_het_2024, tiger_secr_sig_het_2024, tiger_secr_go_sig_het_2024)

# summary for best fit model
summary(tiger_secr_sig_het_2024)

# realised and calculated estimates of tiger number
region.N(tiger_secr_sig_het_2024)

