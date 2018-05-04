library(hsdar)
library(fda) # For Outlier detection
library(fda.usc) # For Outlier detection
library(prospectr) # For spectral binning
library(gdata) #To drop levels of a df
library(caret)
library(reshape2)
library(cowplot)
library(VSURF)
library(colourpicker)
library(tidyverse)


# 1 Load data and remove unwanted subsets and wavebands

data.original <- read.csv("data/All samples_converted.csv", check.names = FALSE)
#data.original <- data.original[apply(data.original[, -1], MARGIN = 1, function(x) all(x < 100)), ] 
#data.original <- subset(data.original, !(Type %in% c('OHWa'))) 
data.original$Type <- drop.levels(data.original$Type)
as.vector(unique(data.original$Type))
data.rmv.noise <- data.original[,match('400', names(data.original)):match('2400', names(data.original))] #this removes the unwanted bands
names(data.rmv.noise[,c(1,2001)]) #Check bands after removal (_superseded names(data.rmv.noise[,c(52,1452)]))

data.wo.noise <- cbind(data.original['Type'],data.rmv.noise)# Final df after noisy end are removed
 
subsets <- split(data.wo.noise, data.wo.noise$Type)#create a list of subsets based on the factor value.

subsets$OHWa <- subsets$OHWa[apply(subsets$OHWa[, -1], MARGIN = 1, function(x) all(x < 100)), ] 
subsets$ASP <- subsets$ASP[apply(subsets$ASP[, -1], MARGIN = 1, function(x) all(x < 100)), ] 
subsets$C <- subsets$C[apply(subsets$C[, -1], MARGIN = 1, function(x) all(x < 100)), ] 
subsets$OXa <- subsets$OXa[apply(subsets$OXa[, -1], MARGIN = 1, function(x) all(x < 100)), ] 
subsets$OXf <- subsets$OXf[apply(subsets$OXf[, -1], MARGIN = 1, function(x) all(x < 100)), ] 

# 2 Screen Spectra manually for outlier
unique(data.wo.noise$Type)
pdf("all.before.manual.2.pdf", width = 32, height =18)
par(mfrow=c(6,9))
for (i in 1:length(subsets)){
        
        labnames <- list(main="Width", xlab="Wavelength [nm]", ylab="Reflectance [%]")  
        set <- subsets[[i]] 
        set.2 <- as.matrix(set[,2:2002])
        f.res <- fdata(set.2, argvals = as.integer(names(set[,2:2002])), names = labnames)
        
       plot(f.res, main = names(subsets)[i])
       
}
dev.off()

# 3 Remove obvious outlier per hand according to plots

#HairyButtercupFlower

which.min(subsets$HairyButtercupFlowerYellow[,'1000'] )
subsets$HairyButtercupFlowerYellow <- subsets$HairyButtercupFlowerYellow[-7,] # use multiple times and check plot if spectra are gone

#Remove specific subsets

pdf("all.after.manual.pdf", width = 32, height =18)
par(mfrow=c(4,7))
for (i in 1:length(subsets)){
        
        labnames <- list(main="Width", xlab="Wavelength [nm]", ylab="Reflectance [%]")  
        set <- subsets[[i]] 
        set.2 <- as.matrix(set[,2:2002])
        f.res <- fdata(set.2, argvals = as.integer(names(set[,2:2002])), names = labnames)
        
        plot(f.res, main = names(subsets)[i])
        
}
dev.off()


# 4 fda pkg outlier detection

source('R/Remove_Functional_Outlier_June2017.R')


cleaned_data <- lapply(subsets, rmv.funct.outlier)

pdf("all.after.automatic.pdf", width = 32, height =18)
par(mfrow=c(4,7))
for (i in 1:length(cleaned_data)){
        
        labnames <- list(main="Width", xlab="Wavelength [nm]", ylab="Reflectance [%]")  
        set <- cleaned_data[[i]] 
        set.2 <- as.matrix(set[,2:902])
        f.res <- fdata(set.2, argvals = as.integer(names(set[,2:902])), names = labnames)
        
        plot(f.res, main = names(cleaned_data)[i])
        
}
dev.off()

cleaned.df <- bind_rows(cleaned_data)
