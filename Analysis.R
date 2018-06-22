# (1) how clean (spectral variation) these profiles are for each species
# (2) then what the accuracy of classification between these species. 
# (3) simulation of the reflectance profile of these species as measured by using multispectral sensors, and to test how distinguishable between the different species again. 
# - Sequoia (drone sensor), landsat-8 and Sentinel-2

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
library(tictoc)

source('R/20170601_FUN_DropCatVar.R')
source('R/20171224_FUN_raw2speclibhsdar.R')
source('R/FUN_Remove_Functional_Outlier.R')
source('R/FUN_prepggwide2long.R')

dir.create('output', FALSE, FALSE)


# 1 Load data and remove unwanted subsets and wavebands

data.original <- read.csv("data/All samples_converted.csv", check.names = FALSE)

data.rmv.noise <- data.original[,match('400', 
    names(data.original)):match('2400', 
    names(data.original))] #Rmv bands above and below given numeric values

names(data.rmv.noise[,c(1,ncol(data.rmv.noise))]) #Check bands after removal

data.wo.noise <- cbind(data.original['Type'],
                       data.rmv.noise) #Build final df
# 1.1 Tidy data

unique(data.wo.noise$Type)
'%!in%' <- function(x,y)!('%in%'(x,y)) #opposite of %in% operator

test <- data.wo.noise[data.wo.noise$Type %!in% c('AlpineSunray_Flower',
                                                 'Xsub_flower'), ]
test$Type <- drop.levels(test$Type)

unique(test$Type)

# 1.2 Append tidied versions by Michael Chang

a <- read.csv("data/AlpineSunray_Flower2.csv", check.names = FALSE)
b <- read.csv("data/XSub_flower2.csv", check.names = FALSE)

Type <- rep("AlpineSunray_Flower",length(a$Wavelength))
Type2 <- rep("XSub_flower",length(b$Wavelength))

a[,1] <- Type
a <- rename(a, Type=Wavelength)
a <- a[,match('400', names(a)):match('2400', names(a))]
a <- cbind(Type, a)
names(a)

b[,1] <- Type2
b <- rename(b, Type=Wavelength)
b <- b[,match('400', names(b)):match('2400', names(b))]
b <- cbind('Type'=Type2, b)
names(b)

df <- rbind(test, a, b)

unique(df$Type)

df <- df[order(df$Type),]
 
subsets <- split(data.wo.noise, data.wo.noise$Type) 

plots.bef <- list()
for (i in 1:length(subsets)){
  
  p <- prep_gg(as.data.frame(subsets[[i]]))
  
  plots.bef[[i]] <- ggplot(p, aes(Wavelength, Reflectance, colour = Type)) +
    geom_point(aes(shape=Type), size = .1)+
    labs(title=paste(names(subsets[i])), x= "WL", y="R")+
    theme_minimal()+
    theme(text=element_text(size=8))+
    theme(legend.position="none")
  
}

p.bef <- plot_grid(plotlist=plots.bef)

ggsave("output/before.pdf",
       plot = p.bef,
       width = 40,
       height = 20,
       units = "cm",
       dpi = 50
)


# 4 fda pkg outlier detection


cleaned_data <- lapply(subsets, rmv.funct.outlier,depth.mode, 1, 0.05, 0.5)


plots <- list()
for (i in 1:length(cleaned_data)){
  
  p <- prep_gg(as.data.frame(cleaned_data[[i]]), agg = FALSE)
  
  plots[[i]] <- ggplot(p, aes(Wavelength, Reflectance, colour = Type)) +
    geom_point(aes(shape=Type), size = .1)+
    labs(title=paste(names(cleaned_data[i])), x= "WL", y="R")+
    theme_minimal()+
    theme(text=element_text(size=8))+
    theme(legend.position="none")
    
        
}


p.aft <- plot_grid(plotlist=plots)

ggsave("output/after.pdf",
       plot = p.aft,
       width = 40,
       height = 20,
       units = "cm",
       dpi = 50
)


cleaned.df <- bind_rows(cleaned_data)

write.csv(cleaned.df, '20180505_clean_spectra.csv')

# Random Forest Classification -------------------------------------------------

classif <- read.csv('20180505_clean_spectra.csv', check.names = FALSE)
names(classif)
classif <- classif[,-1]
# set.seed(20180427)

# Partition Data

inTraining <- createDataPartition(classif$Type, p = .75, list = FALSE)
train <- classif[ inTraining,]
test  <- classif[-inTraining,]

# Tune random forest resample process to create variable samples for each tree

rfControl <- trainControl(
  method = "boot",
  number = 100
)

# Approximate mtry (# of vars randomly sampled as candidates at each split)

rfGrid <- expand.grid(mtry = trunc(sqrt(ncol(train)))) 

# RF Model Training

tic("RF") #Start timing

rfFit <- train(Type ~ ., data = train,
               method = "rf",
               importance = TRUE, ntree=500,
               trControl = rfControl, tuneGrid = rfGrid,
               metric = "Accuracy", maximize = TRUE) #runs ~2.2h

rfFit #Model output based on training data


# RF Model Testing

rfPred <- 
  predict.train(rfFit, test[, !names(test) %in% c("Type")], type = "raw")

toc() # Stop timing

# Collect results in a list and save

res.all <- 
  list(fit = rfFit,
       pred = predict.train(rfFit, test[, !names(test) %in% c("Type")], type = "raw"),
       confusion = confusionMatrix(rfPred, test$Type),
       varImp = varImp(rfFit, scale = FALSE))

sink(file = 'output/I_res_all.txt')
res.all
sink()

saveRDS(res.all, 'output/I_res_all.rds')
res.all <- readRDS("output/I_res_all.rds")

#Classification of leaves only
unique(classif$Type)
levels(classif$Type)

leav <- filter(classif, Type %in% c('AlpineGroundsel_leaf',
                            'AlpineSunray_leaf',
                            'BillyButton_leaf',
                            'BillyButton_leaf_mainrange',
                            'BroadLeaveGrass',
                            'ButtercupFelted_leaf',
                            'CandleHeath_strez_leaf',
                            'Carex_spp_leaf',
                            'CoralHeath_leaf',
                            'Dandelion_leaf',
                            'Flatweed',
                            'HairyButtercupLeaf',
                            'MixedGrass',
                            'MountainCelery_leaf',
                            'MouseearHW_leaf',
                            'MouseEarHW_strez_leaf',
                            'NativeYamDaisy_leaf',
                            'OHW_leaf',
                            'Pale_Everlasting_Strez_leaf',
                            'PineappleGrass_strez_leaf',
                            'PricklySnowGrass',
                            'SheepSorrell_leaf',
                            'SilverSnowDaisy_topsideLeaf',
                            'SilverSnowDaisy_undersideLeaf',
                            'SmallStarPlantain_leaf',
                            'Snow_Gentian_leaf',
                            'SpoonDaisy_strez_leaf',
                            'StJohn_Worf_leaf',
                            'Variable_Eyebright_leaf',
                            'WhiteCover_leaf',
                            'WoollyBillyButton_strez_leaf',
                            'Xsub_leaf'))
leav$Type <- drop.levels(leav$Type)

# Partition Data

inTraining <- createDataPartition(leav$Type, p = .75, list = FALSE)
train <- leav[ inTraining,]
test  <- leav[-inTraining,]

# Tune random forest resample process to create variable samples for each tree

rfControl <- trainControl(
  method = "boot",
  number = 100
)

# Approximate mtry (# of vars randomly sampled as candidates at each split)

rfGrid <- expand.grid(mtry = trunc(sqrt(ncol(train)))) 

# RF Model Training

tic("RF") #Start timing

rfFit <- train(Type ~ ., data = train,
               method = "rf",
               importance = TRUE, ntree=500,
               trControl = rfControl, tuneGrid = rfGrid,
               metric = "Accuracy", maximize = TRUE) #runs ~1.8h

rfFit #Model output based on training data


# RF Model Testing

rfPred <- 
  predict.train(rfFit, test[, !names(test) %in% c("Type")], type = "raw")

toc() # Stop timing

# Collect results in a list and save

res.leav <- 
  list(fit = rfFit,
       pred = predict.train(rfFit, test[, !names(test) %in% c("Type")], type = "raw"),
       confusion = confusionMatrix(rfPred, test$Type),
       varImp = varImp(rfFit, scale = FALSE))

sink(file = 'output/II_res_leav.txt')
res.leav
sink()

saveRDS(res.leav, 'output/II_res_leav.rds')
res.leav <- readRDS("output/II_res_leav.rds")
#What is CarpetHeath_L_F

#Classification of flowers only

flow <- classif[classif$Type %!in% c('AlpineGroundsel_leaf',
                                    'AlpineSunray_leaf',
                                    'BillyButton_leaf',
                                    'BillyButton_leaf_mainrange',
                                    'BroadLeaveGrass',
                                    'ButtercupFelted_leaf',
                                    'CandleHeath_strez_leaf',
                                    'Carex_spp_leaf',
                                    'CoralHeath_leaf',
                                    'Dandelion_leaf',
                                    'Flatweed',
                                    'HairyButtercupLeaf',
                                    'MixedGrass',
                                    'MountainCelery_leaf',
                                    'MouseearHW_leaf',
                                    'MouseEarHW_strez_leaf',
                                    'NativeYamDaisy_leaf',
                                    'OHW_leaf',
                                    'Pale_Everlasting_Strez_leaf',
                                    'PineappleGrass_strez_leaf',
                                    'PricklySnowGrass',
                                    'SheepSorrell_leaf',
                                    'SilverSnowDaisy_topsideLeaf',
                                    'SilverSnowDaisy_undersideLeaf',
                                    'SmallStarPlantain_leaf',
                                    'Snow_Gentian_leaf',
                                    'SpoonDaisy_strez_leaf',
                                    'StJohn_Worf_leaf',
                                    'Variable_Eyebright_leaf',
                                    'WhiteCover_leaf',
                                    'WoollyBillyButton_strez_leaf',
                                    'Xsub_leaf'),]

flow$Type <- drop.levels(flow$Type)

# Partition Data

inTraining <- createDataPartition(flow$Type, p = .75, list = FALSE)
train <- flow[ inTraining,]
test  <- flow[-inTraining,]

# Tune random forest resample process to create variable samples for each tree

rfControl <- trainControl(
  method = "boot",
  number = 100
)

# Approximate mtry (# of vars randomly sampled as candidates at each split)

rfGrid <- expand.grid(mtry = trunc(sqrt(ncol(train)))) 

# RF Model Training

tic("RF") #Start timing

rfFit <- train(Type ~ ., data = train,
               method = "rf",
               importance = TRUE, ntree=500,
               trControl = rfControl, tuneGrid = rfGrid,
               metric = "Accuracy", maximize = TRUE) #runs ~0.14h

rfFit #Model output based on training data


# RF Model Testing

rfPred <- 
  predict.train(rfFit, test[, !names(test) %in% c("Type")], type = "raw")

toc() # Stop timing

# Collect results in a list and save

res.flow <- 
  list(fit = rfFit,
       pred = predict.train(rfFit, test[, !names(test) %in% c("Type")], type = "raw"),
       confusion = confusionMatrix(rfPred, test$Type),
       varImp = varImp(rfFit, scale = FALSE))

sink(file = 'output/III_res_flow.txt')
res.flow
sink()

saveRDS(res.flow, 'output/III_res_flow.rds')
res.flow <- readRDS("output/III_res_flow.rds")

#Resample and classify

hypdata <- classif
#hypdata <- DropClass(hypdata, hypdata$Type, "Healthy")

# Create spectral library to use hsdar pkg

speclib <- raw2speclib(hypdata) #Function needs just numbers as colnames.

# Resample hyperspectral data to Micasense band specifications
#Green 530-570nm, Red 640-680nm, Red Edge 730-740 nm, NIR 770-810nm (Sequoia)

center <-  c(550, 660, 735, 790)
fwhm <- c(20, 20, 5, 20)

sequoia <- as.data.frame(cbind(center, fwhm))

data_seq <- spectralResampling(speclib, sequoia)
data_senti2 <- spectralResampling(speclib, 'Sentinel2')

plot(data_senti2)
plot(speclib)
plot(data_seq)

# Extract reflectance data from seq and senti spectral library for classification

sentidata <- as.data.frame(data_senti2@spectra@spectra_ma)
sentidata <- cbind('Type'=hypdata$Type, sentidata)
str(sentidata)
names(sentidata)

seqdata <- as.data.frame(data_seq@spectra@spectra_ma)
seqdata <- cbind('Type'=hypdata$Type, seqdata)
str(seqdata)
names(seqdata)
# Rename columns
newnamesSenti <- c("Type", "B1Aero", "B2Blue", "B3Green", "B4Red",
                   "B5RE1", "B6RE2", "B7RE3", "B8NIR", "B9WaterVap",
                   "B10SWIR1", "B11SWIR2", "B12SWIR3", "B13")

newnamesSeq <- c("Type", "Green", "Red", "RedEdge", "NIR")

names(sentidata) <- newnamesSenti
names(seqdata) <- newnamesSeq

# Sentinel2 Classification

# Partition Data

inTraining <- createDataPartition(sentidata$Type, p = .75, list = FALSE)
train <- sentidata[ inTraining,]
test  <- sentidata[-inTraining,]

# Tune random forest resample process to create variable samples for each tree

rfControl <- trainControl(
  method = "boot",
  number = 100
)

# Approximate mtry (# of vars randomly sampled as candidates at each split)

rfGrid <- expand.grid(mtry = trunc(sqrt(ncol(train)))) 

# RF Model Training

tic("RF") #Start timing

rfFit <- train(Type ~ ., data = train,
               method = "rf",
               importance = TRUE, ntree=500,
               trControl = rfControl, tuneGrid = rfGrid,
               metric = "Accuracy", maximize = TRUE) #runs ~0.14h

rfFit #Model output based on training data


# RF Model Testing

rfPred <- 
  predict.train(rfFit, test[, !names(test) %in% c("Type")], type = "raw")

toc() # Stop timing

# Collect results in a list and save

res.senti <- 
  list(fit = rfFit,
       pred = predict.train(rfFit, test[, !names(test) %in% c("Type")], type = "raw"),
       confusion = confusionMatrix(rfPred, test$Type),
       varImp = varImp(rfFit, scale = FALSE))

sink(file = 'output/IV_res_senti.txt')
res.senti
sink()

saveRDS(res.senti, 'output/IV_res_senti.rds')
res.senti <- readRDS("output/IV_res_senti.rds")

# Sequoia Classification

# Partition Data

inTraining <- createDataPartition(seqdata$Type, p = .75, list = FALSE)
train <- seqdata[ inTraining,]
test  <- seqdata[-inTraining,]

# Tune random forest resample process to create variable samples for each tree

rfControl <- trainControl(
  method = "boot",
  number = 100
)

# Approximate mtry (# of vars randomly sampled as candidates at each split)

rfGrid <- expand.grid(mtry = trunc(sqrt(ncol(train)))) 

# RF Model Training

tic("RF") #Start timing

rfFit <- train(Type ~ ., data = train,
               method = "rf",
               importance = TRUE, ntree=500,
               trControl = rfControl, tuneGrid = rfGrid,
               metric = "Accuracy", maximize = TRUE) #runs ~0.14h

rfFit #Model output based on training data


# RF Model Testing

rfPred <- 
  predict.train(rfFit, test[, !names(test) %in% c("Type")], type = "raw")

toc() # Stop timing

# Collect results in a list and save

res.seq <- 
  list(fit = rfFit,
       pred = predict.train(rfFit, test[, !names(test) %in% c("Type")], type = "raw"),
       confusion = confusionMatrix(rfPred, test$Type),
       varImp = varImp(rfFit, scale = FALSE))

sink(file = 'output/V_res_seq.txt')
res.seq
sink()

saveRDS(res.seq, 'output/V_res_seq.rds')
res.seq <- readRDS("output/V_res_seq.rds")