# Spectral discriminaton of invasive plant species: Hawkweed
  
# This code analysis spectral profiles of plant species from Mount Kosciuszko
# National Park, NSW, Australia. Initially, outlying spectra from each species 
# will be removed. Then, a Random Forest classifier will be trained and validated. Follwoing questions will be
# addressed:
#   
# Can spectral profiles of all species be accurately classified? 
# Is the classification still accurate if spectral profiles are resampled according to the specifications of following sensors:
# Sequoia (drone sensor) 
# Landsat-8
# Sentinel-2

## Set up working environment

dir.create('output', FALSE, FALSE)

# install.packages(c("dplyr",
#                    "hsdar", 
#                    "fda",
#                    "fda.usc",
#                    "prospectr",
#                    "gdata",
#                    "caret",
#                    "reshape2",
#                    "cowplot",
#                    "ggplot2",
#                    "tictoc"))

library(dplyr)
library(hsdar) #hyperspectral data processing
library(fda) #outlier detection
library(fda.usc) #outlier detection
library(prospectr) #spectral binning
library(gdata) #drop factor levels
library(caret) #classification
library(reshape2) #reformat wide/long
library(cowplot) #customize ggplot2 prints
library(tictoc) #record time 

source('R/FUN_drop_cat_var.R') 
source('R/FUN_raw2speclibhsdar.R')
source('R/FUN_Remove_Functional_Outlier.R')
source('R/FUN_prepggwide2long.R')

## Cleaning data

data.original <- read.csv("data/All samples_converted.csv", check.names = FALSE)

data.rmv.noise <- data.original[,match('400', 
                  names(data.original)):match('2400', 
                  names(data.original))] #Rmv bands above and below given numeric values

names(data.rmv.noise[,c(1,ncol(data.rmv.noise))]) #Check bands after removal

data.wo.noise <- cbind(data.original['Type'],
                       data.rmv.noise) #Build final df


subsets <- split(data.wo.noise, data.wo.noise$Type)#split df by type into list

plots.bef <- list()
for (i in 1:length(subsets)){#plot spectra before outlier detection
  
  p <- prep_gg(as.data.frame(subsets[[i]]))
  
  plots.bef[[i]] <- ggplot(p, aes(Wavelength, Reflectance, colour = Type)) +
    geom_point(aes(shape=Type), size = .1)+
    labs(title=paste(names(subsets[i])), x= "WL", y="R")+
    theme_minimal()+
    theme(text=element_text(size=8))+
    theme(legend.position="none")
}

p.bef <- plot_grid(plotlist=plots.bef)

ggsave("output/spectrabefore.pdf",
       plot = p.bef,
       width = 40,
       height = 20,
       units = "cm",
       dpi = 50
)

cleaned_data <- lapply(subsets, rmv.funct.outlier,depth.mode, 1, 0.05, 0.5)

plots <- list()
for (i in 1:length(cleaned_data)){#plot spectra after outlier detection
  
  p <- prep_gg(as.data.frame(cleaned_data[[i]]), agg = FALSE)
  
  plots[[i]] <- ggplot(p, aes(Wavelength, Reflectance, colour = Type)) +
    geom_point(aes(shape=Type), size = .1)+
    labs(title=paste(names(cleaned_data[i])), x= "WL", y="R")+
    theme_minimal()+
    theme(text=element_text(size=8))+
    theme(legend.position="none")
}


p.aft <- plot_grid(plotlist=plots)

ggsave("output/spectraafter.pdf",
       plot = p.aft,
       width = 40,
       height = 20,
       units = "cm",
       dpi = 50
)


cleaned.df <- bind_rows(cleaned_data)

write.csv(cleaned.df, 'output/20180622_clean_spectra.csv', row.names = FALSE)

rmved <- nrow(data.original)-nrow(cleaned.df)
paste(rmved, "spectral profiles removed.")

## Classification of all spectra

cleaned.df <- read.csv("output/20180622_clean_spectra.csv", check.names = FALSE)

inTraining <- createDataPartition(cleaned.df$Type, p = .75, list = FALSE)
train <- cleaned.df[ inTraining,]
test  <- cleaned.df[-inTraining,]

options(warn=-1)

subsets <- c(5, 50, 500)

ctrl <- rfeControl(functions = rfFuncs,
                   method = "repeatedcv",
                   repeats = 5,
                   verbose = FALSE)

rfProfile <- rfe(x=train[,-1], y=train[,1],
                 sizes = subsets,
                 rfeControl = ctrl)

sink("output/I_rfe.txt")
rfProfile
sink()

saveRDS(rfProfile, 'output/I_rfe.rds')

# Define the training control
fitControl <- trainControl(
  method = 'boot',                   # k-fold cross validation
  number = 5,                      # number of folds
  savePredictions = 'final',       # saves predictions for optimal                                           tuning parameter
  classProbs = T,                  # should class probabilities be                                           returned
  summaryFunction=multiClassSummary  # results summary function
) 

rfGrid <- expand.grid(mtry = trunc(sqrt(ncol(train)))) 

tic("train")
# RF Model Training

rfFit <- train(Type~.,
               train,method = "rf",
               importance = TRUE, 
               ntree=500,
               trControl = fitControl, 
               tuneGrid = rfGrid,
               metric = "Accuracy", 
               maximize = TRUE) 

toc()

# RF Model Testing

rfPred <- 
  predict.train(rfFit, test[, !names(test) %in% c("Type")], type = "raw")

res.all <- 
  list(fit = rfFit,
       pred = predict.train(rfFit, test[, !names(test) %in% c("Type")],        type = "raw"),
       confusion = confusionMatrix(rfPred, test$Type),
       varImp = varImp(rfFit, scale = FALSE))

sink(file = 'output/I_res_all.txt')
res.all
sink()

saveRDS(res.all, 'output/I_res_all.rds')

## Classification of leaves only

leaf <- cleaned.df[grep("eaf", cleaned.df$Type), ]
leaf$Type <- drop.levels(leaf$Type)
unique(leaf$Type)

inTraining <- createDataPartition(leaf$Type, p = .75, list = FALSE)
train.l <- leaf[ inTraining,]
test.l  <- leaf[-inTraining,]

options(warn=-1)

subsets <- c(5, 50, 500)

ctrl <- rfeControl(functions = rfFuncs,
                   method = "repeatedcv",
                   repeats = 5,
                   verbose = FALSE)

rfProfile.l <- rfe(x=train.l[,-1], y=train.l[,1],
                   sizes = subsets,
                   rfeControl = ctrl)

sink("output/II_rfe.l.txt")
rfProfile.l
sink()

saveRDS(rfProfile.l, 'output/II_rfe.l.rds')


rfGrid <- expand.grid(mtry = trunc(sqrt(ncol(train.l)))) 
tic("train.l")

# RF Model Training

rfFit.l <- train(Type~., 
                 train.l,
                 method = "rf",
                 importance = TRUE, 
                 ntree=500,
                 trControl = fitControl, 
                 tuneGrid = rfGrid,
                 metric = "Accuracy", 
                 maximize = TRUE) 
toc()

# RF Model Testing

rfPred.l <- 
  predict.train(rfFit.l, test.l[, !names(test.l) %in% c("Type")], type = "raw")

res.l <- 
  list(fit = rfFit.l,
       pred = predict.train(rfFit.l, test.l[, !names(test.l) %in% c("Type")],        type = "raw"),
       confusion = confusionMatrix(rfPred.l, test.l$Type),
       varImp = varImp(rfFit.l, scale = FALSE))

sink(file = 'output/II_res.l.txt')
res.l
sink()

saveRDS(res.l, 'output/II_res.l.rds')

## Classification flower only

flower <- cleaned.df[grep("ower", cleaned.df$Type), ]
flower$Type <- drop.levels(flower$Type)
unique(flower$Type)


inTraining <- createDataPartition(flower$Type, p = .75, list = FALSE)
train.f <- flower[ inTraining,]
test.f  <- flower[-inTraining,]

options(warn=-1)

subsets <- c(5, 50, 500)

ctrl <- rfeControl(functions = rfFuncs,
                   method = "repeatedcv",
                   repeats = 5,
                   verbose = FALSE)

rfProfile.f <- rfe(x=train.f[,-1], y=train.f[,1],
                   sizes = subsets,
                   rfeControl = ctrl)

sink("output/III_rfe.f.txt")
rfProfile.f
sink()

saveRDS(rfProfile.f, 'output/III_rfe.f.rds')


rfGrid <- expand.grid(mtry = trunc(sqrt(ncol(train.f)))) 
tic("train.f")
# RF Model Training

rfFit.f <- train(Type~., 
                 train.f,
                 method = "rf",
                 importance = TRUE, 
                 ntree=500,
                 trControl = fitControl, 
                 tuneGrid = rfGrid,
                 metric = "Accuracy", 
                 maximize = TRUE) 

toc()


# RF Model Testing

rfPred.f <- 
  predict.train(rfFit.f, test.f[, !names(test.f) %in% c("Type")], type = "raw")

res.f <- 
  list(fit = rfFit.f,
       pred = predict.train(rfFit.f, test.f[, !names(test.f) %in% c("Type")],        type = "raw"),
       confusion = confusionMatrix(rfPred.f, test.f$Type),
       varImp = varImp(rfFit.f, scale = FALSE))

sink(file = 'output/III_res.f.txt')
res.f
sink()

saveRDS(res.f, 'output/III_res.f.rds')

## Resample spectral data to other sensor specifications

hypdata <- cleaned.df

# Create spectral library to use hsdar pkg

speclib <- raw2speclib(hypdata) #Function requires just numbers as colnames.
plot(speclib)

center <-  c(550, 660, 735, 790)
fwhm <- c(20, 20, 5, 20)

sequoia <- as.data.frame(cbind(center, fwhm))

data_seq <- spectralResampling(speclib, sequoia)
plot(data_seq)

data_senti2 <- spectralResampling(speclib, 'Sentinel2')
plot(data_senti2)

seqdata <- as.data.frame(data_seq@spectra@spectra_ma)
seqdata <- cbind('Type'=hypdata$Type, seqdata)
str(seqdata)
names(seqdata)
newnamesSeq <- c("Type", "Green", "Red", "RedEdge", "NIR")
names(seqdata) <- newnamesSeq


sentidata <- as.data.frame(data_senti2@spectra@spectra_ma)
sentidata <- cbind('Type'=hypdata$Type, sentidata)
str(sentidata)
names(sentidata)
newnamesSenti <- c("Type", "B1Aero", "B2Blue", "B3Green", "B4Red",
                   "B5RE1", "B6RE2", "B7RE3", "B8NIR", "B9WaterVap",
                   "B10SWIR1", "B11SWIR2", "B12SWIR3", "B13")

names(sentidata) <- newnamesSenti


## Sentinel2 Classification

inTraining <- createDataPartition(sentidata$Type, p = .75, list = FALSE)
train.sen <- sentidata[ inTraining,]
test.sen  <- sentidata[-inTraining,]

options(warn=-1)

subsets <- c(5, 50, 500)

ctrl <- rfeControl(functions = rfFuncs,
                   method = "repeatedcv",
                   repeats = 5,
                   verbose = FALSE)

rfProfile.sen <- rfe(x=train.sen[,-1], y=train.sen[,1],
                     sizes = subsets,
                     rfeControl = ctrl)

sink("output/IV_rfe.sen.txt")
rfProfile.sen
sink()

saveRDS(rfProfile.sen, 'output/IV_rfe.sen.rds')

rfGrid <- expand.grid(mtry = trunc(sqrt(ncol(train.sen)))) 
tic("train.sen")

# RF Model Training

rfFit.sen <- train(Type~., 
                   train.sen,
                   method = "rf",
                   importance = TRUE, 
                   ntree=500,
                   trControl = fitControl, 
                   tuneGrid = rfGrid,
                   metric = "Accuracy", 
                   maximize = TRUE) 
toc()

# RF Model Testing

rfPred.sen <- 
  predict.train(rfFit.sen, test.sen[, !names(test.sen) %in% c("Type")], type = "raw")

res.sen <- 
  list(fit = rfFit.sen,
       pred = predict.train(rfFit.sen, test.sen[, !names(test.sen) %in% c("Type")],        type = "raw"),
       confusion = confusionMatrix(rfPred.sen, test.sen$Type),
       varImp = varImp(rfFit.sen, scale = FALSE))

sink(file = 'output/IV_res.sen.txt')
res.sen
sink()

saveRDS(res.sen, 'output/IV_res.sen.rds')

## Sequoia Classification

inTraining <- createDataPartition(seqdata$Type, p = .75, list = FALSE)
train.seq <- seqdata[ inTraining,]
test.seq  <- seqdata[-inTraining,]

options(warn=-1)

subsets <- c(5, 50, 500)

ctrl <- rfeControl(functions = rfFuncs,
                   method = "repeatedcv",
                   repeats = 5,
                   verbose = FALSE)

rfProfile.seq <- rfe(x=train.seq[,-1], y=train.seq[,1],
                     sizes = subsets,
                     rfeControl = ctrl)

sink("output/V_rfe.seq.txt")
rfProfile.seq
sink()

saveRDS(rfProfile.seq, 'output/V_rfe.seq.rds')

rfGrid <- expand.grid(mtry = trunc(sqrt(ncol(train.seq)))) 
tic("train.seq")
# RF Model Training

rfFit.seq <- train(Type~., 
                   train.seq,
                   method = "rf",
                   importance = TRUE, 
                   ntree=500,
                   trControl = fitControl, 
                   tuneGrid = rfGrid,
                   metric = "Accuracy", 
                   maximize = TRUE) 
toc()

# RF Model Testing

rfPred.seq <- 
  predict.train(rfFit.seq, test.seq[, !names(test.seq) %in% c("Type")], type = "raw")

res.seq <- 
  list(fit = rfFit.seq,
       pred = predict.train(rfFit.seq, test.seq[, !names(test.seq) %in% c("Type")],        type = "raw"),
       confusion = confusionMatrix(rfPred.seq, test.seq$Type),
       varImp = varImp(rfFit.seq, scale = FALSE))

sink(file = 'output/V_res.seq.txt')
res.seq
sink()

saveRDS(res.seq, 'output/V_res.seq.rds')

#End