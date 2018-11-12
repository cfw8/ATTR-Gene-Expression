

library(readr)
library(dplyr)
library(ggplot2)
library(randomForest) 
library(Hmisc)
library(caret)


## Read in data from directory & merge data files

geneExpressionData <- read_tsv("GSE67784.txt.gz")
colnames(geneExpressionData)[1] <- "SampleID"

metaData <- read_csv("Metadata/Ensembl87_Biomart_GRCh38.p7.csv")
colnames(metaData)[1] <- "Gene"
colnames(metaData)[2] <- "Chromosome"
colnames(metaData)[3] <- "GeneSymbol"

# Remove unnecessary column
metaData <- select(metaData, -4)

metaData <- filter(metaData, metaData$Chromosome %in% c(1:22, "X", "Y"), !is.na(metaData$GeneSymbol))
# filteredMetaData <- filter(metaData, metaData$Chromosome %in% c(1:22))

metaData$Chromosome <- factor(metaData$Chromosome)
# filteredMetaData$Chromosome <- factor(filteredMetaData$Chromosome)

## Write a function to create separate datasets

getClassData <- function(geneExpressionData, classFile) {
  classData <- read_tsv(classFile)
  classData$Class <- factor(classData$Class)
  
  mergedData <- inner_join(classData, geneExpressionData)
  
  return(mergedData)
}


## Create separate datasets, combine them

asymptomaticMergedData <- getClassData(geneExpressionData, "Class/Asymptomatic__gender.txt")
controlsMergedData <- getClassData(geneExpressionData, "Class/Only_Controls__gender.txt")

asymptomaticData <- cbind("A", asymptomaticMergedData)
colnames(asymptomaticData)[1] <- "Type"
controlsData <- cbind("C", controlsMergedData)
colnames(controlsData)[1] <- "Type"

fullData <- rbind(asymptomaticData, controlsData) 




## Splitting the data into training and test data

splitdf <- function(data) {
  index <- 1:nrow(data)
  trainindex <- sample(index, length(index)*.5)
  trainset <- data[trainindex, ]
  testset <- data[-trainindex, ]
  res <- list(trainset=trainset,testset=testset)
  return(res)
}


## RANDOM FOREST ##

## Only predicting between control and asymptomatic (using sig genes by type) with gender and age

# Using significant by type

### Create training and test data
trainandtestData <- splitdf(fullData)
trainData <- trainandtestData[[1]]
testData <- trainandtestData[[2]]

### Training data ###
samp <- colnames(fullData[,4:ncol(fullData)])
predictors <- trainData[, c("Class", samp)]

# Run RF many times for a monte carlo estimate
correctRate1 <- NULL
for (i in 1:1000) {
  fitGE1 <- randomForest(predictors, y=trainData$Type, ntree=600)
  
  ### Test Data ###
  predType1 <- predict(fitGE1, newdata = testData)
  t1 <- table(predType1, as.factor(as.character(testData$Type)))
  correctRate1[i] <- (t1[1,1]+t1[2,2])/sum(t1)
}
mcEst1 <- mean(correctRate1)
CI1 <- quantile(correctRate1, c(0.025,0.975))


## NO GENDER IN MODEL ##

### Training data ###
predictorsNoSex <- trainData[,samp]

# Run RF many times for a monte carlo estimate
correctRate2 <- NULL
for (i in 1:1000) {
  fitGE2 <- randomForest(predictorsNoSex, y=trainData$Type, ntree=600)
  
  ### Test Data ###
  predType2 <- predict(fitGE2, newdata = testData)
  t2 <- table(predType2, as.factor(as.character(testData$Type)))
  correctRate2[i] <- (t2[1,1]+t2[2,2])/sum(t2)
}
mcEst2 <- mean(correctRate2)
CI2 <- quantile(correctRate2, c(0.025,0.975))


## SIMULATED DATA ##

# Permute Type and Class from data
Class <- sample(fullData$Class, 167, replace=FALSE)
Type <- sample(fullData$Type, 167, replace=FALSE)
simData <- as.data.frame(cbind(Type, Class, fullData[,4:ncol(fullData)]))


trainandtestData1 <- splitdf(simData)
trainData1 <- trainandtestData1[[1]]
testData1 <- trainandtestData1[[2]]
predictorsSim <- trainData1[,colnames(simData[,2:ncol(simData)])]

### Train data ###

# Run RF many times for a monte carlo estimate
correctRate3 <- NULL
for (i in 1:1000) {
  fitGEsim <- randomForest(predictorsSim, y=trainData1$Type, ntree=800)
  
  ### Test Data ###
  predType3 <- predict(fitGEsim, newdata = testData1)
  t3 <- table(predType3, as.factor(as.character(testData1$Type)))
  correctRate3[i] <- (t3[1,1]+t3[2,2])/sum(t3)
}
mcSim <- mean(correctRate3)
CI3 <- quantile(correctRate3, c(0.025,0.975))


## Confusion Matrix
confusionMatrix(t1, positive = "A")
confusionMatrix(t2, positive = "A")
confusionMatrix(t3, positive = "A")

## ROC Curve

