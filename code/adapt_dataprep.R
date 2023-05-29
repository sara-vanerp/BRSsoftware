## Bayesian Regularized SEM Software: Data prep ADAPT illustration
## Author: Sara van Erp

library(MASS)
library(lavaan)
library(psych)

set.seed(06032023)

## Data: generated based on the observed covariance matrix from the ADAPT study
covAdpt <- read.csv2("./data/obscov_adapt.csv") 
simDat <- mvrnorm(n = 1000,
                  mu = rep(0, 65),
                  Sigma = covAdpt)
colnames(simDat) <- colnames(covAdpt)

## EFA
# similarly to the ADAPT study, first run an EFA on the training set, followed by a (partial) CFA on the test set
trainDat <- simDat[1:748, ]
fa.parallel(trainDat)

fitEFA <- fa(r = trainDat,
             nfactors = 8,
             rotate = "oblimin",
             fm = "pa")
print(fitEFA, sort = TRUE)

## For subsequent analyses, we only use the first three factors, to make it more practical to run.
## All items are used for a factor that have the highest loading that is still > 0.20. 
## Although a main loading of 0.20 or 0.30 is not ideal, this is the type of real data setting that we wish to test the methods on
## In addition, some items have quite substantial cross-loadings (e.g. AVVB21)

## CFA
# model based on standardized loadings from EFA
# F1: 18; F2: 6; F3: 12
testDat <- simDat[749:nrow(simDat), ]
modCFA <- "F1 =~ AVVB52 + AVVB56 + AVVB54 + AVVB57 + AVVB61 + AVVB59 + AVVB55 + AVVB60 +
           AVVB48 + AVVB49 + AVVB35b + AVVB50 + AVVB63 + AVVB19 + AVVB62 + AVVB58 + AVVB51 + AVVB20 
           F2 =~ AVVB3 + AVVB4 + AVVB2 + AVVB5 + AVVB9 + AVVB8
           F3 =~ AVVB31 + AVVB32 + AVVB30 + AVVB34 + AVVB22 + AVVB53 + AVVB45 + AVVB15 + AVVB43 + AVVB35a + AVVB46 + AVVB21" 

fitCFA <- cfa(modCFA, data = testDat)
summary(fitCFA, fit.measures = TRUE) # model does not fit very well
modindices(fitCFA, sort = TRUE) # 1st MI corresponds to two items (31 and 32) that are very similar, one of which was removed in the original study

## Change names
# change variable names to make the model specification in Mplus and the regularization packages easier
# order variables based on which factor they have the main loading
ord <- c("AVVB52", "AVVB56", "AVVB54", "AVVB57", "AVVB61", "AVVB59", "AVVB55", "AVVB60", "AVVB48", "AVVB49", "AVVB35b", 
  "AVVB50", "AVVB63", "AVVB19", "AVVB62", "AVVB58", "AVVB51", "AVVB20", "AVVB3", "AVVB4", "AVVB2", "AVVB5", "AVVB9", "AVVB8",
   "AVVB31", "AVVB32", "AVVB30", "AVVB34", "AVVB22", "AVVB53", "AVVB45", "AVVB15", "AVVB43", "AVVB35a", "AVVB46", "AVVB21")

testDatOrd <- testDat[, ord]
colnames(testDatOrd) <- paste0("y", 1:ncol(testDatOrd))
trainDatOrd <- trainDat[, ord]
colnames(trainDatOrd) <- paste0("y", 1:ncol(trainDatOrd))

## Standardization
# before applying any type of regularization, it is important to standardize the data so that the penalization affects each path similarly
testDatSD <- scale(testDatOrd, center = TRUE, scale = TRUE)
save(testDatSD, file = "./data/testDatSD_adapt.dat")
MplusAutomation::prepareMplusData(as.data.frame(testDatSD), filename = "./data/testDatSD_adapt_Mplus.dat") # save data for Mplus
