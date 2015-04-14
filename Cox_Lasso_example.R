### AppGenomique: Cox Lasso example
### Jeremy L Hour
### 17 mars 2015

### SET WORKING DIRECTORIES

# Jeremy
# setwd("\\paradis\eleves\JL.HOUR\3A\Apprentissage_genomique") # At ENSAE
setwd("/Users/jeremylhour/Dropbox/AppGenomique")

rm(list=ls())

library("Survival")
library("penalized")


### Example given in the penalized tutorial
data(nki70)
attach(nki70)

pen <- penalized(Surv(time, event), penalized = nki70[,8:77],
                 unpenalized = ~ER+Age+Diam+N+Grade, data = nki70, lambda1 = 10)
show(pen)
coefficients(pen)

coefficients(pen, "penalized")
basehaz(pen)
# A single lasso fit using the clinical risk factors
pen <- penalized(Surv(time, event), penalized = ~ER+Age+Diam+N+Grade,
                 data = nki70, lambda1=10, standardize=TRUE)
# using steps
pen <- penalized(Surv(time, event), penalized = nki70[,8:77],
                 data = nki70, lambda1 = 1,steps = 20)
plotpath(pen)
# A fused lasso fit predicting survival
pen <- penalized(Surv(time, event), penalized = nki70[,8:77], data = nki70,
                 lambda1 = 1, lambda2 = 2, fusedl = TRUE)
plot(coefficients(pen, "all"),type="l",xlab = "probes",ylab = "coefficient value")
plot(predict(pen,penalized=nki70[,8:77]))


### Own example
source("functions/GenerateSurvData.R")
data <- GenerateSurvData(n=500,p=300,s=5)

Cox.model <- penalized(Surv(data$y, data$cens), penalized = data$X,
                      lambda1 = 10, steps=20,
                       model="cox",
                      standardize=T)
show(Cox.model)
coefficients(Cox.model, "penalized")
plotpath(Cox.model)

### Cross-validation
CV.Cox <- optL1(Surv(data$y, data$cens), penalized = data$X,
                    fold = 5,
                    model="cox",
              standardize=T)
coefficients(CV.Cox)
plot(CV.Cox$predictions)

prof <- profL1(Surv(data$y, data$cens), penalized = data$X,
               fold = 5,
               model="cox",
               standardize=T, steps=10)
plot(prof$lambda, prof$cvl, type="l")

### At optimal lambda value
LambdaOpt <- CV.Cox$lambda
Cox.modelOpt <- penalized(Surv(data$y, data$cens), penalized = data$X,
                       lambda1 = LambdaOpt,
                       model="cox",
                       standardize=T)
plot(data$y,lin.pred(Cox.modelOpt))


