library(survival)
library(plyr)
library(DMwR)
library(boot)
library(mlogit)
library(VGAM)
library("kernlab", lib.loc="~/R/win-library/3.1")





## function permettant l'évaluation des NA
CHECK_NA <- function(base)
{  
  names.others <- names(base)
  pct.missing_2 <- matrix(NA,length(base),2)
  for(i in 1:length(base))
  {
    pct.missing_2[i,1] <- names.others[i]
    pct.missing_2[i,2] <- 100*sum(is.na(base[[i]]))/length(base[[i]])
  }
  pct.missing_2=as.data.frame(pct.missing_2)
  pct.missing_2=rename(pct.missing_2, c("V1"="Variable", "V2"="Pourcentage_NA"))
  print(pct.missing_2)
}



