library(survival)
library(plyr)
library(DMwR)
library(boot)
library(mlogit)
library(VGAM)
library("kernlab", lib.loc="~/R/win-library/3.1")

## fonction pour avoir le nombre d'unique d'une colonne
funk_nb_unique <- function(data,esc){
  ## data: data.frame , esc:string
  return(paste0(paste(esc,"nb of unique : "),NROW(unique(data[,esc]))))
}

## fonction pour avoir le nombre de na dans une colonne
funk_nb_na <- function(data,esc){
  ## data: data.frame , esc:string
  return(paste0(paste(esc,"nb of NA : "),sum(is.na(data[,esc]))))
}

## fonction pour avoir le nombre valeur manquante dans une colonne
funk_nb_missings <- function(data,esc){
  ## data: data.frame , esc:string
  return(paste0(paste(esc,"nb of missing values : "),NROW(data[(data[,esc] ==""),esc])))
}

## fonction pour binarizer une colonne pour les char:
funk_binarize_char <-function(data,esc){
  res = data.frame(data[,esc],stringsAsFactors = FALSE)
  names(res) = esc
  res[,1] = sapply(res[,1],as.character)
  unik = unique(res[,1])
  k = 0
  for(iter in unik){
    res[(res[,esc] == iter),esc] = k
    k = k + 1
  }    
  return(res)
}


## function permettant l'évaluation des NA
funk_check_NA <- function(base)
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
  pct.missing_2 = pct.missing_2[with(pct.missing_2, order(Pourcentage_NA,decreasing = T)),]
  print(pct.missing_2)
  pct.missing_2[,2] = sapply(pct.missing_2[,2],as.character)
  pct.missing_2[,2] = sapply(pct.missing_2[,2],as.numeric)
  pct.missing_2[,1] = sapply(pct.missing_2[,1],as.character)
  return(pct.missing_2)
}



## function de passage des factors aux numeriques
funk_num <-  function(data,dico){
  num_name = na.omit(dico[(dico[,"Type"] == "NUM") & (dico[,"Table"] == "CoreTable"),c("Var.Name")])
  num_name = sapply(num_name,as.character)
  res = data[,num_name]
  for(esc in num_name){
    res[,esc]=sapply(res[,esc],as.character)
    res[,esc]=sapply(res[,esc],as.numeric)
  }
  return(res)
}


## function de passage des factors aux numeriques
funk_char <-  function(data,dico){
  num_name = na.omit(dico[(dico[,"Type"] == "CHAR") & (dico[,"Table"] == "CoreTable"),c("Var.Name")])
  num_name = sapply(num_name,as.character)
  res = data[,num_name]
  for(esc in num_name){
    res[,esc]=sapply(res[,esc],as.character)
  }
  return(res)
}



