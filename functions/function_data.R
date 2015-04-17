library(survival)
library(plyr)
library(DMwR)
library(boot)
library(mlogit)
library(VGAM)
library("kernlab", lib.loc="~/R/win-library/4.1")

## fonction pour avoir le nombre d'unique d'une colonne
funk_nb_unique <- function(data,esc){
  ## data: data.frame , esc:string
  # print(paste0(paste(esc,"nb of unique : "),NROW(unique(data[,esc]))))
  return(NROW(unique(data[,esc])))
}

## fonction pour avoir le nombre d'unique par colonne dans une table
funk_check_nb_unique <- function(data){    
  # data : data.frame
  test = data.frame(matrix(0,NROW(names(data)),2))
  names(test) = c("Variable","NB_unique")
  for(esc in c(1:(NROW(names(data))))){
    test[esc,] <- c(names(data)[esc],
                    funk_nb_unique(data,names(data)[esc]))
  }
  test = test[with(test, order(NB_unique,decreasing = T)),]
  return(test)
}

## fonction pour avoir le nombre de na dans une colonne
funk_nb_na <- function(data,esc){
  ## data: data.frame , esc:string
  # print(paste0(paste(esc,"nb of NA : "),sum(is.na(data[,esc]))))
  return(sum(is.na(data[,esc])))
}

## fonction pour avoir le nombre valeur manquante dans une colonne
funk_nb_missings <- function(data,esc){
  ## data: data.frame , esc:string
  # print(paste0(paste(esc,"nb of missing values : "),
  #             NROW(data[(data[,esc] ==""),esc])))
  return(NROW(data[(data[,esc] ==""),esc]))
}

## fonction pour binarizer une colonne pour les char:
funk_binarize_char <-function(data,esc){
  # data : data.frame , esc : string
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
  # base : data.frame
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

# function permettant l'évaluation des missings
funk_check_missings <- function(data){    
  # data : data.frame
  test = data.frame(matrix(0,NROW(names(data)),2))
  names(test) = c("Variable","pourcentage_missings")
  for(esc in c(1:(NROW(names(data))))){
    test[esc,] <- c(names(data)[esc],
                    (funk_nb_missings(data,esc)*100/NROW(data)))
  }
  test = test[with(test, order(pourcentage_missings,decreasing = T)),]
  return(test)
}

## function de passage des factors aux numeriques
funk_num <-  function(data,dico){
  #data : data.frame, dico : data.frame
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
  #data : data.frame, dico : data.frame
  num_name = na.omit(dico[(dico[,"Type"] == "CHAR") & (dico[,"Table"] == "CoreTable"),c("Var.Name")])
  num_name = sapply(num_name,as.character)
  res = data[,num_name]
  for(esc in num_name){
    res[,esc]=sapply(res[,esc],as.character)
  }
  return(res)
}

## Fonction supprimant les colonnes d'un dataframe par string
funk_remove_names <- function(data,name_to_remove){
  #data : data.frame, name_to_remove : list of strings
  tl = c()
  for(esc in name_to_remove){
    tl = c(tl,grep(esc ,fixed = T, names(data)))
  }
  if(length(tl)>0){
    data <- data[,-tl]
  }else{
    print("no column to remove")
  }
  return(data)
}
## fonction nous donnons un tableaux des valeurs des colonnes à binariser
funk_check_binar <- function(data){    
  # data : data.frame
  test = funk_check_nb_unique(data)
  name_to_binarize = test[(test[,2] == 2),1]
  test = data.frame(matrix(0,NROW(name_to_binarize),3))
  names(test) = c("Variable","unik1","unik2")
  for(esc in c(1:(NROW(name_to_binarize)))){
    test_1 = unique(data[,name_to_binarize[esc]])
    test[esc,] <- c(name_to_binarize[esc],
                    test_1[1],
                    test_1[2]                        
    )
  }
  return(test)
}

## fonction nous donnos un tableaux des uniks des colonnes binariser ou non
funk_uniks_show <- function(data){
  max_unik = sapply(char[,test_1],unique)
  max_unik =max(sapply(max_unik,NROW))
  test = data.frame(matrix(0,NCOL(data),max_unik+1))
  nam = c()
  for(esc in c(1:(max_unik))){ nam <- c(nam,paste0("unik_",esc))}
  names(test) = c("Variable",nam)
  for(esc in c(1:(NROW(test)))){      
    test_1 = unique(data[,esc])
    test[esc,1] = names(data)[esc]
    for(esc1 in c(1:(NROW(test_1)))){
      test[esc,esc1+1] <- test_1[esc1]              
    }
  }
  return(test)
}