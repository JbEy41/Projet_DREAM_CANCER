library(plyr)
require(plyr)
library(DMwR)
library(boot)
library(mlogit)
library(VGAM)


## indicatgin path
setwd("C:\\cygwin64\\home\\MDjerrab\\AppGenomique\\data\\leaderboard\\")
## sourcer les fonctions nécaissaires
source("../../functions/function_data.R")


## DATA
data_train  =read.table("../CoreTable_training.csv",header = T,sep = ";")
data  = read.table("CoreTable_leaderboard.csv",header = T, sep = ",")
dico = read.table("../dico_table.csv",header = T, sep = ";")



##########################################################################
##########################################################################
##########################################################################

num = funk_num(data,dico)
char = funk_char(data,dico)
# data = cbind(num,char) ## concatene après tout le travail
## effectuée

#2# 1 ) Suppression des colonnes NA de la table num
test =funk_check_NA(base = num)
name_to_remove = test[(test[,2]>49),1]
num <- funk_remove_names(num,name_to_remove)

#2# 2 ) Suppression des colonnes NA de la table char
test = funk_check_NA(char)
name_to_remove = test[(test[,2] > 49),1]
char <- funk_remove_names(char,name_to_remove)

#3# 1) spoter les colonnes avec des uniques = 1 => les supprimer : num
test = funk_check_nb_unique(num)
name_to_remove = test[(test[,2] == 1),1]
num <- funk_remove_names(num,name_to_remove)

#3# 2) spoter les colonnes avec des uniques = 1 => les supprimer : char
test = funk_check_nb_unique(char)
name_to_remove = test[(test[,2] == 1),1]
char <- funk_remove_names(char,name_to_remove)


# DATA RECONSTRUCTION :
#4# spoter les colonnes avec 2 uniques et binarizer (Seulement pour les char)
# On vérifie que les colonnes avec plus de 2-uniks mais moins de 10
# ne sont pas simplement des 2-uniks
test = funk_check_nb_unique(char)
test[,2]=sapply(test[,2],as.numeric)
test_1 = test[(test[,2]>2) & (test[,2]<11),1]
test = funk_uniks_show(char[,test_1]) ## On remplace les "." point par NA
## On va juste modifier la colonne smoke  => pas de smoke

test = funk_check_binar(char) ## pour voir comment binariser
# le colonnes en Yes et missings vont être binariser avec Yes =1
test_1 = c(grep("Y", test[,2]),grep("YES", test[,2]))
test_1 = unique(c(test_1,grep("Y", test[,3]),grep("YES", test[,3])))
test_1 = test[test_1 ,1] 
for(esc in test_1){
  char[(char[,esc] == "Y"),esc] = 1
  char[(char[,esc] == "YES"),esc] = 1
  char[(char[,esc] == ""),esc] = 0
}
# Il reste 4 colonnes �  traiter �  la main
test = funk_check_binar(char) ## pour voir comment binariser
names_to_binarize = test[(test[,2] != "0")&(test[,2] != "1"),1]
for(esc in names_to_binarize){
  char[,esc] = funk_binarize_char(char,esc) ## la conversion s'effetue
  ## Dans le sens d'apparition
}

# On convertit en numérique
for(esc in test[,1]){
  char[,esc] =   sapply(char[,esc] ,as.numeric)
}

# On reconcatene la table DATA avant la complétion
data = cbind(num,char)

##5# #funk ok# spoter les colonnes avec des missings 
# voir si complétion ou rejet des colonnes
# faire une différence entre char et num
num = data[,grep("numeric",sapply(data,class))]
char = data[,sapply(data,class) =="character"]

# Définir la liste des varaibles pour la complétions
test = funk_check_missings(num)
name_num_complete = test[(test[,2] ==0),1]
test = funk_check_NA(num)
name_num_complete = unique(c(name_num_complete,test[(test[,2] ==0),1]))

# Définir les variables �  compléter
test = funk_check_missings(data)
name_to_complete = test[(test[,2] !=0),1]
test = funk_check_NA(data)
name_to_complete = unique(c(name_to_complete,test[(test[,2] !=0),1]))

# on boucle sur les variables �  compléter: complétion pour les num
name_char=c()
for(esc in name_to_complete){
  #esc = name_to_complete[5]
  if(is.numeric(data[,esc])){
    kant=c()
    for(esc1 in name_num_complete){
      
      kant = c(kant,cor(data[,esc],
                        data[,esc1],
                        use="complete.obs"))
    }
    test_reg = name_num_complete[sort.list(kant,decreasing = T)[1]]
    test1 = summary(glm(data[,esc] ~ data[,test_reg],na.action = na.omit))
    test1$coefficients[1]
    test1$coefficients[2]
    test1= round(data[(data[,esc] == "")|(is.na(data[,esc])),
                      test_reg]*test1$coefficients[2]+
                   test1$coefficients[1])
    #test1[test1>1]=1 ## voir si cela pose problem
    data[(data[,esc] == "")|(is.na(data[,esc])),esc] = test1
  }else{  
    name_char = c(name_char,esc)
  }
}



# on vérifie qu'il n'y a plus de missings et de na
funk_check_NA(data)
funk_check_missings(data)

head(data)

write.csv2(data,"DATA_CLEAN_PREDICT.csv")






























