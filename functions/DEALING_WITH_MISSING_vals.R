## packages 
# install.packages("data.table") ## if need be
library(data.table)
library(survival)
library(plyr)
library(DMwR)
library(boot)
library(mlogit)
library(VGAM)


## indicatgin path
setwd("C:\\cygwin64\\home\\MDjerrab\\AppGenomique\\data")

## sourcer les fonctions nécaissaires
source("../functions/function_data.R")



## DATA
data  = read.table("CoreTable_training.csv",header = T,sep = ";")
dico = read.table("dico_table.csv",header = T, sep = ";")
#########################################################################
######################## TRAVAIL A EFFECTUER ############################
#########################################################################

# DATA CLEANING :
#1# #funk ok# convertir grace au dico: num en num et char en char.
#2# #funk ok# spoter les colonnes avec des na et rejeter celle au
            # dessus de 50% na
#3# #funk ok# spoter les colonnes avec des uniques = 1 => les supprimer

# DATA RECONSTRUCTION :
#4# #funk ok# spoter les colonnes avec des uniques = 2 => binarize 
            # (0 ou 1 mais réfléchir sur l'ordre implicite)
#5# #funk ok# spoter les colonnes avec des missings 
            # voir si complétion ou rejet des colonnes
#6# # Remplacer les missings et les NA avec les méthodes adéquates
#6# bis # Supprimer les colonnes intraitables

# SAVE CLEAN DATA:
#########################################################################

# DATA CLEANING :
#1# Convercion factor en character et unique
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
                                      ## On va juste modifier la colonne smoke
for(esc in test[,1]){char[(char[,esc] == "."),esc] = NA}
char$SMOKE[(char$SMOKE == "YES")] = 1
char$SMOKE[(char$SMOKE == "NO")] = 0
#Don't smoke 1:Smke

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
    # Il reste 4 colonnes à traiter à la main
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

    # Définir les variables à compléter
test = funk_check_missings(data)
name_to_complete = test[(test[,2] !=0),1]
test = funk_check_NA(data)
name_to_complete = unique(c(name_to_complete,test[(test[,2] !=0),1]))

    # on boucle sur les variables à compléter: complétion pour les num
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
    # On s'occupe des cols char non complete
    ## funk_check_missings(data[,name_char])
    ##   Variable pourcentage_missings
    ## 1 SMOKSTAT               96.125
    ## 2 SMOKFREQ              92.8125
    ## 3 TSTAG_DX              72.3125
    ## 4    SMOKE               70.375 ## We'll complete only this one

    # Remove the other
data <- funk_remove_names(data,name_char[1:3])
    #Complete SMOKE
test_reg = name_num_complete[sort.list(kant,decreasing = T)[1]]
test =data[,c("SMOKE",test_reg[1])]
test = test[(test$SMOKE =="1")|(test$SMOKE == "0"),]
test$SMOKE = as.numeric(test$SMOKE)
test1 = summary(glm(SMOKE~.,test,family = "binomial"))
test_bis = data[(data$SMOKE !="1")&(data$SMOKE != "0"),c("SMOKE",test_reg[1])]
test2 = test1$coefficients[1,1]+ 
        test1$coefficients[2,1]*test_bis[,2] 
test2[test2 == unique(test2[!is.na(test2)])[1]]<-0
test2[test2 == unique(test2[!is.na(test2)])[2]]<-1
test2[is.na(test2)]<-0
data$SMOKE = as.numeric(data$SMOKE)
data$SMOKE[is.na(data$SMOKE)] = test2

# on vérifie qu'il n'y a plus de missings et de na
funk_check_NA(data)
funk_check_missings(data)


## SAVE DATA CLEANED
write.csv2(data,"DATA_CLEAN.csv")




