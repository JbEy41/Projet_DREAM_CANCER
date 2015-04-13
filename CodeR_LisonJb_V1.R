#########################################################################################
#                                 DREAM CANCER PROJECT                                  #
#########################################################################################

# Packages
# install.packages("survival") (si nécessaire)
library(survival)
library(plyr)
library(DMwR)

#Chemin Lison : //paradis/eleves/LGrappin/Bureau
# Chemin JB : /Users/eymeoudjean-benoit/Dropbox
CoreTable_training_raw <- read.csv("/Users/eymeoudjean-benoit/Dropbox/AppGenomique/data/CoreTable_training.csv",header=TRUE)  
CoreTable_training <- read.csv("/Users/eymeoudjean-benoit/Dropbox/AppGenomique/data/CoreTable_training.csv",header=TRUE,na.strings=c("","."),dec=".",stringsAsFactors =FALSE)  


#########################################################################################
# 1 PARTIE DATA

# La base est structurée de la manière suivante:
# var 1:29: description de l'individu (taille,poids,cigarette,placebo,..) 
# var 30:54 données biologiques (num)
# var 55:131 historique médical au sens large (logical)



sapply(CoreTable_training_raw,function(x) class(x)) # On regarde le type des données initiales
sapply(CoreTable_training,function(x) class(x)) # On regarde le type une fois que l'on a transformé

#########################
# 1.1 Recodages variables


sapply(CoreTable_training,function(x) class(x)) # Beaucoup de variables sont à changer en numeric

# Death variable
table(CoreTable_training$DEATH, useNA = "always") # Données initiales
CoreTable_training$DEATH[is.na(CoreTable_training$DEATH)] <- 0 # On transforme les NA en NO
CoreTable_training$DEATH[CoreTable_training$DEATH=="YES"] <- 1
table(CoreTable_training$DEATH, useNA = "always") # On regarde que ça marche
CoreTable_training$DEATH=as.numeric(CoreTable_training$DEATH)
table(CoreTable_training$DEATH)

# Region
CoreTable_training$REGION_C[which(CoreTable_training$REGION_C=="MISSING")] <- NA

# Les variables  WGTBLCAT,HGTBLCAT sont des catégories des var WEIGHTBL et HEIGHTBL 
# Je les recode en utilisant des quartiles 
CoreTable_training$WGTBLCAT<-cut(CoreTable_training$WEIGHTBL,quantile(subset(CoreTable_training$WEIGHTBL,is.na(CoreTable_training$WEIGHTBL)==FALSE),probs = seq(0, 1, 0.25)), right=TRUE,include.lowest =TRUE, labels=c(1:4))
CoreTable_training$HGTBLCAT<-cut(CoreTable_training$HEIGHTBL,quantile(subset(CoreTable_training$HEIGHTBL,is.na(CoreTable_training$HEIGHTBL)==FALSE),probs = seq(0, 1, 0.25)), right=TRUE,include.lowest =TRUE, labels=c(1:4))


# On recode les indicatrices de consommation de cigarette
CoreTable_training$SMOKE=as.numeric(as.factor(CoreTable_training$SMOKE))-1 # 0: Don't smoke 1:Smke
CoreTable_training$SMOKFREQ=as.numeric(as.factor(CoreTable_training$SMOKFREQ)) # 1:GREATER THAN OR EQUAL TO 1 PACK PER DAY; 2: LESS THAN 1 PACK PER DAY
CoreTable_training$SMOKSTAT=as.numeric(as.factor(CoreTable_training$SMOKSTAT)) # 1: ONGOING 2: WITHIN THE PAST 12 MONTHS

# On recode d'autres variables mal codées
CoreTable_training$ALP=as.numeric(CoreTable_training$ALP)
CoreTable_training$LDH=as.numeric(CoreTable_training$LDH)
CoreTable_training$AGEGRP=as.numeric(CoreTable_training$AGEGRP) 
CoreTable_training$AGEGRP[is.na(CoreTable_training$AGEGRP)] <-85 # On code les 85 ans et +, 85
CoreTable_training$AGEGRP2<-cut(CoreTable_training$AGEGRP,quantile(CoreTable_training$AGEGRP,probs = seq(0, 1, 1/3)), right=TRUE,include.lowest =TRUE, labels=c(1:3))

table(CoreTable_training$AGEGRP,useNA="always")
table(test,useNA="always")



# On regarde qu'il n'y a pas de valeurs manquantes dans la base initiale pour les variables historiques
apply(CoreTable_training_raw[,55:length(CoreTable_training_raw)],2,function(x) sum(is.na(x)))
# OK


# On recode en binaire les variables indicatrices.
names <- c("DEATH",names(CoreTable_training)[55:length(names(CoreTable_training))])
for(i in names)
{
  CoreTable_training[[i]][is.na(CoreTable_training[[i]])] <- 0
  CoreTable_training[[i]][CoreTable_training[[i]]=="YES"] <- 1
  CoreTable_training[[i]][CoreTable_training[[i]]=="Y"] <- 1
  CoreTable_training[[i]] <- as.numeric(CoreTable_training[[i]])
}


#########################
# 1.2 DETECTIONS DES VALEURS MANQUANTES


################## FONCTION
# Fonctionne 
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
  pct.missing_2=rename(pct.missing_2, c("V1"="Variable", "V2"="Pourcentage NA"))
  print(pct.missing_2)
}
CHECK_NA(CoreTable_training)


# Verif sur BMI
sum(is.na(CoreTable_training$BMI))/length(CoreTable_training$BMI)*100
pct.missing_2[which(pct.missing_2$V1=="BMI"),] # OK


#########################
# 1.3 TRAITEMENT DES VALERUS MANQUANTES

## INDICATRICE DE CONSOMMATION DE TABAC

# Pour l'indicatrice de fumer, on peut regarder les corrélations avec les données biologiques qui ne sont pas trop censurées (>40%)
NA_BIO <- CHECK_NA(CoreTable_training[,30:54]) # On regarde le % de missing sur les données bio
sapply(NA_BIO,function(x) class(x))

name_var=matrix()
for(i in 1:length(CoreTable_training[,30:54]))
{
  if(as.numeric(levels(NA_BIO[i,2]))[NA_BIO[i,2]]<40){name_var[i] <- as.character(NA_BIO[i,1])}
}
name_var # Les variables que l'on prend pas sont marqué NA
name_var=name_var[complete.cases(name_var)==TRUE] # Celles que l'on prend
name_var=c("SMOKE","BMI","RACE_C","AGEGP","WEIGHTBL",name_var)


density(as.numeric(CoreTable_training[["AGEGRP"]]))

#Modèle Logit 

# On crée la base de données sur laquelle on calle le logit
base_logit <- c()
for(i in name_var){base_logit[[i]] <- CoreTable_training[[i]]}

base_logit=as.data.frame(base_logit)
CHECK_NA(base_logit)

base_logit=subset(base_logit,is.na(base_logit$SMOKE)==FALSE)
CHECK_NA(base_logit)

base_logit=base_logit[1:15]
CHECK_NA(base_logit)

base_logit=knnImputation(base_logit, k =10, meth = "median")
CHECK_NA(base_logit)


smoke.glm = glm(SMOKE ~ BMI+RACE_C+WEIGHTBL+ALP+ALT+AST+CA+CREAT+HB+LDH+NEU+PLT+PSA+TBILI, data = base_logit, family = "binomial")
summary(smoke.glm)






######################################################


# Kaplan-Meier Estimates
my.surv_cancer <- Surv(CoreTable_training$LKADT_P,as.numeric(CoreTable_training$DEATH))
my.fit_cancer <-survfit(my.surv_cancer ~ 1)
plot(my.fit_cancer,conf.int=TRUE,mark.time=FALSE,col=4, main="Kaplan-Meier estimate with 95% confidence bounds",xlab="days", ylab="survival function")

# Attention: Il y a deux dates de référence pour la variable LKADT (DATE,RANDOMIZATION).
# Je n'ai pas trouvé ce que ça voulait dire sur le site, j'ai tracé les deux et il ne semble pas y avoir de diff
# En fait, l'une des sous-populations semble avoir moins d'individus mais les fonctions de survie n'ont pas l'air d'être différentes.
#table(CoreTable_training$LKADT_REF)
#my.fit_cancer_check <-survfit(my.surv_cancer ~ CoreTable_training$LKADT_REF)
#plot(my.fit_cancer_check,conf.int=TRUE,mark.time=FALSE,col=2, main="Kaplan-Meier estimate with 95% confidence bounds",xlab="days", ylab="survival function")

# Study Race
tbl <- table(CoreTable_training$RACE_C,useNA = "always")
cbind(tbl,round(prop.table(tbl),2)) # 86% de blancs, 5% de noirs, 3% d'asiatiques (3% missing)

# Study Age
tbl <- table(CoreTable_training$AGEGRP2,useNA = "always")
cbind(tbl,round(prop.table(tbl),2)) # 31% 18-64 ans, 44% 65-74, 24% >=75
# Idée Jb: redécouper la variable âge en quintles pour équilibrer le nombre d'observations par classes

# Study World Region
tbl <- table(CoreTable_training$REGION_C,useNA = "always")
cbind(tbl,round(prop.table(tbl),2)) # E.Europe 13%, W.Europe 29%, N.America 14%, S.America 5%,  (30% missing)

# Two-way table, race and region: Pour regarder comment son distribuées les missing values
tbl = table(CoreTable_training$RACE_C,CoreTable_training$REGION_C)
prop.table(tbl,2) # Les missing values ont l'air d'être à peu près distribué uniformément sur la population

# Gleason Score
tbl <- table(CoreTable_training$GLEAS_DX,useNA = "always")
cbind(tbl,round(prop.table(tbl),2)) # Gros problème sur le Gleason Score, il y a 75% de NA. 
#C'est un pb général aux variables que l'on trouve que dans le ASCENT2

# Study Smoking
tbl <- table(CoreTable_training$SMOKE,useNA = "always")
cbind(tbl,round(prop.table(tbl),2)) # Tjs pb de NA (70%); attention: pour cette variable, les NA sont codées sans valeur.

# Study Tumor-stage Score at Initial diagnosis
tbl <- table(CoreTable_training$TSTAG_DX,useNA = "always")
cbind(tbl,round(prop.table(tbl),2)) # Tjs pb de NA (72%); attention: pour cette variable, les NA sont codées sans valeur.

# Study BMI
summary(CoreTable_training$BMI) # 10 NA, Rappel: >18 = Maigre, 18< & >25 Normal; 25< & >30 Surpoids; > 30 obesité 
d <- density(CoreTable_training$BMI[is.na(CoreTable_training$BMI)==FALSE]) 
plot(d,main="Distribution BMI")
# Il y a cinq variables relatives à la taille/poids, il faudra selectionner la meilleure

# Study of the baseline performance status and Kaplan-Meier 
tbl <- table(CoreTable_training$ECOG_C,useNA = "always")
cbind(tbl,round(prop.table(tbl),2)) 
my.fit_cancer_ECOG <-survfit(my.surv_cancer ~CoreTable_training$ECOG_C)
plot(my.fit_cancer_ECOG,conf.int=FALSE,mark.time=TRUE,col=8, main="Kaplan-Meier estimates with different performance status",xlab="days", ylab="survival function")



