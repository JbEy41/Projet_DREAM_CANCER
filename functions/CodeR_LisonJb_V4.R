#########################################################################################
#                                 DREAM CANCER PROJECT                                  #
#########################################################################################

# Packages
# install.packages("survival") (si nÃ©cessaire)
library(survival)
library(plyr)
library(DMwR)
library(boot)
library(mlogit)
library(VGAM)
library(class)

#Chemin Lison : //paradis/eleves/LGrappin/Bureau/AppGenomique/data
# Autre chemin Lison : C:/Users/Lison/Desktop/Travail 3A/AppGen
# Chemin JB : /Users/eymeoudjean-benoit/Dropbox
CoreTable_training_raw <- read.csv("C:/Users/Lison/Desktop/Travail 3A/AppGen/CoreTable_training.csv",header=TRUE)  
CoreTable_training <- read.csv("C:/Users/Lison/Desktop/Travail 3A/AppGen/CoreTable_training.csv",header=TRUE,na.strings=c("","."),dec=".",stringsAsFactors =FALSE)  


#########################################################################################
# 1 PARTIE DATA

# La base est structurée de la manière suivante:
# var 1:29: description de l'individu (taille,poids,cigarette,placebo,..) 
# var 30:54 données biologiques (num)
# var 55:131 historique médical au sens large (logical)

# On regarde le type des données initiales
# sapply(CoreTable_training_raw,function(x) class(x)) 
# On regarde le type une fois que l'on a transformé
# sapply(CoreTable_training,function(x) class(x)) 

#########################
# 1.1 Recodages variables

# Beaucoup de variables sont Ã  changer en numeric
#sapply(CoreTable_training,function(x) class(x)) 

# Death variable
#table(CoreTable_training$DEATH, useNA = "always") # DonnÃ©es initiales
CoreTable_training$DEATH[is.na(CoreTable_training$DEATH)] <- 0 # On transforme les NA en NO
CoreTable_training$DEATH[CoreTable_training$DEATH=="YES"] <- 1
#table(CoreTable_training$DEATH, useNA = "always") # On regarde que Ã§a marche
CoreTable_training$DEATH=as.numeric(CoreTable_training$DEATH)
#table(CoreTable_training$DEATH)

# Region
CoreTable_training$REGION_C[which(CoreTable_training$REGION_C=="MISSING")] <- NA

# Les variables  WGTBLCAT,HGTBLCAT sont des catÃ©gories des var WEIGHTBL et HEIGHTBL 
# Je les recode en utilisant des quartiles 
CoreTable_training$WGTBLCAT<-cut(CoreTable_training$WEIGHTBL,quantile(subset(CoreTable_training$WEIGHTBL,is.na(CoreTable_training$WEIGHTBL)==FALSE),probs = seq(0, 1, 0.25)), right=TRUE,include.lowest =TRUE, labels=c(1:4))
CoreTable_training$HGTBLCAT<-cut(CoreTable_training$HEIGHTBL,quantile(subset(CoreTable_training$HEIGHTBL,is.na(CoreTable_training$HEIGHTBL)==FALSE),probs = seq(0, 1, 0.25)), right=TRUE,include.lowest =TRUE, labels=c(1:4))

# On recode les indicatrices de consommation de cigarette
CoreTable_training$SMOKE=as.numeric(as.factor(CoreTable_training$SMOKE))-1 # 0: Don't smoke 1:Smke
CoreTable_training$SMOKFREQ=as.numeric(as.factor(CoreTable_training$SMOKFREQ)) # 1:GREATER THAN OR EQUAL TO 1 PACK PER DAY; 2: LESS THAN 1 PACK PER DAY
CoreTable_training$SMOKSTAT=as.numeric(as.factor(CoreTable_training$SMOKSTAT)) # 1: ONGOING 2: WITHIN THE PAST 12 MONTHS

# On recode d'autres variables mal codes
CoreTable_training$ALP=as.numeric(CoreTable_training$ALP)
CoreTable_training$LDH=as.numeric(CoreTable_training$LDH)
CoreTable_training$AGEGRP=as.numeric(CoreTable_training$AGEGRP) 
CoreTable_training$AGEGRP[is.na(CoreTable_training$AGEGRP)] <-85 # On code les 85 ans et +, 85
CoreTable_training$AGEGRP2<-cut(CoreTable_training$AGEGRP,quantile(CoreTable_training$AGEGRP,probs = seq(0, 1, 1/3)), right=TRUE,include.lowest =TRUE, labels=c(1:3))

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

# Table intermédiaire avec les bons types de variables
CoreTable_training_typeOK <- CoreTable_training


### DETECTION DES VALEURS MANQUANTES
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

# Retour à al table initiale avec le bon formatage de variable
CoreTable_training <- CoreTable_training_typeOK


## COMPLETION DES VARIABLES NUMERIQUES A FAIBLE POURCENTAGE DE MISSING
#-----------------------------------------------------------------------
# Liste des variables numériques avec un faible taux de missing
names.faible <- c("BMI", "HEIGHTBL", "WEIGHTBL", "ALP", "ALT", "AST", "CA", "CREAT", "HB", "NEU", "PLT", "PSA", "TBILI", "WBC")
# Fonction de remplacement par la moyenne
replace.by.mean <- function(var.name)
{
  var <- CoreTable_training[,var.name]
  m <- mean(var, na.rm=TRUE)
  var[(is.na(var)==TRUE)] <- m
  CoreTable_training[,var.name] <- var
  return(CoreTable_training)
}
# Remplacement par la moyenne
for (i in names.faible)
{
  print(i)
  CoreTable_training <- replace.by.mean(i)
}
# Check
CHECK_NA(CoreTable_training)


## COMPLETION DES VARIABLES CARACTERE A FAIBLE POURCENTAGE DE MISSING
#---------------------------------------------------------------------
# Seule variable concernée : ECOG_C
table(CoreTable_training$ECOG_C, useNA="always")
# Une seule modalité manquante : on lui affecte la modalité 1 (sur représentée dans l'éch : fully)
CoreTable_training[is.na(CoreTable_training$ECOG_C)==TRUE,] 
CoreTable_training[is.na(CoreTable_training$ECOG_C)==TRUE,"ECOG_C"] <- 0
CHECK_NA(CoreTable_training)


## COMPLETION DES VARIABLES NUMERIQUES A FORT POURCENTAGE DE MISSING
#---------------------------------------------------------------------
# Liste des variables à compléter
names.fort <- c("LDH", "TESTO", "NA.", "MG", "PHOS", "ALB", "TPRO", "CCRC", "GLU")

# Modeles de regression pour analyser le pouvoir explicatif de la variable sur les deux variables qu'on doit prédire
glm.death <- function(var.name, vector.completed=NULL, completed=FALSE)
{
  if (completed==FALSE)
  {
    return(summary(glm(CoreTable_training$DEATH~CoreTable_training[,var.name], family="binomial")))
  } else if (completed==TRUE) {
    return(summary(glm(CoreTable_training$DEATH~vector.completed, family="binomial")))
  }
  
}
lm.lkadt <- function(var.name, vector.completed=NULL, completed=FALSE)
{
  if (completed==FALSE)
  {
    return(summary(lm(CoreTable_training$LKADT_P~CoreTable_training[,var.name])))
  } else if (completed==TRUE) {
    return(summary(lm(CoreTable_training$LKADT_P~vector.completed)))
  }
}
cor.lkadt <- function(var.name, vector.completed=NULL, completed=FALSE)
{
  if (completed==FALSE)
  {
    return(cor(CoreTable_training$DEATH,CoreTable_training[,var.name], use="complete.obs"))
  } else if (completed==TRUE) {
    return(cor(CoreTable_training$DEATH,vector.completed, use="complete.obs"))
  }
}
get.avail <- function(var.name)
{
  return(CoreTable_training[!is.na(CoreTable_training[,var.name]),])
}
get.to.complete <- function(var.name)
{
  return(CoreTable_training[is.na(CoreTable_training[,var.name]),])
}

# Analyse variable par variable

# 1 - LDH
#--------

# Avant complétion
glm.death("LDH") 
# estimateur proche de 0, AIC = 1058.5
lm.lkadt("LDH")
# estimateur proche de 0, R2 = 1.2%
cor.lkadt("LDH")
# corrélation de 0.1289 avec LKADT_P
# Peu de pouvoir explicatif a priori

# Completion 1 : par la moyenne
m.LDH <- replace.by.mean("LDH")[,"LDH"]
glm.death("LDH", vector.completed=m.LDH, completed=TRUE) 
# estrimateur proche de 0, AIC = 2162.5
lm.lkadt("LDH", vector.completed=m.LDH, completed=TRUE)
# estimateur proche de 0, R2 = 0.2%
cor.lkadt("LDH", vector.completed=m.LDH, completed=TRUE)
# corrélation de 0.0867 avec LKADT_P

# Completion 2 : par une régression sur les autres variables biologiques complètes
LDH.avail <- get.avail("LDH")
names.covar <- c("ALP", "ALT", "AST", "CA", "CREAT", "HB", "NEU", "PLT", "PSA", "TBILI", "WBC")
covar.LDH <- as.matrix(LDH.avail[,names.covar])
reg.LDH <- lm(LDH.avail$LDH~covar.LDH)
summary(reg.LDH)
# R2 = 33.7%
# On ne garde que les covariables significatives
names.LDH <- c("ALP", "ALT", "AST", "CREAT", "HB")
covar.LDH <- as.matrix(LDH.avail[,names.LDH])
reg.LDH <- lm(LDH.avail$LDH~covar.LDH)
summary(reg.LDH)
beta.LDH <- summary(reg.LDH)$coefficients[,1]
# R2 = 33.4%
# Completion du vecteur à l'aide de la regression faite ici
LDH.to.complete <- get.to.complete("LDH")
LDH.to.complete[,"LDH"] <- as.matrix(cbind(1,LDH.to.complete[,names.LDH]))%*%t(t(beta.LDH))
r.LDH <- CoreTable_training[,"LDH"]
r.LDH[(is.na(r.LDH)==TRUE)] <- LDH.to.complete[,"LDH"]
# Vérification de l'impact sur death
glm.death("LDH", vector.completed=r.LDH, completed=TRUE) 
# estrimateur proche de 0, AIC = 2146.5
lm.lkadt("LDH", vector.completed=r.LDH, completed=TRUE)
# estimateur proche de 0, R2 = 1.5%
cor.lkadt("LDH", vector.completed=r.LDH, completed=TRUE)
# corrélation de 0.1284 avec LKADT_P

# On complète avec la régression car ne modifie pas trop les rapports entre LDH et DEATH ou LKADT_P
CoreTable_training[,"LDH"] <- r.LDH
CHECK_NA(CoreTable_training)


# 2 - TESTO
#------------


# Avant complétion
glm.death("TESTO") 
# estimateur proche de -0.30, AIC = 988.4
lm.lkadt("TESTO")
# estimateur proche de -10, R2 = 0.7%
cor.lkadt("TESTO")
# corrélation de -0.069 avec LKADT_P
# A un pouvoir explicatif a priori

# Completion 1 : par la moyenne
m.TESTO <- replace.by.mean("TESTO")[,"TESTO"]
glm.death("TESTO", vector.completed=m.TESTO, completed=TRUE) 
# estrimateur proche de -0.62, AIC = 2149.8
lm.lkadt("TESTO", vector.completed=m.TESTO, completed=TRUE)
# estimateur proche de -10, R2 = 0.4%
cor.lkadt("TESTO", vector.completed=m.TESTO, completed=TRUE)
# corrélation de -0.047 avec LKADT_P

# Completion 2 : par une régression sur les autres variables biologiques complètes
TESTO.avail <- get.avail("TESTO")
covar.TESTO <- as.matrix(TESTO.avail[,names.covar])
reg.TESTO <- lm(TESTO.avail$TESTO~covar.TESTO)
summary(reg.TESTO)
# R2 = 0%
# On ne garde que les variables significatives
names.TESTO <- "PLT"
covar.TESTO <- as.matrix(TESTO.avail[,names.TESTO])
reg.TESTO <- lm(TESTO.avail$TESTO~covar.TESTO)
summary(reg.TESTO)
# R2 = 0.1%

# R2 vraiment beaucoup trop faible, on complète par la moyenne car ne distord pas trop les relations
CoreTable_training[,"TESTO"] <- m.TESTO
CHECK_NA(CoreTable_training)


# 3 - NA. (niveau de sodium)
#--------------------------------

# Avant complétion
glm.death("NA.") 
# estimateur proche de 0, AIC = 1544.5
lm.lkadt("NA.")
# estimateur proche de 8, R2 = 0.7%
cor.lkadt("NA.")
# corrélation de -0.075 avec LKADT_P
# Peu de pouvoir explicatif a priori

# Completion 1 : par la moyenne
m.NA <- replace.by.mean("NA.")[,"NA."]
glm.death("NA.", vector.completed=m.NA, completed=TRUE) 
# estrimateur proche de 0, AIC = 2168.4
lm.lkadt("NA.", vector.completed=m.NA, completed=TRUE)
# estimateur proche de 8, R2 = 0.67%
cor.lkadt("NA.", vector.completed=m.NA, completed=TRUE)
# corrélation de -0.063 avec LKADT_P

# Completion 2 : par une régression sur les autres variables biologiques complètes
NA.avail <- get.avail("NA.")
covar.NA <- as.matrix(NA.avail[,names.covar])
reg.NA <- lm(NA.avail[,"NA."]~covar.NA)
summary(reg.NA)
# R2 = 12.37%
# On ne garde que les variables significatives
names.NA <- c("CA","HB","NEU")
covar.NA <- as.matrix(NA.avail[,names.NA])
reg.NA <- lm(NA.avail[,"NA."]~covar.NA)
summary(reg.NA)
beta.NA <- summary(reg.NA)$coefficients[,1]
# R2 = 11.9%
# Completion du vecteur à l'aide de la regression faite ici
NA.to.complete <- get.to.complete("NA.")
NA.to.complete[,"NA."] <- as.matrix(cbind(1,NA.to.complete[,names.NA]))%*%t(t(beta.NA))
r.NA <- CoreTable_training[,"NA."]
r.NA[(is.na(r.NA)==TRUE)] <- NA.to.complete[,"NA."]
# Vérification de l'impact sur death
glm.death("NA.", vector.completed=r.NA, completed=TRUE) 
# estrimateur proche de 0, AIC = 2161
lm.lkadt("NA.", vector.completed=r.NA, completed=TRUE)
# estimateur proche de 8, R2 = 0.68%
cor.lkadt("NA.", vector.completed=r.NA, completed=TRUE)
# corrélation de -0.092 avec LKADT_P

# On complète par la moyenne car R2 trop faible
CoreTable_training[,"NA."] <- m.NA
CHECK_NA(CoreTable_training)


# 4 - MG
#--------

# Avant complétion
glm.death("MG") 
# estimateur proche de 0.2, AIC = 1506.8
lm.lkadt("MG")
# estimateur proche de 56, R2 = 0.24%
cor.lkadt("MG")
# corrélation de 0.031 avec LKADT_P
# A un pouvoir explicatif a priori

# Completion 1 : par la moyenne
m.MG <- replace.by.mean("MG")[,"MG"]
glm.death("MG", vector.completed=m.MG, completed=TRUE) 
# estimateur proche de 0.2, AIC = 2173.8
lm.lkadt("MG", vector.completed=m.MG, completed=TRUE)
# estimateur proche de 56, R2 = 0.2%
cor.lkadt("MG", vector.completed=m.MG, completed=TRUE)
# corrélation de 0.026 avec LKADT_P

# Completion 2 : par une régression sur les autres variables biologiques complètes
MG.avail <- get.avail("MG")
covar.MG <- as.matrix(MG.avail[,names.covar])
reg.MG <- lm(MG.avail[,"MG"]~covar.MG)
summary(reg.MG)
# R2 = 0.3%
# On ne garde que les variables significatives
names.MG <- c("TBILI")
covar.MG <- as.matrix(MG.avail[,names.MG])
reg.MG <- lm(MG.avail[,"MG"]~covar.MG)
summary(reg.MG)
# R2 = 1%

# R2 beaucoup trop faible, on complète par la moyenne
CoreTable_training[,"MG"] <- m.MG
CHECK_NA(CoreTable_training)


# 5 - PHOS
#----------

# Avant complétion
glm.death("PHOS") 
# estimateur proche de 0.1, AIC = 1516.2
lm.lkadt("PHOS")
# estimateur proche de 117, R2 = 1.5%
cor.lkadt("PHOS")
# corrélation de 0.02 avec LKADT_P
# Peut avoir un pouvoir explicatif a priori

# Completion 1 : par la moyenne
m.PHOS <- replace.by.mean("PHOS")[,"PHOS"]
glm.death("PHOS", vector.completed=m.PHOS, completed=TRUE) 
# estimateur proche de 0.1, AIC = 2174.5
lm.lkadt("PHOS", vector.completed=m.PHOS, completed=TRUE)
# estimateur proche de 117, R2 = 0.13%
cor.lkadt("PHOS", vector.completed=m.PHOS, completed=TRUE)
# corrélation de 0.017 avec LKADT_P

# Completion 2 : par une régression sur les autres variables biologiques complètes
PHOS.avail <- get.avail("PHOS")
covar.PHOS <- as.matrix(PHOS.avail[,names.covar])
reg.PHOS <- lm(PHOS.avail[,"PHOS"]~covar.PHOS)
summary(reg.PHOS)
# R2 = 7.8%
# On ne garde que les variables significatives
names.PHOS <- c("ALP", "CA", "NEU", "WBC")
covar.PHOS <- as.matrix(PHOS.avail[,names.PHOS])
reg.PHOS <- lm(PHOS.avail[,"PHOS"]~covar.PHOS)
summary(reg.PHOS)
# R2 = 6.9%

# R2 trop faible, on complète par la moyenne car ne distord pas trop
CoreTable_training[,"PHOS"] <- m.PHOS
CHECK_NA(CoreTable_training)


# 6 - ALB
#----------

# Avant complétion
glm.death("ALB") 
# estimateur proche de -0.1, AIC = 1469.8
lm.lkadt("ALB")
# estimateur proche de -1.5, R2 = 0%
cor.lkadt("ALB")
# corrélation de -0.23 avec LKADT_P
# Pas de pouvoir explicatif a priori

# Completion 1 : par la moyenne
m.ALB <- replace.by.mean("ALB")[,"ALB"]
glm.death("ALB", vector.completed=m.ALB, completed=TRUE) 
# estimateur proche de -0.1, AIC = 2111
lm.lkadt("ALB", vector.completed=m.ALB, completed=TRUE)
# estimateur proche de -1.5, R2 = 0%
cor.lkadt("ALB", vector.completed=m.ALB, completed=TRUE)
# corrélation de -0.19 avec LKADT_P

# Completion 2 : par une régression sur les autres variables biologiques complètes
ALB.avail <- get.avail("ALB")
covar.ALB <- as.matrix(ALB.avail[,names.covar])
reg.ALB <- lm(ALB.avail[,"ALB"]~covar.ALB)
summary(reg.ALB)
# R2 = 24.5%
# On ne garde que les variables significatives
names.ALB <- c("CA", "HB", "PLT", "TBILI")
covar.ALB <- as.matrix(ALB.avail[,names.ALB])
reg.ALB <- lm(ALB.avail[,"ALB"]~covar.ALB)
summary(reg.ALB)
beta.ALB <- summary(reg.ALB)$coefficients[,1]
# R2 = 24.7%

# Completion du vecteur à l'aide de la regression faite ici
ALB.to.complete <- get.to.complete("ALB")
ALB.to.complete[,"ALB"] <- as.matrix(cbind(1,ALB.to.complete[,names.ALB]))%*%t(t(beta.ALB))
r.ALB <- CoreTable_training[,"ALB"]
r.ALB[(is.na(r.ALB)==TRUE)] <- ALB.to.complete[,"ALB"]
# Vérification de l'impact sur death
glm.death("ALB", vector.completed=r.ALB, completed=TRUE) 
# estrimateur proche de -0.13, AIC = 2089
lm.lkadt("ALB", vector.completed=r.ALB, completed=TRUE)
# estimateur proche de -1.17, R2 = 0%
cor.lkadt("ALB", vector.completed=r.ALB, completed=TRUE)
# corrélation de -0.23 avec LKADT_P

# On complète par la régression car meilleure que complétion par la moyenne et ne distord pas trop le modèle
CoreTable_training[,"ALB"] <- r.ALB
CHECK_NA(CoreTable_training)

# 7 - TPRO
#-----------

# Avant complétion
glm.death("TPRO") 
# estimateur proche de 0.03, AIC = 1506.6
lm.lkadt("TPRO")
# estimateur proche de 7.7, R2 = 2%
cor.lkadt("TPRO")
# corrélation de 0.095 avec LKADT_P
# Pas de pouvoir explicatif a priori

# Completion 1 : par la moyenne
m.TPRO <- replace.by.mean("TPRO")[,"TPRO"]
glm.death("TPRO", vector.completed=m.TPRO, completed=TRUE) 
# estimateur proche de 0.03, AIC = 2164.711
lm.lkadt("TPRO", vector.completed=m.TPRO, completed=TRUE)
# estimateur proche de 7.7, R2 = 1.7%
cor.lkadt("TPRO", vector.completed=m.TPRO, completed=TRUE)
# corrélation de 0.08 avec LKADT_P

# Completion 2 : par une régression sur les autres variables biologiques complètes
TPRO.avail <- get.avail("TPRO")
covar.TPRO <- as.matrix(TPRO.avail[,names.covar])
reg.TPRO <- lm(TPRO.avail[,"TPRO"]~covar.TPRO)
summary(reg.TPRO)
# R2 = 16.5%
# On ne garde que les variables significatives
names.TPRO <- c("CA", "HB", "CREAT", "NEU", "WBC")
covar.TPRO <- as.matrix(TPRO.avail[,names.TPRO])
reg.TPRO <- lm(TPRO.avail[,"TPRO"]~covar.TPRO)
summary(reg.TPRO)
beta.TPRO <- summary(reg.TPRO)$coefficients[,1]
# R2 = 15.3%

# Completion du vecteur à l'aide de la regression faite ici
TPRO.to.complete <- get.to.complete("TPRO")
TPRO.to.complete[,"TPRO"] <- as.matrix(cbind(1,TPRO.to.complete[,names.TPRO]))%*%t(t(beta.TPRO))
r.TPRO <- CoreTable_training[,"TPRO"]
r.TPRO[(is.na(r.TPRO)==TRUE)] <- TPRO.to.complete[,"TPRO"]
# Vérification de l'impact sur death
glm.death("TPRO", vector.completed=r.TPRO, completed=TRUE) 
# estrimateur proche de 0.02, AIC = 2170.8
lm.lkadt("TPRO", vector.completed=r.TPRO, completed=TRUE)
# estimateur proche de 6.9, R2 = 1.4%
cor.lkadt("TPRO", vector.completed=r.TPRO, completed=TRUE)
# corrélation de 0.05 avec LKADT_P

# On complète avec la moyenne car distord moins la relation entre death et TPRO non completee
CoreTable_training[,"TPRO"] <- m.TPRO
CHECK_NA(CoreTable_training)


# 8 - CCRC
#----------

# Avant complétion
glm.death("CCRC") 
# estimateur proche de 0, AIC = 1282.6
lm.lkadt("CCRC")
# estimateur proche de 0.8, R2 = 0.3%
cor.lkadt("CCRC")
# corrélation de 0.03 avec LKADT_P
# Pas de pouvoir explicatif a priori

# Completion 1 : par la moyenne
m.CCRC <- replace.by.mean("CCRC")[,"CCRC"]
glm.death("CCRC", vector.completed=m.CCRC, completed=TRUE) 
# estimateur proche de 0, AIC = 2173.9
lm.lkadt("CCRC", vector.completed=m.CCRC, completed=TRUE)
# estimateur proche de 0.8, R2 = 0.3%
cor.lkadt("CCRC", vector.completed=m.CCRC, completed=TRUE)
# corrélation de 0.024 avec LKADT_P

# Completion 2 : par une régression sur les autres variables biologiques complètes
CCRC.avail <- get.avail("CCRC")
covar.CCRC <- as.matrix(CCRC.avail[,names.covar])
reg.CCRC <- lm(CCRC.avail[,"CCRC"]~covar.CCRC)
summary(reg.CCRC)
# R2 = 49.6% : GOOD !
# On ne garde que les variables significatives
names.CCRC <- c("ALT", "CREAT")
covar.CCRC <- as.matrix(CCRC.avail[,names.CCRC])
reg.CCRC <- lm(CCRC.avail[,"CCRC"]~covar.CCRC)
summary(reg.CCRC)
beta.CCRC <- summary(reg.CCRC)$coefficients[,1]
# R2 = 49.5% : GOOD !

# Completion du vecteur à l'aide de la regression faite ici
CCRC.to.complete <- get.to.complete("CCRC")
CCRC.to.complete[,"CCRC"] <- as.matrix(cbind(1,CCRC.to.complete[,names.CCRC]))%*%t(t(beta.CCRC))
r.CCRC <- CoreTable_training[,"CCRC"]
r.CCRC[(is.na(r.CCRC)==TRUE)] <- CCRC.to.complete[,"CCRC"]
# Vérification de l'impact sur death
glm.death("CCRC", vector.completed=r.CCRC, completed=TRUE) 
# estrimateur proche de 0, AIC = 2173.7
lm.lkadt("CCRC", vector.completed=r.CCRC, completed=TRUE)
# estimateur proche de 0.8, R2 = 0.6%
cor.lkadt("CCRC", vector.completed=r.CCRC, completed=TRUE)
# corrélation de 0.027 avec LKADT_P

# On compléte par la régression parce que bon R2 et ne distord pas les relations
CoreTable_training[,"CCRC"] <- r.CCRC
CHECK_NA(CoreTable_training)


# 9 - GLU
#---------

# Avant complétion
glm.death("GLU") 
# estimateur proche de 0, AIC = 1524.5
lm.lkadt("GLU")
# estimateur proche de -14, R2 = 1.5%
cor.lkadt("GLU")
# corrélation de -0.1 avec LKADT_P
# Peut avoir un pouvoir explicatif a priori

# Completion 1 : par la moyenne
m.GLU <- replace.by.mean("GLU")[,"GLU"]
glm.death("GLU", vector.completed=m.GLU, completed=TRUE) 
# estimateur proche de 0, AIC = 2163.1
lm.lkadt("GLU", vector.completed=m.GLU, completed=TRUE)
# estimateur proche de -14, R2 = 1.2%
cor.lkadt("GLU", vector.completed=m.GLU, completed=TRUE)
# corrélation de -0.08 avec LKADT_P

# Completion 2 : par une régression sur les autres variables biologiques complètes
GLU.avail <- get.avail("GLU")
covar.GLU <- as.matrix(GLU.avail[,names.covar])
reg.GLU <- lm(GLU.avail[,"GLU"]~covar.GLU)
summary(reg.GLU)
# R2 = 11.8%
# On ne garde que les variables significatives
names.GLU <- c("NEU", "WBC", "ALT")
covar.GLU <- as.matrix(GLU.avail[,names.GLU])
reg.GLU <- lm(GLU.avail[,"GLU"]~covar.GLU)
summary(reg.GLU)
# R2 = 10.9%

# R2 trop faible, on complète par la moyenne (ne distord pas trop)
CoreTable_training[,"GLU"] <- m.GLU
CHECK_NA(CoreTable_training)


## RESTE A COMPLETER REGION_C ET SMOKE
# On les complète par kNN

# 1 - SMOKE
#-----------
table(CoreTable_training[,"SMOKE"], CoreTable_training[,"STUDYID"])
# Variable qui n'est renseignée que pour l'étude ASCENT2

# On voudrait compléter par kNN sur les variables numériques (variables biologiques)
glm.smoke <- glm(CoreTable_training$SMOKE~as.matrix(CoreTable_training[,names.faible]), family="binomial")
summary(glm.smoke)
# PB : aucune des variables biologiques n'est significative !!

# TEST KNN qui fonctionne meme si on ne peut pas vraiment l'appliquer ...
kNN.smoke.training <- CoreTable_training[!is.na(CoreTable_training[,"SMOKE"]),c(names.faible)]
kNN.smoke.test <- CoreTable_training[is.na(CoreTable_training[,"SMOKE"]),c(names.faible)]
classify <- as.factor(CoreTable_training[!is.na(CoreTable_training[,"SMOKE"]),"SMOKE"])
kNN.smoke.predict <- knn(kNN.smoke.training, kNN.smoke.test, classify, k=5)

# PB : Comment compléter cette variable ?

# REGION_C
#----------
table(CoreTable_training[,"REGION_C"], CoreTable_training[,"STUDYID"])
# Variable qui n'est pas renseignée pour l'étude ASCENT2
# PAs d'idées pour compléter


