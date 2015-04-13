CoreTable <- read.csv("//paradis/eleves/LGrappin/Bureau/AppGenomique/data/CoreTable_training.csv",header=TRUE,na.strings=c("","."),dec=".",stringsAsFactors =FALSE)  
CoreTable <- as.data.frame(CoreTable)
attach(CoreTable)

CoreTable.raw <- read.csv("//paradis/eleves/LGrappin/Bureau/AppGenomique/data/CoreTable_training.csv")  

# Change modalities
# Sortie des modalités
pct.missing <- function(var)
{
  return(100*sum(is.na(var))/length(var))
}
pct.missing(DEATH)

noquote(paste("CoreTable$",names(CoreTable)[5],sep=""))
