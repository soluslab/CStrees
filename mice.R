setwd("/Users/eduarte/Documents/gitHubRepos/CStrees/")
setwd("/Users/liamsolus/Dropbox/Research/github-repositories/CStrees/")


# A function for producing fitted interventional CStrees.  The function returns a list of lists in which 
# the first element is the interventionalCStree, the second is its penalized log-likelihood given the data,
# where the penalization term is adjusted to account for the magnitude of the causal effect of each intervention
# target learned, and the third is the sum of the magnitudes of these causal effects.  
# currently entries 4, 5 6, and 7 also record the log-likelihood of the interventional CStree, the log-likelihood 
# of its observational subtree, the log-likelihood of its interventional subtree and its BIC, respectively. 
fittedInterventionalCStrees <- function(p,d,L,S,Dat) {
  Trees <- list()
  Targets <- list(list(),list("1"))
  s <- length(S)
  if (p!=2) {
    for (i in c(2:s)) {
      P <- list()
      LevelTargs <- list(list())
      for (j in c(1:length(S[[i]]))) {
        LevelTargs <- c(LevelTargs,combn(S[[i]],j,simplify=FALSE))
      }
      for (I in Targets) {
        for (M in LevelTargs) {
          if (i!=2) {
            P<-c(P,list(c(I,list(M))))
          } else {
            P<-c(P,list(list(I,M)))
          }
        }
        Targets <- P
      }
    }
  }
  for (I in Targets) {
    ICST <- sevt_fit(toInterventionalCStree(p,2,L,I),Dat,lambda=1)
    obsICST <- subtree(ICST,c("Obs"))
    intICST <- subtree(ICST,c("Int"))
    stagesByLevel<- lapply(stages(ICST),unique)
    stagesByLevel<- lapply(stagesByLevel, length)
    numStages<- Reduce("+",stagesByLevel) +1
    CE <- c()
    for (j in c(2:(p+1))) {
      U <- as.list(unique(stages(ICST,paste("V",j,sep = ""))))
      V <- U[as.integer(U)<=2^(j-2)]
      stageeffects <-c()
      count <- 1
      for (i in V) {
        if (is.element(toString((as.integer(i)+(2^(j-2)))),U) == TRUE) {
          count <- count+1
          probObs <- eval(parse(text = paste("obsICST$prob$V",j,"$`",i,"`",sep = "")))
          probInt <- eval(parse(text = paste("intICST$prob$V",j,"$`",toString((as.integer(i)+(2^(j-2)))),"`",sep = "")))
          ce <- abs((probInt[1]+2*probInt[2])-(probObs[1]+2*probObs[2]))
          stageeffects <-c(stageeffects,ce)
        }
        CE <- c(CE,sum(stageeffects))
        #      CE <- c(CE,sum(stageeffects)/(length(eval(parse(text = paste("obsICST$stages$V",j,sep = ""))))))
      }
    }
    causalEffect <- sum(CE)
    s <- (logLik(obsICST) + logLik(intICST))*2-(numStages-(causalEffect))*log(nrow(Dat))
    b <- (logLik(obsICST) + logLik(intICST))*2-(numStages)*log(nrow(Dat))
    Trees <- c(Trees,list(list(list(ICST),list(s),list((logLik(obsICST) + logLik(intICST))),list(causalEffect),list(logLik(obsICST)),list((logLik(intICST))),list(b))))
  }
  return(Trees)
}


#---pCAMKII pPKCG NR1 pS6 analysis-----
# Mice data: importing data set from UCI Machine Learning Repository.
library(gdata)
myCortexdf <- read.xls("Data_Cortex_Nuclear.xls", sheet = 1, header = TRUE)
summary(myCortexdf)

# We select four variables according to column two of Table 3, which are noted to be targeted by memantine. 
myObsCortex<- subset(myCortexdf,myCortexdf$class == "c-SC-s", select = c(pCAMKII_N,pPKCG_N,NR1_N,pS6_N,class))
summary(myObsCortex)
nrow(myObsCortex)

# Then we discretize by quantile
myObsCortex$pCAMKII_Ncat <- cut(myObsCortex$pCAMKII_N,quantile(myObsCortex$pCAMKII_N, probs = c(0,0.5,1)),include.lowest =TRUE)
myObsCortex$pPKCG_Ncat <- cut(myObsCortex$pPKCG_N,quantile(myObsCortex$pPKCG_N, probs = c(0,0.5,1)),include.lowest =TRUE)
myObsCortex$NR1_Ncat  <- cut(myObsCortex$NR1_N,quantile(myObsCortex$NR1_N, probs = c(0,0.5,1)),include.lowest =TRUE)
myObsCortex$pS6_Ncat  <- cut(myObsCortex$pS6_N,quantile(myObsCortex$pS6_N, probs = c(0,0.5,1)),include.lowest =TRUE)

# After discretization, we learn the BIC-optimal CStree on the observed data; 
# that is, the samples with class c-SC-s.  The interventional data will be those with class c-SC-m.
myObsdata <- data.frame(myObsCortex$pCAMKII_Ncat,myObsCortex$pPKCG_Ncat,myObsCortex$NR1_Ncat,myObsCortex$pS6_Ncat)                    
names(myObsdata) <- c("V1","V2","V3","V4") #This relabels the variables in the dataset with the variable names produced by CStrees(p,2).
levels(myObsdata$V1) <- c(1:2) #This relabels the outcomes to match the outcome names used by CStrees(p,2).
levels(myObsdata$V2) <- c(1:2)
levels(myObsdata$V3) <- c(1:2)
levels(myObsdata$V4) <- c(1:2)
summary(myObsdata)
# Generate all CStrees on four binary variables to be scored.
CStreesList <- CStrees(4,2)
# Score all CStrees on four binary varibales given the data myObsdata.
M <- CStreesList[[100]]
plot(M)
M.fit <- sevt_fit(M,data=myObsdata,lambda=1)
t <- BIC(M.fit)
MM <- M.fit
for (N in CStreesList) {
  N.fit <- sevt_fit(N,data=myObsdata,lambda=1)
  s <- BIC(N.fit)
  if (t > s) {
    MM <- N.fit
    t <- s
  }
}
# The BIC-optimal CStree given the data myObsdata is then:
optObsTree <- MM
plot(optObsTree)
summary(optObsTree)
# Its BIC score is 
obsScore <- t
obsScore
# Its log-likelihood is
obsLogLik <- logLik(optObsTree)
obsLogLik


#-------------- Tree #1 (optObsTree)--------------------------------------------
## Causal ordering = 1 4 2 3  
# We first build the observational and interventional data set where class = c-SC-s is the observed condition and
# class = c-SC-m is the intervention.  Note that we already have the discretized observational data set, but need to 
# record its class as well:
myObsdata <- data.frame(myObsCortex$class,myObsCortex$pCAMKII_Ncat,myObsCortex$pPKCG_Ncat,myObsCortex$NR1_Ncat,myObsCortex$pS6_Ncat)
names(myObsdata) <- c("V1","V2","V4","V5","V3")
levels(myObsdata$V2) <- c(1:2) #This relabels the outcomes to match the outcome names used by CStrees(p,2).
levels(myObsdata$V3) <- c(1:2)
levels(myObsdata$V4) <- c(1:2)
levels(myObsdata$V5) <- c(1:2) 
summary(myObsdata)
#Then we need to build the interventional discretized data set.
myIntCortex<- subset(myCortexdf,myCortexdf$class == "c-SC-m", select = c(pCAMKII_N,pPKCG_N,NR1_N,pS6_N,class))
summary(myIntCortex)
myIntCortex$pCAMKII_Ncat <- cut(myIntCortex$pCAMKII_N,quantile(myIntCortex$pCAMKII_N, probs = c(0,0.5,1)),include.lowest =TRUE)
myIntCortex$pPKCG_Ncat <- cut(myIntCortex$pPKCG_N,quantile(myIntCortex$pPKCG_N, probs = c(0,0.5,1)),include.lowest =TRUE)
myIntCortex$NR1_Ncat  <- cut(myIntCortex$NR1_N,quantile(myIntCortex$NR1_N, probs = c(0,0.5,1)),include.lowest =TRUE)
myIntCortex$pS6_Ncat  <- cut(myIntCortex$pS6_N,quantile(myIntCortex$pS6_N, probs = c(0,0.5,1)),include.lowest =TRUE)
myIntdata <- data.frame(myIntCortex$class,myIntCortex$pCAMKII_Ncat,myIntCortex$pPKCG_Ncat,myIntCortex$NR1_Ncat,myIntCortex$pS6_Ncat)                    
names(myIntdata) <- c("V1","V2","V4","V5","V3") #This relabels the variables in the dataset with the variable names produced by CStrees(p,2).
levels(myIntdata$V2) <- c(1:2) #This relabels the outcomes to match the outcome names used by CStrees(p,2).
levels(myIntdata$V3) <- c(1:2)
levels(myIntdata$V4) <- c(1:2)
levels(myIntdata$V5) <- c(1:2) 
summary(myIntdata)
# We then combine the discretizeed observed and interventional data sets into a single data frame.
myMixdata <- rbind(myObsdata,myIntdata)
myMixdata$V1 <- factor(myMixdata$V1)
levels(myMixdata$V1) <- factor(c("Int","Obs"))
myMixdata$V1<-relevel(myMixdata$V1,"Obs")
summary(myMixdata)

# We now construct the stage partition of the observed tree which is needed as input for fittedInterventionalCStrees.
mecParts<- list(list(list(list("1","2")),
                     list(list("3","4")),
                     list(list("1","2","3","4"),list("5","7"))
))

# Producing the different fitted interventional CStrees possible from intervening on the tree with staging mecParts.
listInterventionalCStrees<-list()
for (part in mecParts){
  ObsTree<-toCStree(4,2,part)
  Smore <- stages(ObsTree)
  S<-lapply(Smore, unique)
  S<-c(list(list(1)),S)
  newInterventionalTrees <-fittedInterventionalCStrees(4,2,part,S,myMixdata)
  listInterventionalCStrees<-c(listInterventionalCStrees,newInterventionalTrees)
}

# Find the fitted interventional CStree with the highest BIC score.
MM <- listInterventionalCStrees[[1]]
plot(MM[[1]][[1]])
t <- listInterventionalCStrees[[1]][[7]][[1]][1]
for (N in listInterventionalCStrees) {
  s <- N[[7]][[1]][1]
  if (s > t) {
    MM <- N
    t <- s
  }
}

# The optimal scoring interventional CStree is
optInt1 <- MM
plot(optInt1[[1]][[1]])
# Its score is
optIntScore1<-t
optIntScore1
# We begin a list containing the optimal trees learned from each member of the equvialence class
optIntTrees <- list(list(optInt1))


#-------------- Tree #2 --------------------------------------------
## Causal ordering = 4 1 2 3  
# We first build the observational and interventional data set where class = c-SC-s is the observed condition and
# class = c-SC-m is the intervention.  Note that we already have the discretized observational data set, but need to 
# record its class as well:
myObsdata <- data.frame(myObsCortex$class,myObsCortex$pCAMKII_Ncat,myObsCortex$pPKCG_Ncat,myObsCortex$NR1_Ncat,myObsCortex$pS6_Ncat)
names(myObsdata) <- c("V1","V3","V4","V5","V2")
levels(myObsdata$V2) <- c(1:2) #This relabels the outcomes to match the outcome names used by CStrees(p,2).
levels(myObsdata$V3) <- c(1:2)
levels(myObsdata$V4) <- c(1:2)
levels(myObsdata$V5) <- c(1:2) 
summary(myObsdata)
#Then we need to build the interventional discretized data set.
myIntCortex<- subset(myCortexdf,myCortexdf$class == "c-SC-m", select = c(pCAMKII_N,pPKCG_N,NR1_N,pS6_N,class))
summary(myIntCortex)
myIntCortex$pCAMKII_Ncat <- cut(myIntCortex$pCAMKII_N,quantile(myIntCortex$pCAMKII_N, probs = c(0,0.5,1)),include.lowest =TRUE)
myIntCortex$pPKCG_Ncat <- cut(myIntCortex$pPKCG_N,quantile(myIntCortex$pPKCG_N, probs = c(0,0.5,1)),include.lowest =TRUE)
myIntCortex$NR1_Ncat  <- cut(myIntCortex$NR1_N,quantile(myIntCortex$NR1_N, probs = c(0,0.5,1)),include.lowest =TRUE)
myIntCortex$pS6_Ncat  <- cut(myIntCortex$pS6_N,quantile(myIntCortex$pS6_N, probs = c(0,0.5,1)),include.lowest =TRUE)
myIntdata <- data.frame(myIntCortex$class,myIntCortex$pCAMKII_Ncat,myIntCortex$pPKCG_Ncat,myIntCortex$NR1_Ncat,myIntCortex$pS6_Ncat)                    
names(myIntdata) <- c("V1","V3","V4","V5","V2") #This relabels the variables in the dataset with the variable names produced by CStrees(p,2).
levels(myIntdata$V2) <- c(1:2) #This relabels the outcomes to match the outcome names used by CStrees(p,2).
levels(myIntdata$V3) <- c(1:2)
levels(myIntdata$V4) <- c(1:2)
levels(myIntdata$V5) <- c(1:2) 
summary(myIntdata)
# We then combine the discretizeed observed and interventional data sets into a single data frame.
myMixdata <- rbind(myObsdata,myIntdata)
myMixdata$V1 <- factor(myMixdata$V1)
levels(myMixdata$V1) <- factor(c("Int","Obs"))
myMixdata$V1<-relevel(myMixdata$V1,"Obs")
summary(myMixdata)

# We now construct the stage partition of the observed tree which is needed as input for fittedInterventionalCStrees.
mecParts<- list(list(list(list("1","2")),
                     list(list("2","4")),
                     list(list("1","2","5","6"),list("3","7"))
))

# Producing the different fitted interventional CStrees possible from intervening on the tree with staging mecParts.
listInterventionalCStrees<-list()
for (part in mecParts){
  ObsTree<-toCStree(4,2,part)
  Smore <- stages(ObsTree)
  S<-lapply(Smore, unique)
  S<-c(list(list(1)),S)
  newInterventionalTrees <-fittedInterventionalCStrees(4,2,part,S,myMixdata)
  listInterventionalCStrees<-c(listInterventionalCStrees,newInterventionalTrees)
}

# Find the fitted interventional CStree with the highest BIC score.
MM <- listInterventionalCStrees[[1]]
plot(MM[[1]][[1]])
t <- listInterventionalCStrees[[1]][[7]][[1]][1]
for (N in listInterventionalCStrees) {
  s <- N[[7]][[1]][1]
  if (s > t) {
    MM <- N
    t <- s
  }
}

# The optimal scoring interventional CStree is
optInt2 <- MM
plot(optInt2[[1]][[1]])
# Its score is
optIntScore2<-t
optIntScore2
# Notice that optIntScore1 = optIntScore2, which agrees with the observation that the resulting 
# interventional CStrees are statistically equivalent:
optIntScore1 == optIntScore2





#We now do the same analysis for a second set of four proteins from the same data set.
# Mice data: importing data set from UCI Machine Learning Repository.
library(gdata)
myCortexdf <- read.xls("Data_Cortex_Nuclear.xls", sheet = 1, header = TRUE)
summary(myCortexdf)

# We select four variables according to column two of Table 3, which are noted to be targeted by memantine. 
myObsCortex<- subset(myCortexdf,myCortexdf$class == "c-SC-s", select = c(pPKCG_N,pNUMB_N,pNR1_N,pCAMKII_N,class))
summary(myObsCortex)
nrow(myObsCortex)

# Then we discretize by quantile:
myObsCortex$pPKCG_Ncat <- cut(myObsCortex$pPKCG_N,quantile(myObsCortex$pPKCG_N, probs = c(0,0.5,1)),include.lowest =TRUE)
myObsCortex$pNUMB_Ncat <- cut(myObsCortex$pNUMB_N,quantile(myObsCortex$pNUMB_N, probs = c(0,0.5,1)),include.lowest =TRUE)
myObsCortex$pNR1_Ncat  <- cut(myObsCortex$pNR1_N,quantile(myObsCortex$pNR1_N, probs = c(0,0.5,1)),include.lowest =TRUE)
myObsCortex$pCAMKII_Ncat  <- cut(myObsCortex$pCAMKII_N,quantile(myObsCortex$pCAMKII_N, probs = c(0,0.5,1)),include.lowest =TRUE)

# After discretization, we learn the BIC-optimal CStree for the observed data; that is, the samples with 
# class c-SC-s.  The interventional data will be those with class c-SC-m.
myObsdata <- data.frame(myObsCortex$pPKCG_Ncat,myObsCortex$pNUMB_Ncat,myObsCortex$pNR1_Ncat,myObsCortex$pCAMKII_Ncat)                    
names(myObsdata) <- c("V1","V2","V3","V4") # This relabels the variables in the dataset with the variable names produced by CStrees(p,2).
levels(myObsdata$V1) <- c(1:2) # This relabels the outcomes to match the outcome names used by CStrees(p,2).
levels(myObsdata$V2) <- c(1:2)
levels(myObsdata$V3) <- c(1:2)
levels(myObsdata$V4) <- c(1:2)
summary(myObsdata)
# Generate all CStrees on four binary variables to be scored.
CStreesList <- CStrees(4,2)
# Score all CStrees on four binary varibales given the data myObsdata.
M <- CStreesList[[100]]
plot(M)
M.fit <- sevt_fit(M,data=myObsdata,lambda=1)
t <- BIC(M.fit)
MM <- M.fit
for (N in CStreesList) {
  N.fit <- sevt_fit(N,data=myObsdata,lambda=1)
  s <- BIC(N.fit)
  if (t > s) {
    MM <- N.fit
    t <- s
  }
}
# The BIC-optimal CStree given the data myObsdata is then:
optObsTree <- MM
plot(optObsTree)
summary(optObsTree)
# Its BIC score is 
obsScore <- t
obsScore
# Its log-likelihood is
obsLogLik <- logLik(optObsTree)
obsLogLik

#----Learning an optimal interventional CStree for the data (c-SC-s,c-SC-m)------------
# We then learn an optimal element of the equivalence class of optObsTree.  To do so, we score each
# interventional CStree that arises from any intervention on an element of the equivalence class of optObsTree, 
# where the score is the penalized log-likelihood is the second entry of the output of the fittedInterventionalCStrees
#function above.  This score tries to balance between learning a model that optimizes the
# likelihood of the data while minimizing the number of parameters so that the new parameters introduced to
# capture the intervention targets should have as large of a causal effect as possible. This score is the second 
# element in each list produced by the function fittedInterventionalCStrees. 
# There are four CStrees in the equivalence class of optObsTree, each distingushed by its causal ordering.
# We now learn the optimal interventional CStree for each class, and take the highest scoring of the four
# learned trees.  Altrnatively, we can run the same code with BIC, which is the 7-th entry in the lists output by
# fittedInterventionalCStrees.  

#-------------- Tree #1 (optObsTree)--------------------------------------------
## Causal ordering = 1 2 3 4 
# We first build the observational and interventional data set where class = c-SC-s is the observed condition and
# class = c-SC-m is the intervention.  Note that we already have the discretized observational data set, but need to 
# record its class as well:
myObsdata <- data.frame(myObsCortex$class,myObsCortex$pPKCG_Ncat,myObsCortex$pNUMB_Ncat,myObsCortex$pNR1_Ncat,myObsCortex$pCAMKII_Ncat)
names(myObsdata) <- c("V1","V2","V3","V4","V5")
levels(myObsdata$V2) <- c(1:2) #This relabels the outcomes to match the outcome names used by CStrees(p,2).
levels(myObsdata$V3) <- c(1:2)
levels(myObsdata$V4) <- c(1:2)
levels(myObsdata$V5) <- c(1:2) 
summary(myObsdata)
#Then we need to build the interventional discretized data set.
myIntCortex<- subset(myCortexdf,myCortexdf$class == "c-SC-m", select = c(pPKCG_N,pNUMB_N,pNR1_N,pCAMKII_N,class))
summary(myIntCortex)
myIntCortex$pPKCG_Ncat <- cut(myIntCortex$pPKCG_N,quantile(myIntCortex$pPKCG_N, probs = c(0,0.5,1)),include.lowest =TRUE)
myIntCortex$pNUMB_Ncat <- cut(myIntCortex$pNUMB_N,quantile(myIntCortex$pNUMB_N, probs = c(0,0.5,1)),include.lowest =TRUE)
myIntCortex$pNR1_Ncat  <- cut(myIntCortex$pNR1_N,quantile(myIntCortex$pNR1_N, probs = c(0,0.5,1)),include.lowest =TRUE)
myIntCortex$pCAMKII_Ncat  <- cut(myIntCortex$pCAMKII_N,quantile(myIntCortex$pCAMKII_N, probs = c(0,0.5,1)),include.lowest =TRUE)
myIntdata <- data.frame(myIntCortex$class,myIntCortex$pPKCG_Ncat,myIntCortex$pNUMB_Ncat,myIntCortex$pNR1_Ncat,myIntCortex$pCAMKII_Ncat)                    
names(myIntdata) <- c("V1","V2","V3","V4","V5") #This relabels the variables in the dataset with the variable names produced by CStrees(p,2).
levels(myIntdata$V2) <- c(1:2) #This relabels the outcomes to match the outcome names used by CStrees(p,2).
levels(myIntdata$V3) <- c(1:2)
levels(myIntdata$V4) <- c(1:2)
levels(myIntdata$V5) <- c(1:2) 
summary(myIntdata)
# We then combine the discretizeed observed and interventional data sets into a single data frame.
myMixdata <- rbind(myObsdata,myIntdata)
myMixdata$V1 <- factor(myMixdata$V1)
levels(myMixdata$V1) <- factor(c("Int","Obs"))
myMixdata$V1<-relevel(myMixdata$V1,"Obs")
summary(myMixdata)

# We now construct the stage partition of the observed tree which is needed as input for fittedInterventionalCStrees.
mecParts<- list(list(list(list("1","1")),
                     list(list("1","3"),list("2","4")),
                     list(list("1","3"),list("2","4"),list("5","6","7","8"))
))

# Producing the different fitted interventional CStrees possible from intervening on the tree with staging mecParts.
listInterventionalCStrees<-list()
for (part in mecParts){
  ObsTree<-toCStree(4,2,part)
  Smore <- stages(ObsTree)
  S<-lapply(Smore, unique)
  S<-c(list(list(1)),S)
  newInterventionalTrees <-fittedInterventionalCStrees(4,2,part,S,myMixdata)
  listInterventionalCStrees<-c(listInterventionalCStrees,newInterventionalTrees)
}

# Find the fitted interventional CStree with the highest score.
MM <- listInterventionalCStrees[[1]]
t <- listInterventionalCStrees[[1]][[2]][[1]][1]
for (N in listInterventionalCStrees) {
  s <- N[[2]][[1]][1]
  if (s > t) {
    MM <- N
    t <- s
  }
}

# The optimal scoring interventional CStree is
optInt1 <- MM
plot(optInt1[[1]][[1]])
# Its score is
optIntScore1<-t
optIntScore1
# We begin a list containing the optimal trees learned from each member of the equvialence class
optIntTrees <- list(list(optInt1))



#-------------- Tree #2 --------------------------------------------
## Causal ordering = 3 2 1 4
myObsdata <- data.frame(myObsCortex$class,myObsCortex$pPKCG_Ncat,myObsCortex$pNUMB_Ncat,myObsCortex$pNR1_Ncat,myObsCortex$pCAMKII_Ncat)
names(myObsdata) <- c("V1","V4","V3","V2","V5")
levels(myObsdata$V2) <- c(1:2) #This relabels the outcomes to match the outcome names used by CStrees(p,2).
levels(myObsdata$V3) <- c(1:2)
levels(myObsdata$V4) <- c(1:2)
levels(myObsdata$V5) <- c(1:2) 
summary(myObsdata)
#Then we need to build the interventional discretized data set.
myIntCortex<- subset(myCortexdf,myCortexdf$class == "c-SC-m", select = c(pPKCG_N,pNUMB_N,pNR1_N,pCAMKII_N,class))
summary(myIntCortex)
myIntCortex$pPKCG_Ncat <- cut(myIntCortex$pPKCG_N,quantile(myIntCortex$pPKCG_N, probs = c(0,0.5,1)),include.lowest =TRUE)
myIntCortex$pNUMB_Ncat <- cut(myIntCortex$pNUMB_N,quantile(myIntCortex$pNUMB_N, probs = c(0,0.5,1)),include.lowest =TRUE)
myIntCortex$pNR1_Ncat  <- cut(myIntCortex$pNR1_N,quantile(myIntCortex$pNR1_N, probs = c(0,0.5,1)),include.lowest =TRUE)
myIntCortex$pCAMKII_Ncat  <- cut(myIntCortex$pCAMKII_N,quantile(myIntCortex$pCAMKII_N, probs = c(0,0.5,1)),include.lowest =TRUE)
myIntdata <- data.frame(myIntCortex$class,myIntCortex$pPKCG_Ncat,myIntCortex$pNUMB_Ncat,myIntCortex$pNR1_Ncat,myIntCortex$pCAMKII_Ncat)                    
names(myIntdata) <- c("V1","V4","V3","V2","V5") #This relabels the variables in the dataset with the variable names produced by CStrees(p,2).
levels(myIntdata$V2) <- c(1:2) #This relabels the outcomes to match the outcome names used by CStrees(p,2).
levels(myIntdata$V3) <- c(1:2)
levels(myIntdata$V4) <- c(1:2)
levels(myIntdata$V5) <- c(1:2) 
summary(myIntdata)
# We then combine the discretizeed observed and interventional data sets into a single data frame.
myMixdata <- rbind(myObsdata,myIntdata)
myMixdata$V1 <- factor(myMixdata$V1)
levels(myMixdata$V1) <- factor(c("Int","Obs"))
myMixdata$V1<-relevel(myMixdata$V1,"Obs")
summary(myMixdata)

# We now construct the stage partition of the observed tree which is needed as input for fittedInterventionalCStrees.
mecParts<- list(list(list(list("1","1")),
                     list(list("1","3"),list("2","4")),
                     list(list("1","3"),list("5","7"),list("2","4","6","8"))
))

# Producing the different fitted interventional CStrees possible from intervening on the tree with staging mecParts.
listInterventionalCStrees<-list()
for (part in mecParts){
  ObsTree<-toCStree(4,2,part)
  Smore <- stages(ObsTree)
  S<-lapply(Smore, unique)
  S<-c(list(list(1)),S)
  newInterventionalTrees <-fittedInterventionalCStrees(4,2,part,S,myMixdata)
  listInterventionalCStrees<-c(listInterventionalCStrees,newInterventionalTrees)
}

# Find the fitted interventional CStree with the highest score.
MM <- listInterventionalCStrees[[1]]
t <- listInterventionalCStrees[[1]][[2]][[1]][1]
for (N in listInterventionalCStrees) {
  s <- N[[2]][[1]][1]
  if (s > t) {
    MM <- N
    t <- s
  }
}

# The optimal scoring interventional CStree is
optInt2 <- MM
plot(optInt2[[1]][[1]])
# Its score is
optIntScore2<-t
optIntScore2
# We add the learned tree to the list of optimal trees learned.
optIntTrees <- c(optIntTrees,list(optInt2))


#-------------- Tree #3 --------------------------------------------
## Causal ordering = 2 1 3 4
myObsdata <- data.frame(myObsCortex$class,myObsCortex$pPKCG_Ncat,myObsCortex$pNUMB_Ncat,myObsCortex$pNR1_Ncat,myObsCortex$pCAMKII_Ncat)
names(myObsdata) <- c("V1","V3","V2","V4","V5")
levels(myObsdata$V2) <- c(1:2) #This relabels the outcomes to match the outcome names used by CStrees(p,2).
levels(myObsdata$V3) <- c(1:2)
levels(myObsdata$V4) <- c(1:2)
levels(myObsdata$V5) <- c(1:2) 
summary(myObsdata)
#Then we need to build the interventional discretized data set.
myIntCortex<- subset(myCortexdf,myCortexdf$class == "c-SC-m", select = c(pPKCG_N,pNUMB_N,pNR1_N,pCAMKII_N,class))
summary(myIntCortex)
myIntCortex$pPKCG_Ncat <- cut(myIntCortex$pPKCG_N,quantile(myIntCortex$pPKCG_N, probs = c(0,0.5,1)),include.lowest =TRUE)
myIntCortex$pNUMB_Ncat <- cut(myIntCortex$pNUMB_N,quantile(myIntCortex$pNUMB_N, probs = c(0,0.5,1)),include.lowest =TRUE)
myIntCortex$pNR1_Ncat  <- cut(myIntCortex$pNR1_N,quantile(myIntCortex$pNR1_N, probs = c(0,0.5,1)),include.lowest =TRUE)
myIntCortex$pCAMKII_Ncat  <- cut(myIntCortex$pCAMKII_N,quantile(myIntCortex$pCAMKII_N, probs = c(0,0.5,1)),include.lowest =TRUE)
myIntdata <- data.frame(myIntCortex$class,myIntCortex$pPKCG_Ncat,myIntCortex$pNUMB_Ncat,myIntCortex$pNR1_Ncat,myIntCortex$pCAMKII_Ncat)                    
names(myIntdata) <- c("V1","V3","V2","V4","V5") #This relabels the variables in the dataset with the variable names produced by CStrees(p,2).
levels(myIntdata$V2) <- c(1:2) #This relabels the outcomes to match the outcome names used by CStrees(p,2).
levels(myIntdata$V3) <- c(1:2)
levels(myIntdata$V4) <- c(1:2)
levels(myIntdata$V5) <- c(1:2) 
summary(myIntdata)
# We then combine the discretizeed observed and interventional data sets into a single data frame.
myMixdata <- rbind(myObsdata,myIntdata)
myMixdata$V1 <- factor(myMixdata$V1)
levels(myMixdata$V1) <- factor(c("Int","Obs"))
myMixdata$V1<-relevel(myMixdata$V1,"Obs")
summary(myMixdata)

# We now construct the stage partition of the observed tree which is needed as input for fittedInterventionalCStrees.
mecParts<- list(list(list(list("1","1")),
                     list(list("1","2"),list("3","4")),
                     list(list("1","5"),list("2","6"),list("3","4","7","8"))
))

# Producing the different fitted interventional CStrees possible from intervening on the tree with staging mecParts.
listInterventionalCStrees<-list()
for (part in mecParts){
  ObsTree<-toCStree(4,2,part)
  Smore <- stages(ObsTree)
  S<-lapply(Smore, unique)
  S<-c(list(list(1)),S)
  newInterventionalTrees <-fittedInterventionalCStrees(4,2,part,S,myMixdata)
  listInterventionalCStrees<-c(listInterventionalCStrees,newInterventionalTrees)
}

# Find the fitted interventional CStree with the highest score.
MM <- listInterventionalCStrees[[1]]
t <- listInterventionalCStrees[[1]][[2]][[1]][1]
for (N in listInterventionalCStrees) {
  s <- N[[2]][[1]][1]
  if (s > t) {
    MM <- N
    t <- s
  }
}

# The optimal scoring interventional CStree is
optInt3 <- MM
plot(optInt3[[1]][[1]])
# Its score is
optIntScore3<-t
optIntScore3
# We add the learned tree to the list of optimal trees learned.
optIntTrees <- c(optIntTrees,list(optInt3))


#-------------- Tree #4 --------------------------------------------
## Causal ordering = 2 3 1 4
myObsdata <- data.frame(myObsCortex$class,myObsCortex$pPKCG_Ncat,myObsCortex$pNUMB_Ncat,myObsCortex$pNR1_Ncat,myObsCortex$pCAMKII_Ncat)
names(myObsdata) <- c("V1","V4","V2","V3","V5")
levels(myObsdata$V2) <- c(1:2) #This relabels the outcomes to match the outcome names used by CStrees(p,2).
levels(myObsdata$V3) <- c(1:2)
levels(myObsdata$V4) <- c(1:2)
levels(myObsdata$V5) <- c(1:2) 
summary(myObsdata)
#Then we need to build the interventional discretized data set.
myIntCortex<- subset(myCortexdf,myCortexdf$class == "c-SC-m", select = c(pPKCG_N,pNUMB_N,pNR1_N,pCAMKII_N,class))
summary(myIntCortex)
myIntCortex$pPKCG_Ncat <- cut(myIntCortex$pPKCG_N,quantile(myIntCortex$pPKCG_N, probs = c(0,0.5,1)),include.lowest =TRUE)
myIntCortex$pNUMB_Ncat <- cut(myIntCortex$pNUMB_N,quantile(myIntCortex$pNUMB_N, probs = c(0,0.5,1)),include.lowest =TRUE)
myIntCortex$pNR1_Ncat  <- cut(myIntCortex$pNR1_N,quantile(myIntCortex$pNR1_N, probs = c(0,0.5,1)),include.lowest =TRUE)
myIntCortex$pCAMKII_Ncat  <- cut(myIntCortex$pCAMKII_N,quantile(myIntCortex$pCAMKII_N, probs = c(0,0.5,1)),include.lowest =TRUE)
myIntdata <- data.frame(myIntCortex$class,myIntCortex$pPKCG_Ncat,myIntCortex$pNUMB_Ncat,myIntCortex$pNR1_Ncat,myIntCortex$pCAMKII_Ncat)                    
names(myIntdata) <- c("V1","V4","V2","V3","V5") #This relabels the variables in the dataset with the variable names produced by CStrees(p,2).
levels(myIntdata$V2) <- c(1:2) #This relabels the outcomes to match the outcome names used by CStrees(p,2).
levels(myIntdata$V3) <- c(1:2)
levels(myIntdata$V4) <- c(1:2)
levels(myIntdata$V5) <- c(1:2) 
summary(myIntdata)
# We then combine the discretizeed observed and interventional data sets into a single data frame.
myMixdata <- rbind(myObsdata,myIntdata)
myMixdata$V1 <- factor(myMixdata$V1)
levels(myMixdata$V1) <- factor(c("Int","Obs"))
myMixdata$V1<-relevel(myMixdata$V1,"Obs")
summary(myMixdata)

# We now construct the stage partition of the observed tree which is needed as input for fittedInterventionalCStrees.
mecParts<- list(list(list(list("1","1")),
                     list(list("1","2"),list("3","4")),
                     list(list("1","5"),list("3","7"),list("2","4","6","8"))
))

# Producing the different fitted interventional CStrees possible from intervening on the tree with staging mecParts.
listInterventionalCStrees<-list()
for (part in mecParts){
  ObsTree<-toCStree(4,2,part)
  Smore <- stages(ObsTree)
  S<-lapply(Smore, unique)
  S<-c(list(list(1)),S)
  newInterventionalTrees <-fittedInterventionalCStrees(4,2,part,S,myMixdata)
  listInterventionalCStrees<-c(listInterventionalCStrees,newInterventionalTrees)
}

# Find the fitted interventional CStree with the highest score.
MM <- listInterventionalCStrees[[1]]
t <- listInterventionalCStrees[[1]][[2]][[1]][1]
for (N in listInterventionalCStrees) {
  s <- N[[2]][[1]][1]
  if (s > t) {
    MM <- N
    t <- s
  }
}

# The optimal scoring interventional CStree is
optInt4 <- MM
plot(optInt4[[1]][[1]])
# Its score is
optIntScore4<-t
optIntScore4
# We add the learned tree to the list of optimal trees learned.
optIntTrees <- c(optIntTrees,list(optInt4))

# We finally select the interventional CStree from the list optIntTrees with the highest score.
for (i in c(1: length(optIntTrees))) {
  if (eval(parse(text = paste("optIntScore",i,sep = ""))) == max(optIntScore1,optIntScore2,optIntScore3,optIntScore4)) {
    print(i)
  }
}

# Hence the optimal interventional CStree is optInt2 which has causal ordering 3 2 1 4.
# Thus, in the optimal tree V2 = pNR1, V3 = pNUMB, V4 = pPKCG, and V5 = pCAMKII.
plot(optInt2[[1]][[1]])


















