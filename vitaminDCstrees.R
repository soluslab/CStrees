# VitaminD
# Causal discovery with CStrees and interventional CStrees


#----------------- Interventional CStrees ---------------%
#Learning an interventionalCStree on four levels using the VitD data
#-------------- Tree #1
myObsvData <- myTotalData[which(myTotalData$filaggrin == 0),]
x<-c(200:400) 
myObsvData<- myObsvData[x,]
#---
myData <- VitD
myData <- subset(myData,select = -c(filaggrin,death)) #truncate the data set to remove all discrete rows.
myData <- discretize(myData,method="quantile",breaks=2) #discretize all non-discrete variables.
myData <- cbind(myData,VitD$death,VitD$filaggrin) #append non-interventional discrete variable back into data frame.
#----
myTotalData <- VitD
myObsvData <- myTotalData[which(myTotalData$filaggrin == 0),]
x<-c(200:400) 
myObsvData<- myObsvData[x,]
myInvData<- myTotalData[which(myTotalData$filaggrin==1),]
myShortData<-rbind(myObsvData,myInvData)
summary(myShortData)
nrow(myShortData)
myData <- subset(myShortData,select = -c(filaggrin,death)) #truncate the data set to remove all discrete rows.
myData <- discretize(myData,method="quantile",breaks=2) #discretize all non-discrete variables.
myData <- cbind(myData,myShortData$death,myShortData$filaggrin)
#-----
names(myData)
names(myData) <- c("V3","V4","V5","V2","V1") #Relabel the variables in the dataset with the variable names produced by CStrees(p,d)
levels(myData$V1) <- factor(c("Obs","Int"))
levels(myData$V2) <- c(1:2)
levels(myData$V3) <- c(1:2)
levels(myData$V4) <- c(1:2)
levels(myData$V5) <- c(1:2)
# From here on we construct the interventional CStrees
mecParts<- list(list(list(list("1","1")),
                     list(list("1","2","3","4")),
                     list(list("1","3"),list("2","4"),list("5","6","7","8"))
))
##----## Linear extension = 4 1 2 3
# V1=age, V2=vitaminD, V3= time, V4=Death
L1 = list(list("1","1"))
L2 = list(list("1","2","3","4"))
L3 = list(list("1","3"),list("2","4"),list("5","6","7","8"))
#--------------------------------------------------------------
listInterventionalCStrees<-list()
for (part in mecParts){
  ObsTree<-toCStree(4,2,part)
  Smore <- stages(ObsTree)
  S<-lapply(Smore, unique)
  S<-c(list(list(1)),S)
  newInterventionalTrees <-interventionalCStrees(4,2,part,S)
  listInterventionalCStrees<-c(listInterventionalCStrees,newInterventionalTrees)
}

M<-listInterventionalCStrees[[14]] # initialize with some model
plot(M)
M.fit <- sevt_fit(M,myData,lambda=1)
obsM <- subtree(M.fit,c("Obs"))
intM <- subtree(M.fit,c("Int"))
stagesByLevel<- lapply(stages(M.fit),unique)
stagesByLevel<- lapply(stagesByLevel, length)
numStages<- Reduce("+",stagesByLevel) +1
numStages
# Now we evaluate the loglik on each of the subtrees
t<- (logLik(obsM) + logLik(intM))*2
#t<- logLik(M.fit)*2-numStages*log(nrow(myData))
MM <- M.fit
for (N in listInterventionalCStrees) {
  N.fit <- sevt_fit(N,myData,lambda=1)
  obsM <- subtree(N.fit,c("Obs"))
  intM <- subtree(N.fit,c("Int"))
  stagesByLevel<- lapply(stages(N.fit),unique)
  stagesByLevel<- lapply(stagesByLevel, length)
  numStages<- Reduce("+",stagesByLevel) +1
#  s<- logLik(N.fit)*2-numStages*log(nrow(myData))
  s <- (logLik(obsM) + logLik(intM))*2
  if (s > t) {
    MM <- N.fit
    t <- s
  }
}
intTT1 <- MM
plot(intTT1)
t1<-t
t
logLik(obsM)
logLik(intM)
#-------------------------------------

#-------------- Tree #2
#Learning an interventionalCStree on three levels using the VitD data
myData <- VitD
myData <- subset(myData,select = -c(filaggrin,death)) #truncate the data set to remove all discrete rows.
myData <- discretize(myData,method="quantile",breaks=2) #discretize all non-discrete variables.
myData <- cbind(myData,VitD$death,VitD$filaggrin) #append non-interventional discrete variable back into data frame.
names(myData)
names(myData) <- c("V4","V3","V5","V2","V1") #Relabel the variables in the dataset with the variable names produced by CStrees(p,d)
levels(myData$V1) <- factor(c("Int","Obs"))
levels(myData$V2) <- c(1:2)
levels(myData$V3) <- c(1:2)
levels(myData$V4) <- c(1:2)
levels(myData$V5) <- c(1:2)
# From here on we construct the interventional CStrees
mecParts<- list(list(list(list("1","2")),
                     list(list("1","2"),list("3","4")),
                     list(list("1","2"),list("3","4"),list("5","6","7","8"))
))
## Linear extension = 4 2 1 3
# V1=age, V2=vitaminD, V3= time, V4=Death
L1 = list(list("1","2"))
L2 = list(list("1","2"),list("3","4"))
L3 = list(list("1","2"),list("3","4"),list("5","6","7","8"))
#--------------------------------------------------------------
listInterventionalCStrees<-list()
for (part in mecParts){
  ObsTree<-toCStree(4,2,part)
  Smore <- stages(ObsTree)
  S<-lapply(Smore, unique)
  S<-c(list(list(1)),S)
  newInterventionalTrees <-interventionalCStrees(4,2,part,S)
  listInterventionalCStrees<-c(listInterventionalCStrees,newInterventionalTrees)
}

M<-listInterventionalCStrees[[1]] # initialize with some model
plot(M)
M.fit <- sevt_fit(M,myData,lambda=1)
obsM <- subtree(M.fit,c("Obs"))
intM <- subtree(M.fit,c("Int"))
stagesByLevel<- lapply(stages(M.fit),unique)
stagesByLevel<- lapply(stagesByLevel, length)
numStages<- Reduce("+",stagesByLevel) +1
numStages
# Now we evaluate the loglik on each of the subtrees
t<- (logLik(obsM) + logLik(intM))*2-numStages*log(nrow(myData))
MM <- M.fit
for (N in listInterventionalCStrees) {
  N.fit <- sevt_fit(N,myData,lambda=1)
  obsM <- subtree(N.fit,c("Obs"))
  intM <- subtree(N.fit,c("Int"))
  stagesByLevel<- lapply(stages(N.fit),unique)
  stagesByLevel<- lapply(stagesByLevel, length)
  numStages<- Reduce("+",stagesByLevel) +1
  s <- (logLik(obsM) + logLik(intM))*2-numStages*log(nrow(myData))
  if (s > t) {
    MM <- N.fit
    t <- s
  }
}
intTT2 <- MM
plot(intTT2)
t2<-t
t2
#---- Tree #3
myData <- VitD
myData <- subset(myData,select = -c(filaggrin,death)) #truncate the data set to remove all discrete rows.
myData <- discretize(myData,method="quantile",breaks=2) #discretize all non-discrete variables.
myData <- cbind(myData,VitD$death,VitD$filaggrin) #append non-interventional discrete variable back into data frame.
names(myData)
names(myData) <- c("V4","V2","V5","V3","V1") #Relabel the variables in the dataset with the variable names produced by CStrees(p,d)
levels(myData$V1) <- factor(c("Int","Obs"))
levels(myData$V2) <- c(1:2)
levels(myData$V3) <- c(1:2)
levels(myData$V4) <- c(1:2)
levels(myData$V5) <- c(1:2)
# From here on we construct the interventional CStrees
mecParts<- list(list(list(list("1","2")),
                     list(list("1","3"),list("2","4")),
                     list(list("1","2"),list("5","6"),list("3","4","7","8"))
))
# V1=age, V2=vitaminD, V3= time, V4=Death
## Linear extension = 2 4 1 3
L1 = list(list("1","2"))
L2 = list(list("1","3"),list("2","4"))
L3 = list(list("1","2"),list("5","6"),list("3","4","7","8"))
TT3 <- toCStree(4,2,list(L1,L2,L3))
plot(TT3)
#--------------------------------------------------------------
listInterventionalCStrees<-list()
for (part in mecParts){
  ObsTree<-toCStree(4,2,part)
  Smore <- stages(ObsTree)
  S<-lapply(Smore, unique)
  S<-c(list(list(1)),S)
  newInterventionalTrees <-interventionalCStrees(4,2,part,S)
  listInterventionalCStrees<-c(listInterventionalCStrees,newInterventionalTrees)
}
length(listInterventionalCStrees)
M<-listInterventionalCStrees[[120]] # initialize with some model
plot(M)
M.fit <- sevt_fit(M,myData,lambda=1)
obsM <- subtree(M.fit,c("Obs"))
intM <- subtree(M.fit,c("Int"))
stagesByLevel<- lapply(stages(M.fit),unique)
stagesByLevel<- lapply(stagesByLevel, length)
numStages<- Reduce("+",stagesByLevel) +1
numStages
# Now we evaluate the loglik on each of the subtrees
t<- (logLik(obsM) + logLik(intM))*2-numStages*log(nrow(myData))
MM <- M.fit
for (N in listInterventionalCStrees) {
  N.fit <- sevt_fit(N,myData,lambda=1)
  obsM <- subtree(N.fit,c("Obs"))
  intM <- subtree(N.fit,c("Int"))
  stagesByLevel<- lapply(stages(N.fit),unique)
  stagesByLevel<- lapply(stagesByLevel, length)
  numStages<- Reduce("+",stagesByLevel) +1
  s <- (logLik(obsM) + logLik(intM))*2-numStages*log(nrow(myData))
  if (s > t) {
    MM <- N.fit
    t <- s
  }
}
intTT3 <- MM
plot(intTT3)
t3<-t

#---
#---- Tree #4
myData <- VitD
myData <- subset(myData,select = -c(filaggrin,death)) #truncate the data set to remove all discrete rows.
myData <- discretize(myData,method="quantile",breaks=2) #discretize all non-discrete variables.
myData <- cbind(myData,VitD$death,VitD$filaggrin) #append non-interventional discrete variable back into data frame.
names(myData)
names(myData) <- c("V2","V4","V5","V3","V1") #Relabel the variables in the dataset with the variable names produced by CStrees(p,d)
levels(myData$V1) <- factor(c("Int","Obs"))
levels(myData$V2) <- c(1:2)
levels(myData$V3) <- c(1:2)
levels(myData$V4) <- c(1:2)
levels(myData$V5) <- c(1:2)
# From here on we construct the interventional CStrees
mecParts<- list(list(list(list("1","2")),
                     list(list("1","3"),list("2","4")),
                     list(list("1","5"),list("2","6"),list("3","4","7","8"))
))
# V1=age, V2=vitaminD, V3= time, V4=Death
# Fourth tree in the equivalence class
## Linear extension = 1 4 2 3
L1 = list(list("1","1"))
L2 = list(list("1","3","2","4"))
L3 = list(list("1","5"),list("2","6"),list("3","4","7","8"))
TT4 <- toCStree(4,2,list(L1,L2,L3))
plot(TT4)
#--------------------------------------------------------------
listInterventionalCStrees<-list()
for (part in mecParts){
  ObsTree<-toCStree(4,2,part)
  Smore <- stages(ObsTree)
  S<-lapply(Smore, unique)
  S<-c(list(list(1)),S)
  newInterventionalTrees <-interventionalCStrees(4,2,part,S)
  listInterventionalCStrees<-c(listInterventionalCStrees,newInterventionalTrees)
}

M<-listInterventionalCStrees[[10]] # initialize with some model
plot(M)
M.fit <- sevt_fit(M,myData,lambda=1)
obsM <- subtree(M.fit,c("Obs"))
intM <- subtree(M.fit,c("Int"))
stagesByLevel<- lapply(stages(M.fit),unique)
stagesByLevel<- lapply(stagesByLevel, length)
numStages<- Reduce("+",stagesByLevel) +1
numStages
# Now we evaluate the loglik on each of the subtrees
t<- (logLik(obsM) + logLik(intM))*2-numStages*log(nrow(myData))
MM <- M.fit
for (N in listInterventionalCStrees) {
  N.fit <- sevt_fit(N,myData,lambda=1)
  obsM <- subtree(N.fit,c("Obs"))
  intM <- subtree(N.fit,c("Int"))
  stagesByLevel<- lapply(stages(N.fit),unique)
  stagesByLevel<- lapply(stagesByLevel, length)
  numStages<- Reduce("+",stagesByLevel) +1
  s <- (logLik(obsM) + logLik(intM))*2-numStages*log(nrow(myData))
  if (s > t) {
    MM <- N.fit
    t <- s
  }
}
intTT4 <- MM
plot(intTT4)
t4<-t
#----

#---- Tree #5
myData <- VitD
myData <- subset(myData,select = -c(filaggrin,death)) #truncate the data set to remove all discrete rows.
myData <- discretize(myData,method="quantile",breaks=2) #discretize all non-discrete variables.
myData <- cbind(myData,VitD$death,VitD$filaggrin) #append non-interventional discrete variable back into data frame.
names(myData)
names(myData) <- c("V2","V3","V5","V4","V1") #Relabel the variables in the dataset with the variable names produced by CStrees(p,d)
levels(myData$V1) <- factor(c("Int","Obs"))
levels(myData$V2) <- c(1:2)
levels(myData$V3) <- c(1:2)
levels(myData$V4) <- c(1:2)
levels(myData$V5) <- c(1:2)
# From here on we construct the interventional CStrees
mecParts<- list(list(list(list("1","2")),
                     list(list("1","2"),list("3","4")),
                     list(list("1","5"),list("3","7"),list("2","4","6","8"))
))
# V1=age, V2=vitaminD, V3= time, V4=Death
# Fifth tree in the partition
## Linear extension = 1 2 4 3
L1 = list(list("1","2"))
L2 = list(list("1","2"),list("3","4"))
L3 = list(list("1","5"),list("3","7"),list("2","4","6","8"))
TT5 <- toCStree(4,2,list(L1,L2,L3))
plot(TT5)
t5<-t
#--------------------------------------------------------------
listInterventionalCStrees<-list()
for (part in mecParts){
  ObsTree<-toCStree(4,2,part)
  Smore <- stages(ObsTree)
  S<-lapply(Smore, unique)
  S<-c(list(list(1)),S)
  newInterventionalTrees <-interventionalCStrees(4,2,part,S)
  listInterventionalCStrees<-c(listInterventionalCStrees,newInterventionalTrees)
}

M<-listInterventionalCStrees[[10]] # initialize with some model
plot(M)
M.fit <- sevt_fit(M,myData,lambda=1)
obsM <- subtree(M.fit,c("Obs"))
intM <- subtree(M.fit,c("Int"))
stagesByLevel<- lapply(stages(M.fit),unique)
stagesByLevel<- lapply(stagesByLevel, length)
numStages<- Reduce("+",stagesByLevel) +1
numStages
# Now we evaluate the loglik on each of the subtrees
t<- (logLik(obsM) + logLik(intM))*2-numStages*log(nrow(myData))
MM <- M.fit
for (N in listInterventionalCStrees) {
  N.fit <- sevt_fit(N,myData,lambda=1)
  obsM <- subtree(N.fit,c("Obs"))
  intM <- subtree(N.fit,c("Int"))
  stagesByLevel<- lapply(stages(N.fit),unique)
  stagesByLevel<- lapply(stagesByLevel, length)
  numStages<- Reduce("+",stagesByLevel) +1
  s <- (logLik(obsM) + logLik(intM))*2-numStages*log(nrow(myData))
  if (s > t) {
    MM <- N.fit
    t <- s
  }
}
intTT5 <- MM
plot(intTT5)
t6<-t
#----


#---- Tree #6
myData <- VitD
myData <- subset(myData,select = -c(filaggrin,death)) #truncate the data set to remove all discrete rows.
myData <- discretize(myData,method="quantile",breaks=2) #discretize all non-discrete variables.
myData <- cbind(myData,VitD$death,VitD$filaggrin) #append non-interventional discrete variable back into data frame.
names(myData)
names(myData) <- c("V3","V2","V5","V4","V1") #Relabel the variables in the dataset with the variable names produced by CStrees(p,d)
levels(myData$V1) <- factor(c("Int","Obs"))
levels(myData$V2) <- c(1:2)
levels(myData$V3) <- c(1:2)
levels(myData$V4) <- c(1:2)
levels(myData$V5) <- c(1:2)
# From here on we construct the interventional CStrees
mecParts<- list(list(list(list("1","2")),
                     list(list("1","3"),list("2","4")),
                     list(list("1","3"),list("5","7"),list("2","4","6","8"))
))
# V1=age, V2=vitaminD, V3= time, V4=Death
#Sixth tree in the partition
## Linear extension = 2 1 4 3
L1 = list(list("1","2"))
L2 = list(list("1","3"),list("2","4"))
L3 = list(list("1","3"),list("5","7"),list("2","4","6","8"))
TT6 <- toCStree(4,2,list(L1,L2,L3))
plot(TT6)
#--------------------------------------------------------------
listInterventionalCStrees<-list()
for (part in mecParts){
  ObsTree<-toCStree(4,2,part)
  Smore <- stages(ObsTree)
  S<-lapply(Smore, unique)
  S<-c(list(list(1)),S)
  newInterventionalTrees <-interventionalCStrees(4,2,part,S)
  listInterventionalCStrees<-c(listInterventionalCStrees,newInterventionalTrees)
}

M<-listInterventionalCStrees[[20]] # initialize with some model
plot(M)
M.fit <- sevt_fit(M,myData,lambda=1)
obsM <- subtree(M.fit,c("Obs"))
intM <- subtree(M.fit,c("Int"))
stagesByLevel<- lapply(stages(M.fit),unique)
stagesByLevel<- lapply(stagesByLevel, length)
numStages<- Reduce("+",stagesByLevel) +1
numStages
# Now we evaluate the loglik on each of the subtrees
t<- (logLik(obsM) + logLik(intM))*2-numStages*log(nrow(myData))
MM <- M.fit
for (N in listInterventionalCStrees) {
  N.fit <- sevt_fit(N,myData,lambda=1)
  obsM <- subtree(N.fit,c("Obs"))
  intM <- subtree(N.fit,c("Int"))
  stagesByLevel<- lapply(stages(N.fit),unique)
  stagesByLevel<- lapply(stagesByLevel, length)
  numStages<- Reduce("+",stagesByLevel) +1
  s <- (logLik(obsM) + logLik(intM))*2-numStages*log(nrow(myData))
  if (s > t) {
    MM <- N.fit
    t <- s
  }
}
intTT6 <- MM
plot(intTT6)
t
#----

plot(listInterventionalCStrees[[10]])

# Summary for results for interventional learning
# V1=age, V2=vitaminD, V3= time, V4=Death, 
#--------------------------------------------------
# Summary of Obsv model: 
# N = 2377
# BIC = 11240.71
# logLik = -5593.146 (df=7)
# Tree#1 ---------------------------------------------
# 4,1,2,3,BIC = -12217.49 (df=7), -12171.06 (df=7)
plot(intTT1)
# Tree#2 ---------------------------------------------
# 4 2 1 3, BIC = -12171.06 (df=7), -12171.06 (df=7)
plot(intTT2)
# Tree#3 ---------------------------------------------
# 2 4 1 3, BIC = -12171.06 (df=7), -12171.06 (df=7)
plot(intTT3)
# Tree#4  --------------------------------------------
# 1 4 2 3, BIC = -12582.61 (df=7)
plot(intTT4)
# Tree#5 ---------------------------------------------
# 1 2 4 3, BIC = -12171.06 (df=7)
plot(intTT5)
# Tree#6 ---------------------------------------------
# 2 1 4 3, BIC = -12171.06 (df=7)
plot(intTT56)
#---------------------------------------------------

# Partitions for statistically equivalent trees with their linear
# extensions

## Trees in the equivalence class of the learned tree.
## First tree in the equivalence class
## Linear extension = 4 1 2 3
L1 = list(list("1","1"))
L2 = list(list("1","2","3","4"))
L3 = list(list("1","3"),list("2","4"),list("5","6","7","8"))
TT1 <- toCStree(4,2,list(L1,L2,L3))
plot(newTT)
## Second tree in the equivalence class
## Linear extension = 4 2 1 3
L1 = list(list("1","2"))
L2 = list(list("1","2"),list("3","4"))
L3 = list(list("1","2"),list("3","4"),list("5","6","7","8"))
TT2 <- toCStree(4,2,list(L1,L2,L3))
## Third tree in the equivalence class
## Linear extension = 2 4 1 3
L1 = list(list("1","2"))
L2 = list(list("1","3"),list("2","4"))
L3 = list(list("1","2"),list("5","6"),list("3","4","7","8"))
TT3 <- toCStree(4,2,list(L1,L2,L3))
# Fourth tree in the equivalence class
## Linear extension = 1 4 2 3
L1 = list(list("1","1"))
L2 = list(list("1","3","2","4"))
L3 = list(list("1","5"),list("2","6"),list("3","4","7","8"))
TT4 <- toCStree(4,2,list(L1,L2,L3))
# Fifth tree in the partition
## Linear extension = 1 2 4 3
L1 = list(list("1","2"))
L2 = list(list("1","2"),list("3","4"))
L3 = list(list("1","5"),list("3","7"),list("2","4","6","8"))
TT5 <- toCStree(4,2,list(L1,L2,L3))
plot(TT5)
#Sixth tree in the partition
## Linear extension = 2 1 4 3
L1 = list(list("1","2"))
L2 = list(list("1","3"),list("2","4"))
L3 = list(list("1","3"),list("5","7"),list("2","4","6","8"))
TT6 <- toCStree(4,2,list(L1,L2,L3))



