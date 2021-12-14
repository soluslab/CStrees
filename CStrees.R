library("stagedtrees")
library("combinat")


# Producing the possible stagings for each level in a CStree. 
# We first generate the partitions of the 1 and 2 dimensional cubes.
OneDpartitions <- list(list(list("1","1")),
                       list(list("1","2")))
TwoDpartitions <- list(list(list("1","1")),
                       list(list("1","2")),
                       list(list("3","4")),
                       list(list("1","3")),
                       list(list("2","4")),
                       list(list("1","3"),list("2","4")),
                       list(list("1","2"),list("3","4")),
                       list(list("1","2","3","4")))

# We now produce the full set of partitions of the 3-dimensional cubes.  
# The following functions are used to read in the partitions of higher dimensional cubes into the correct format.
# Takes a string of characters and breaks it down into a list using the delimiter.
h <- function(mylongstring){
  strsplit(mylongstring,",")
}

# Takes a list and applies to it the function
# g which takes a string and breaks it down to a list by using a space as a delimiter
# then turns the character into an integer
f <- function(mylist){
  lapply(mylist,g)
}
g <- function(mystring){
  pieces <-(strsplit(mystring," "))
  lapply(pieces,as.list)
  #lapply(pieces, strtoi)
}

# flattens one level of the lists of lists
k <- function(mylist){
  return(unlist(mylist,recursive=FALSE))
}

# This part creates the ThreeDpartitions
z<-readLines(con = "threeDpartitions.txt")
listOfchars <- lapply(z, h)
ThreeDpartitions <-lapply(listOfchars,f)
ThreeDpartitions <-lapply(ThreeDpartitions,k)
ThreeDpartitions <- c(ThreeDpartitions,list(list(list("1","1"))))
length(ThreeDpartitions)

# We also produce a subset of the partitions of the 4-dimensional cube.  These partitions are those that can be
# recursively constructed by first partitioning the whole cube into two disjoint facets and then partitioning 
# these two facets.  
z<-readLines(con = "fourDpartitions.txt")
listOfchars <- lapply(z, h)
FourDpartitions <-lapply(listOfchars,f)
FourDpartitions <-lapply(FourDpartitions,k)
FourDpartitions <- c(FourDpartitions,list(list(list("1","1"))))
length(FourDpartitions)

# A list of lists where the i-th entry is the partitions of the i-dimensional cube.
ComputedPartitions = list(OneDpartitions,TwoDpartitions,ThreeDpartitions,FourDpartitions)

# A function for computing the cross product of partitions of the cubes of dimensions 1,...,p-1.  
# The function calls upon the set ComputedPartitions and can only handle values for p up to its cardinality plus 1.
partitions <- function(p) {
  Parts<-OneDpartitions
  if (p!=2) {
    for (i in c(2:(p-1))) {
      P <- list()
      for (L in Parts) {
        for (M in ComputedPartitions[[i]]) {
          if (i!=2) {
            P<-c(P,list(c(L,list(M))))
          } else {
            #if (i!=3) {
            # P<-c(P,c(L,list(M)))
            #} else {
            P<-c(P,list(list(L,M)))
          }
        }
        Parts <- P
      }
    }
  }
  return(Parts)
}

# A function that takes in a number of variables p , a size of outcome space of the variables d, 
# and a list of partitions L, one for each level of the tree, ordered by level from smallest-to-largest level.  
# It returns a CStree with the staging corresponding to the partitions given in L and causal ordering the
# identity permutation of [p].
toCStree <-function(p,d,L) {
  V <-as.list(c(1:p))
  Vars <- list()
  for (i in V) {
    Vars[[i]] <- c(1:d)
  }
  model.CS <- sevt(Vars,full=TRUE)
  for ( j in c(1:(p-1))) {
    t <- length(L[[j]])
    for (s in c(1:t)) {
      for (r in c(2:length(L[[j]][[s]]))) {
        if (L[[j]][[s]][[1]]!= L[[j]][[s]][[2]]) {
          k<-j+1
          model.CS <- join_stages(model.CS,paste("V",(j+1),sep = ""),L[[j]][[s]][[1]],L[[j]][[s]][[r]])
        }
      }
    }
  }
  return(model.CS)
}

# Example on four variables
plot(toCStree(4,2,partitions(4)[[1870]]))
partitions(4)[[1870]]

# Example on five variables.  Here we do not compute all of partitions(5) for computational efficiency.
# Instead, we build the stage list by directly specifying an element of the set partitions(5).
plot(toCStree(5,2,list(OneDpartitions[[2]],TwoDpartitions[[6]],ThreeDpartitions[[130]],FourDpartitions[[50000]])))
list(OneDpartitions[[2]],TwoDpartitions[[6]],ThreeDpartitions[[130]],FourDpartitions[[50000]])


# A function that produces all CStrees on p variables with each variable having d outcomes.  
# It requires that you have computed the set of all partitions of the i-dimensional cube for i< p-1, and added
# them to the list ComputedPartitions.
# For now, d = 2, as the set ComputedPartitions is defined for only binary variables.  
CStrees <- function(p,d) {
  Parts <- partitions(p) 
  Perms <- permn(as.list(c(1:p)))
  models <- list()
  for (i in c(1:length(Perms))) {
    Vars <- list()
    N <- list()
    for (j in c(1:p)) {
      Vars[[j]] <- c(1:d)
      N[[j]]<- paste("V",Perms[[i]][[j]],sep = "")
      names(Vars) <- N
    }
    model.CS <- sevt(Vars,full=TRUE)
    for (L in Parts) {
      M <- model.CS
 #     for ( j in c(1:(p-1))) {
 #       if (class(L[[j]][[1]]) == "numeric") {
 #         if (L[[j]][[1]] != L[[j]][[2]]) {
 #           for (r in c(2:length(L[[j]]))) {
 #             k<-j+1
 #             M <- join_stages(M,paste("V",Perms[[i]][[k]],sep = ""),L[[j]][[1]],L[[j]][[r]])
 #           }
 #         }
 #       } else {
        for (j in c(1:(p-1))) {
          t <- length(L[[j]])
          for (s in c(1:t)) {
            for (r in c(2:length(L[[j]][[s]]))) {
              if (L[[j]][[s]][[1]]!= L[[j]][[s]][[2]]) {
                k<-j+1
                M <- join_stages(M,paste("V",Perms[[i]][[k]],sep = ""),L[[j]][[s]][[1]],L[[j]][[s]][[r]])
              }
            }
          }
        }
    
      models <-c(models,list(M))
    }
  }
  return(models)
}



# Example:
plot(CStrees(3,2)[[1]])
summary(CStrees(3,2)[[1]])


# An application to real data:
library("bnlearn")

data(lizards) #dataset available within the bnlearn library
myData <- lizards
names(myData) <- c("V1","V2","V3") #This relabels the variables in the dataset with the variable names produced by CStrees(p,2).
levels(myData$V1) <- c(1:2) #This relabels the outcomes to match the outcome names used by CStrees(p,2).
levels(myData$V2) <- c(1:2)
levels(myData$V3) <- c(1:2)
CStreesList <- CStrees(3,2)

# The following will recover the BIC-optimal CStree for the dataset myData:
M <- CStreesList[[1]]
M.fit <- sevt_fit(M,data=myData,lambda=1)
t <- BIC(M.fit)
MM <- M.fit
summary(M.fit)
optTrees = list(1)
for (i in c(1:length(CStreesList))) {
  N.fit <- sevt_fit(CStreesList[[i]],data=myData,lambda=1)
  s <- BIC(N.fit)
  if (t > s) {
    MM <- N.fit
    t <- s
    optTrees <- list(i)
  } else {
    if (t==s) {
      optTrees <- c(optTrees,i)
    }
  }
}
# The first BIC-optimal tree in the list CStreesList is plotted, it and described below.  t is its BIC. 
plot(MM)
summary(MM)
t

# The following set contains the indices i for which CStreesList[[i]] is a BIC-optimal tree w.r.t. the myData dataset.
optTrees

# The original variable names and outcomes for each variable are given by the following:
names(lizards)
levels(lizards$Species)
levels(lizards$Diameter)
levels(lizards$Height)


# Fitting a CStree to a data set for coronary heart data set.
data(coronary) # dataset available within the bnlearn library
myData <- coronary
names(myData) <- c("V1","V2","V3","V4","V5","V6") # This function relabels the variables in the dataset with the variable names produced by CStrees(p,d)
levels(myData$V1) <- c(1:2)
levels(myData$V2) <- c(1:2)
levels(myData$V3) <- c(1:2)
levels(myData$V4) <- c(1:2)
levels(myData$V5) <- c(1:2)
levels(myData$V6) <- c(1:2)
myData <- subset(myData,select = -c(V4,V5,V6)) # Truncate to only variables S, M Work, and P Work.
# We look for context-specific structure amongst these three environmental factors.
CStreesList3 <- CStrees(3,2) 
M <- CStreesList3[[1]]
M.fit <- sevt_fit(M,data=myData,lambda=1)
t <- BIC(M.fit)
MM <- M.fit
for (N in CStreesList3) {
  N.fit <- sevt_fit(N,data=myData,lambda=1)
  s <- BIC(N.fit)
  if (t > s) {
    MM <- N.fit
    t <- s
  }
}
plot(MM)
summary(MM)
t


# Fitting a CStree to a data set for coronary heart data set for a set of different four variables.  
# Here we select a single one of the environment factors considered above and compare it with the other three biological indicators.
data(coronary) #dataset available within the bnlearn library
myData <- coronary
names(myData) <- c("V1","V5","V6","V2","V3","V4") # Relabels the variables in the dataset with the variable names produced by CStrees(p,d)
levels(myData$V1) <- c(1:2)
levels(myData$V2) <- c(1:2)
levels(myData$V3) <- c(1:2)
levels(myData$V4) <- c(1:2)
levels(myData$V5) <- c(1:2)
levels(myData$V6) <- c(1:2)
myData <- subset(myData,select = -c(V5,V6)) # Truncate the data set to include a single covariate from S, M Work, and P Work. In this example, we include S.
names(myData) <- c("V1","V2","V3","V4")
CStreesList <- CStrees(4,2) # Although the data set has 6 variables we fit a model to only the first four.
M <- CStreesList[[1]]
M.fit <- sevt_fit(M,data=myData,lambda=1)
t <- BIC(M.fit)
MM <- M.fit
for (N in CStreesList) {
  N.fit <- sevt_fit(N,data=myData,lambda=1)
  s <- BIC(N.fit)
  if (t > s) {
    MM <- N.fit
    t <- s
  }
}
plot(MM)
summary(MM)
t





# Interventional CStrees for one observed distribution and one interventional distribution.

# A function that takes in a number of variables p, a size of outcome space of the variables d, 
# a list of partitions L, one for each level of the tree, ordered by level from smallest-to-largest level, and a
# list of stages targeted for intervention.  I is broken into a list of lists, one list of stages for each level.
# Each stage in I must be a number between 1 and 2^(i-1) if the stage is associated to level i.
# It returns an interventional CStree with the staging corresponding to the partitions given in L and causal 
# ordering the identity permutation of [p].
toInterventionalCStree <-function(p,d,L,I) {
  ObsL <- L
  IntL <- L
  for (i in c(1:(p-1))) {
    m <- length(L[[i]])
    for (j in c(1:m)) {
      t <- length(L[[i]][[j]]) 
      for (k in c(1:t)) {
        IntL[[i]][[j]][[k]] <- toString(as.integer(IntL[[i]][[j]][[k]])+(2^(i)))
      }
    }
  }
  V <-as.list(c(1:(p+1)))
  Vars <- list(c("Obs","Int"))
  for (i in c(1:(p))) {
    Vars[[i+1]] <- c(1:d)
  }
  model.CS <- sevt(Vars,full=TRUE)
  for ( j in c(1:(p-1))) {
    t <- length(ObsL[[j]])
    for (s in c(1:t)) {
      for (r in c(2:length(ObsL[[j]][[s]]))) {
        if (ObsL[[j]][[s]][[1]]!= ObsL[[j]][[s]][[2]]) {
          model.CS <- join_stages(model.CS,paste("V",(j+2),sep = ""),ObsL[[j]][[s]][[1]],ObsL[[j]][[s]][[r]])
        }
      }
    }
  }
  for ( j in c(1:(p-1))) {
    t <- length(IntL[[j]])
    for (s in c(1:t)) {
      for (r in c(2:length(IntL[[j]][[s]]))) {
        if (IntL[[j]][[s]][[1]]!= IntL[[j]][[s]][[2]]) {
          model.CS <- join_stages(model.CS,paste("V",(j+2),sep = ""),IntL[[j]][[s]][[1]],IntL[[j]][[s]][[r]])
        }
      }
    }
  }
  for (j in c(2:(p+1))) {
    U <- as.list(unique(stages(model.CS,paste("V",j,sep = ""))))
    U <- U[as.integer(U)<=2^(j-2)]
    for (i in U) {
      if (is.element(i,I[[j-1]])==FALSE) {
        model.CS <- join_stages(model.CS,paste("V",(j),sep = ""),i,toString((as.integer(i)+(2^(j-2)))))
      }
    }
  }
  return(model.CS)
}


# Example on three variables
Parts <- partitions(3)
TestL <- Parts[[4]]
ObsTree <- toCStree(3,2,TestL)
S2 <- unique(stages(ObsTree,paste("V",2,sep= "")))
S3 <- unique(stages(ObsTree,paste("V",3,sep= "")))
# The first entry in I should be either list() or list(1). For k>1, the k-th entry in I should be a 
# a subset of Sk, the stages in level k.
S2
S3

# Choose I to be the set of targets
I <- list(list(1),list(1),list(1,4))

#This produces the interventional CStree
TestIntTree <- toInterventionalCStree(3,2,TestL,I)
plot(TestIntTree)
summary(TestIntTree)

# Example on four variables
Parts <- partitions(4)
TestL <- Parts[[1000]]
ObsTree <- toCStree(4,2,TestL)
S2 <- unique(stages(ObsTree,paste("V",2,sep= "")))
S3 <- unique(stages(ObsTree,paste("V",3,sep= "")))
S4 <- unique(stages(ObsTree,paste("V",4,sep= "")))
# The first entry in I should be either list() or list(1). For k>1, the k-th entry in I should be a 
# a subset of Sk, the stages in level k.
S2
S3
S4

# Choose I to be the set of targets
I <- list(list(1),list(1),list(3),list(3,5))

#This produces the interventional CStree
TestIntTree <- toInterventionalCStree(4,2,TestL,I)
plot(TestIntTree)
summary(TestIntTree)

# Example on five variables
TestL <- list(OneDpartitions[[2]],TwoDpartitions[[6]],ThreeDpartitions[[130]],FourDpartitions[[50000]])
ObsTree <- toCStree(5,2,TestL)
S2 <- unique(stages(ObsTree,paste("V",2,sep= "")))
S3 <- unique(stages(ObsTree,paste("V",3,sep= "")))
S4 <- unique(stages(ObsTree,paste("V",4,sep= "")))
S5 <- unique(stages(ObsTree,paste("V",5,sep= "")))
# The first entry in I should be either list() or list(1). For k>1, the k-th entry in I should be a 
# a subset of Sk, the stages in level k.
S2
S3
S4
S5

# Choose I to be the set of targets
I <- list(list(),list(),list(2),list(4,7),list(3,5,13))

#This produces the interventional CStree
TestIntTree <- toInterventionalCStree(5,2,TestL,I)
plot(TestIntTree)
summary(TestIntTree)

# A function for producing all interventional CStrees with one intervention target over all possible targets.
# Inputs are number of variables p, size of state spaces d, a list of stages for the observational CStree L,
# and a list of sets from which the targets in level k should be drawn from element k of the list.
interventionalCStrees <- function(p,d,L,S) {
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
    Trees <- c(Trees,list(toInterventionalCStree(p,2,L,I)))
  }
  return(Trees)
}

# Example on three variables.
# First construct an observed CStree
Parts <- partitions(3)
TestL <- Parts[[4]]
ObsTree <- toCStree(3,2,TestL)
plot(ObsTree)

# Then construct the set S:
S2 <- unique(stages(ObsTree,paste("V",2,sep= "")))
S3 <- unique(stages(ObsTree,paste("V",3,sep= "")))
S <- list(list(1),S2,S3)
S

# Then construct all interventional CStrees where ObsTree is the observational tree and there is one other 
# intervention target. All possible subsets for the stages of ObsTree are considered for the second target.
TestTrees <- interventionalCStrees(3,2,TestL,S)
length(TestTrees)
plot(TestTrees[[1]])
plot(TestTrees[[40]])


# Example on five variables.
# First construct an observed CStree
TestL <- list(OneDpartitions[[2]],TwoDpartitions[[6]],ThreeDpartitions[[50]],FourDpartitions[[16]])
ObsTree <- toCStree(5,2,TestL)
plot(ObsTree)

# Then construct the set S:
S2 <- unique(stages(ObsTree,paste("V",2,sep= "")))
S3 <- unique(stages(ObsTree,paste("V",3,sep= "")))
S4 <- unique(stages(ObsTree,paste("V",4,sep= "")))
S5 <- unique(stages(ObsTree,paste("V",5,sep= "")))
S <- list(list(),S2,S3,S4,S5)
S

# Then construct all interventional CStrees where ObsTree is the observational tree and there is one other 
# intervention target. All possible subsets for the stages of ObsTree are considered for the second target.
TestTrees <- interventionalCStrees(5,2,TestL,S)
length(TestTrees) 
plot(TestTrees[[1]])
plot(TestTrees[[2871]])
summary(TestTrees[[1]])
summary(TestTrees[[2871]])



###########################################


# The following function is not used.  Instead a version of it that is used appears in the file mice.R.

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








