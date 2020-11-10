library("stagedtrees")
library("combinat")


######### Producing the possible stagings for each level in a CStree.  
######### (For higher dimensions we will need Eliana's code.)
OneDpartitions <- list(list(list(1,1)),
                       list(list(1,2)))
TwoDpartitions <- list(list(list(1,1)),
                       list(list(1,2)),
                       list(list(3,4)),
                       list(list(1,3)),
                       list(list(2,4)),
                       list(list(1,3),list(2,4)),
                       list(list(1,2),list(3,4)),
                       list(list(1,2,3,4)))
# start of creating 3dpartitions, reading from file Note: load the functions f,g,h from the Read3DCStrees.R
z<-readLines(con = "threeDpartitions.txt")
listOfchars <- lapply(z, h)
ThreeDpartitions <-lapply(listOfchars,f)
# start of creating 4dpartitions, reading from file

##### A list of lists where the i-th entry is the partitions of the i-dimensional cube.
ComputedPartitions = list(OneDpartitions,TwoDpartitions,TwoDpartitions,TwoDpartitions)

##### A function for computing the cross product of partitions of the cubes of dimensions 1,...,p-1.  
##### The function calls upon the set ComputedPartitions and can only handle values for p up to its cardinality plus 1.
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

#### A function that takes in a number of variables p , a size of outcome space of the variables d, 
#### and a list of partitions L, one for each level of the tree, ordered by level from smallest-to-largest level.  
#### It returns a CStree with the staging corresponding to the partitions given in L and causal ordering the
#### identity permutation of [p].
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

#### A function that produces all CStrees on p variables with each variable having d outcomes.  
#### It requires that you have computed the set of all partitions of the i-dimensional cube for i< p-1, and added
#### them to the list ComputedPartitions.
#### For now, d = 2, as the set ComputedPartitions is defined for only binary variables.  
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
      for ( j in c(1:(p-1))) {
        if (class(L[[j]][[1]]) == "numeric") {
          if (L[[j]][[1]] != L[[j]][[2]]) {
            for (r in c(2:length(L[[j]]))) {
              k<-j+1
              M <- join_stages(M,paste("V",Perms[[i]][[k]],sep = ""),L[[j]][[1]],L[[j]][[r]])
            }
          }
        } else {
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
      }
      models <-c(models,list(M))
    }
  }
  return(models)
}



##### Example:
plot(CStrees(3,2)[[1]])
summary(CStrees(3,2)[[1]])


##### Application to Real Data Example:
library("bnlearn")

data(lizards) #dataset available within the bnlearn library
myData <- lizards
names(myData) <- c("V1","V2","V3") #This relabels the variables in the dataset with the variable names produced by CStrees(p,2).
levels(myData$V1) <- c(1:2) #This relabels the outcomes to match the outcome names used by CStrees(p,2).
levels(myData$V2) <- c(1:2)
levels(myData$V3) <- c(1:2)
CStreesList <- CStrees(3,2)

# The following will recover the BIC optimal CStree for the dataset myData:
M <- CStreesList[[1]]
M.fit <- sevt_fit(M,data=myData,lambda=1)
t <- BIC(M.fit)
MM <- M.fit
summary(M.fit)
optTrees = list(1)
for (i in c(1:length(CStreesList))) {
  N.fit <- sevt_fit(CStreesList[[i]],data=myData,lambda=1)
  s <- BIC(N.fit)
  if (t < s) {
    MM <- N.fit
    t <- s
    optTrees <- list(i)
  } else {
    if (t==s) {
      optTrees <- c(optTrees,i)
    }
  }
}
# The first BIC optimal tree in the list CStreesList is plotted, it and described below.  t is its BIC score. 
plot(MM)
summary(MM)
t

# The following set contains the indices i for which CStreesList[[i]] is a BIC optimal tree w.r.t. the myData dataset.
optTrees

# The original variable names and outcomes for each variable are given by the following:
names(lizards)
levels(lizards$Species)
levels(lizards$Diameter)
levels(lizards$Height)





