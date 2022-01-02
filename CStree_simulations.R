library("stagedtrees")
library("purrr")
library("useful")
library("DescTools")


# Identify the stage-defining context of a stage s in the level for variable k of a CStree t.  
# Returns a dataframe with one row specifying the minimal context and column entries the context
# variables.
defContext <- function(t,k,s) { 
  p <- get_path(t,k,s) #returns a data frame of paths to this stage
  contextVars <- c()
  for (j in colnames(p)) {
    if (length(unique(p[,j]))==1) {
      contextVars <- append(contextVars,j)
    }
  } 
  defCon <- as.data.frame(p[,contextVars])
  colnames(defCon) <- contextVars
  rownames(defCon) <- 1:length(rownames(defCon))
  return(head(defCon,1))
}

# A function that takes in two stages s and r in the level associated to variable k in a CStree t
# and returns a the CStree where all stages containing a node that appears in the stage defined by
# the context given by the stage union of s and r have been merged into a single stage.
CSstageJoin <- function(t,k,s,r) {
  overlapStages <- c(r)
  tt <- t
  while (length(overlapStages) != 0) {
    tt <- join_stages(tt,k,s,overlapStages[1])
    Stages <- unique(tt$stages[[k]])
    Stages <- Stages[Stages!=s]
    SU <- defContext(tt,k,s)
    nSU <- colnames(SU)
    if (length(nSU)==0) {
      overlapStages <- Stages
    } else {
      overlapStages <- overlapStages[-1]
      if (length(overlapStages) == 0) {
        for (j in Stages) {
          p <- get_path(tt,k,j)
          p <- as.data.frame(p[,nSU])
          colnames(p) <- nSU
          for (i in rownames(p)) {
            path <- p[i,]
            if (identical(as.list(path),as.list(SU[1,]))==TRUE) {
              overlapStages <- append(overlapStages,j)
            }
          }
        }
        overlapStages <- overlapStages
      } 
    }
    overlapStages <- unique(overlapStages)
  }
  return(tt)
}

# The following function randomly merges stages in a dependence model on n variables to produce a 
# CStree. The probability of merging two stages within the same level is p.  Sparser models have 
# fewer stages, hence a high value of p produces sparse models. The value q is a parameter 
# determining the fraction of stage-merges that can take place in each level, where a higher q means
# fewer possible merges.  Instead of using q, one can also adjust the number of stage-merges either 
# logarithmically or polynomially.  In these cases, run the function with q=1, as it is irrelevant.  
# The polynomial adjustment appears to give the best approximation of the average number of stages 
# decreasing linearly as p increases.  Hence, we use polynomial adjustment in our simulations.  To use
# q adjustment or logarithmic adjustment, just uncomment the appropriate definition of 'merges' and comment
# out the others.
randomCStree <- function(n,p,q) {
  Vars<-list()
  for (j in 1:n) {
    Vars[[j]] <- c(1:2)
  }
  t <- sevt(Vars,full=TRUE)
  for (j in names(t$tree[-1])) {
    Stages <- unique(t$stages[[j]])
    k <- log(length(Stages),2)
    # merges <- rbernoulli(floor((length(Stages)/q)*p),p) #q adjusted bernoulli trials
    # merges <- rbernoulli(floor((length(Stages)/log(n-k+2))*p),p) #logarithmically adjusted
    merges <- rbernoulli(floor(length(Stages)/(1+4*k*(p-p^2))),p) #polynomially adjusted
    for (i in merges) {
      if (i == TRUE) {
        if (length(Stages)!=1) {
          mergingStages <- rdunif(2,1,length(Stages))
          while (mergingStages[1]==mergingStages[2]) {
            mergingStages <- rdunif(2,1,length(Stages))
          }
          t <- CSstageJoin(t,j,Stages[mergingStages[1]],Stages[mergingStages[2]])
          Stages <- unique(t$stages[[j]]) 
        }
      }
    }
  }
  return(t)
}

### Example:
CS4 <- randomCStree(4,3/10,1)
plot(CS4)

# Generates a random stage tree on n levels with merging probability p and merge bound q.
randomStageTree <- function(n,p,q) {
  Vars<-list()
  for (j in 1:n) {
    Vars[[j]] <- c(1:2)
  }
  t <- sevt(Vars,full=TRUE)
  for (j in names(t$tree[-1])) {
    Stages <- unique(t$stages[[j]])
    merges <- rbernoulli(floor(length(Stages)/q),p)
    for (i in merges) {
      if (i == TRUE) {
        if (length(Stages)!=1) {
          mergingStages <- rdunif(2,1,length(Stages))
          while (mergingStages[1]==mergingStages[2]) {
            mergingStages <- rdunif(2,1,length(Stages))
          }
          t <- join_stages(t,j,Stages[mergingStages[1]],Stages[mergingStages[2]])
          Stages <- unique(t$stages[[j]]) 
        }
      }
    }
  }
  return(t)
}

### Example:
ST4 <- randomStageTree(4,3/10,1)
plot(ST4)

# The following function random bernoulli distributions to be assigned to each stage in a 
# binary staged tree t.
probabilities <-function(t) {
  Stages <-list()
  for (j in names(t$stages)) {
    Stages[[j]] <- unique(t$stages[[j]])
  }
  p <- runif(1,0,1)
  m <- names(t$tree)[1]
  Probs <- c()
  Probs[[names(t$tree)[1]]][["1"]] <- c(`1`=p,`2`=1-p) 
  for (j in names(t$stages)) {
    for (s in Stages[j]) {
      for (i in s) {
        p <- runif(1,0,1)
        Probs[[j]][[i]] <- c(`1`=p,`2`=1-p)
      }
    }
  }
  return(Probs)
}

# Generate m binary random CStrees on n variables with probability of merging p and merging bound q.
# Each stage in each tree is randomly assigned a probability distribution via the 'probabilities' function.
randomCStrees <- function(m,n,p,q) {
  L <- c()
  for (i in 1:m) {
    t <- randomCStree(n,p,q)
    t$prob <- probabilities(t)
    L[[i]] <- t
  }
  return(L)
}

### Example:
randCS4 <- randomCStrees(2,4,3/10,1)
plot(randCS4[[1]])
plot(randCS4[[2]])
randCS4[[1]]$prob
randCS4[[2]]$prob

# Generate m binary random staged trees on n variabes with probability of merging p and merging bound
# q with each stage randomly assigned a distribution via the 'probabilities' function.
randomStageTrees <- function(m,n,p,q) {
  L <- c()
  for (i in 1:m) {
    t <- randomStageTree(n,p,q)
    t$prob <- probabilities(t)
    L[[i]] <- t
  }
  return(L)
}

### Example:
randST4 <- randomStageTrees(2,4,3/10,1)
plot(randST4[[1]])
plot(randST4[[2]])
randST4[[1]]$prob
randST4[[2]]$prob

# We then want to generate data sets of varying sizes from each of our randomly generated trees.

# Generates m random staged trees on n variables with merging probability q for each merging probability 
# p = i/10 where i ranges from 1 to 10.  For each tree, we produce a dataset of 10, 1,000, and 10,000 
# samples.
randomStreeDataMult <- function(m,n,q) {
  tsequence <- c()
  for (i in 1:10) {
    tsequence[[i]] <- randomStageTrees(m,n,i/10,q)
  }
  datasets <- c()
  for (i in 1:10) {
    di <- c()
    for (j in 1:m) {
      sj <- sample_from(tsequence[[i]][[j]],nsim=10000,seed=NULL)
      sj1000 <- head(sj,1000)
      di[[j]] <-list(sj,
                     sample_from(tsequence[[i]][[j]],nsim=10,seed=NULL),
                     sj1000)
    }
    datasets[[i]] <- di
  }
  R <- c()
  for (i in 1:10) {
    Ri <- c()
    for (j in 1:m) {
      Ri[[j]] <- list(tsequence[[i]][[j]],datasets[[i]][[j]])
    }
    R[[i]] <-Ri
  }
  return(R)
}

### Example:
randST4mult <- randomStreeDataMult(2,4,1)
length(randST4mult) #There is one element in the output for every p = i/10 for i from 1 to 10.
plot(randST4mult[[2]][[1]][[1]])
plot(randST4mult[[2]][[2]][[1]])
head(randST4mult[[7]][[1]][[2]][[2]])
nrow(randST4mult[[7]][[1]][[2]][[1]])
nrow(randST4mult[[7]][[1]][[2]][[2]])
nrow(randST4mult[[7]][[1]][[2]][[3]])

# Generates m random CStrees on n variables with merging probability q for each merging probability
# p = i/10 where i ranges from 1 to 10.  For each tree, we produce a dataset of 10, 1,000, 10,000, 
# and 100,000 samples.
randomCStreeDataMult <- function(m,n,q) {
  tsequence <- c()
  for (i in 1:10) {
    tsequence[[i]] <- randomCStrees(m,n,i/10,q)
  }
  datasets <- c()
  for (i in 1:10) {
    di <- c()
    for (j in 1:m) {
      sj <- sample_from(tsequence[[i]][[j]],nsim=100000,seed=NULL)
      sj10000 <- head(sj,10000)
      sj1000 <- head(sj,1000)
      di[[j]] <-list(sj10000,
                     sample_from(tsequence[[i]][[j]],nsim=10,seed=NULL),
                     sj1000,
                     sj)
    }
    datasets[[i]] <- di
  }
  R <- c()
  for (i in 1:10) {
    Ri <- c()
    for (j in 1:m) {
      Ri[[j]] <- list(tsequence[[i]][[j]],datasets[[i]][[j]])
    }
    R[[i]] <-Ri
  }
  return(R)
}

### Example:
randCS4mult <- randomCStreeDataMult(2,4,1)
length(randCS4mult) #There is one element in the output for every p = i/10 for i from 1 to 10.
plot(randCS4mult[[5]][[1]][[1]]) #for each i, there are m elements, each a list, the first entry the data-generating tree
plot(randCS4mult[[5]][[2]][[1]])
head(randCS4mult[[7]][[1]][[2]][[2]]) #the second entry contains the data sets, for CStrees there are four.
nrow(randCS4mult[[7]][[1]][[2]][[1]]) #nrow() shows which data set contains how many samples.
nrow(randCS4mult[[7]][[1]][[2]][[2]])
nrow(randCS4mult[[7]][[1]][[2]][[3]])
nrow(randCS4mult[[7]][[1]][[2]][[4]]) #note we have one more data set as compared to the previous function.

# We would then like an algorithm for learning a BIC-optimal CStree with a fixed causal ordering given a
# dataset, as all of our randomly generated trees have the same causal ordering:

# An algorithm that, given the dataframe d, goes through each level, iteratively merging stages according 
# to which merge is BIC-optimal given the data, and only moves to the next level when no merge that 
# decreases BIC can be found.  Each merge is done to produce a CStree making this a CStree analogue of the 
# backward hill-climbing algorithm for staged trees available in the stagedtrees R package.  
CStree_bhc <- function(d) {
  t <- sevt_fit(full(d))
  N <- names(t$tree)
  b <- BIC(t)
  for (j in N[-1]) {
    repeat {
      B <- b
      print(j)
      leveltrees <- c()
      if (length(unique(t$stages[[j]]))>1) {
        stagepairs <- as.matrix(combn(unique(t$stages[[j]]),2))
        for (i in 1:dim(stagepairs)[2]) {
          tm <- sevt_fit(CSstageJoin(t,j,stagepairs[1,i],stagepairs[2,i]),d)
          leveltrees[[i]] <- tm
        }
      }
      for (Tm in leveltrees) {
        Bm <- BIC(Tm)
        if (Bm < b) {
          t <- Tm
          b <- Bm
        }
      }
      if (B == b) {
        break
      }
    }
  }
  return(t)
}

# Given a staged tree and a CStree learned from the same data set, we would like to compare their predictive
# accuracies.

# Computing the predictive accuracy of a learned tree t based on a validation data set D.  
predictiveAccuracy <- function(t,D) {
  n <- nrow(D)
  N <- names(t$tree)
  trues <- c()
  for (i in N) {
    trues <- append(trues,length(which((D[,i] == predict(t,D,class=i))==TRUE)))
  }
  return(sum(trues)/(n*length(N)))
}

# Comparing models for multiple datasets from same trees.  Data sets and trees, D, are input as 
# the output of the previous function along with its parameter values, m,n,q.  The parameter k
# determines how many samples you consider: for the output of the randomCStreeDataMult function k =1,2,3,4 
# corresponds to 10 000, 10, 1 000, and 100 000, and samples, respectively. (see above example.)
averagePredAccuraciesMult <- function(D,k,m,n,q) {
  learnedStrees <- c()
  for (i in 1:10) {
    li <- c()
    for (j in 1:m) {
      li[[j]] <- stages_bhc(full(D[[i]][[j]][[2]][[k]]))
    }
    learnedStrees[[i]] <- li
  }
  averagesStrees <- c()
  for (i in 1:10) {
    accuracies <- c()
    for (j in 1:m) {
      accuracies[[j]] <- predictiveAccuracy(learnedStrees[[i]][[j]],D[[i]][[j]][[2]][[2]])
    }
    averagesStrees[[i]] <- sum(accuracies)/m
  }
  learnedCStrees <- c()
  for (i in 1:10) {
    li <- c()
    for (j in 1:m) {
      li[[j]] <- CStree_bhc(D[[i]][[j]][[2]][[k]])  #stages_bhc(full(D[[i]][[j]][[1]]))
      print(i)
      print(j)
    }
    learnedCStrees[[i]] <- li
  }
  averagesCStrees <- c()
  for (i in 1:10) {
    accuracies <- c()
    for (j in 1:m) {
      accuracies[[j]] <- predictiveAccuracy(learnedCStrees[[i]][[j]],D[[i]][[j]][[2]][[2]])
    }
    averagesCStrees[[i]] <- sum(accuracies)/m
  }
  R <- c()
  R[[1]] <- D
  R[[2]] <- learnedStrees
  R[[3]] <- learnedCStrees
  R[[4]] <- averagesStrees
  R[[5]] <- averagesCStrees
  return(R)
}

### Example:
# Computing the average predictive accuracies based on 10,000 samples (hence, k=2)
aPA_CS_2_4_1_10000 <- averagePredAccuraciesMult(randCS4mult,1,2,4,1) 
# aPA_CS_10_4_1_10000[[1]] returns the input D
# aPA_CS_10_4_1_10000[[2]] returns the list of learned staged trees via backwards hill climbing (stages_bhc function in stagedtrees package)
# aPA_CS_10_4_1_10000[[3]] returns the list of learned CStrees via CStree_bhc
aPA_CS_2_4_1_10000[[4]] #list of average predictive accuracies on the 10 validation samples for the learned staged trees with p = i/10, i going from 1 to 10.
aPA_CS_2_4_1_10000[[5]] #list of average predictive accuracies on the 10 validation samples for the learned CStrees with p = i/10, i going from 1 to 10.
# We can plot the results for both the learned staged trees and CStrees to compare:
x <- c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1)
plot(x,aPA_CS_2_4_1_10000[[4]],type = "b",frame=FALSE,pch = 19, col = "red", 
     xlab = "Probability of Merging (p)", ylab = "Average Predictive Accuracy", 
     xlim = c(0,1),ylim = c(0,1))
lines(x, aPA_CS_2_4_1_10000[[5]], pch = 18, col = "blue", type = "b", lty = 2,ylim=c(1,100))
legend("bottomright", legend=c("Staged Trees", "CStrees"),
       col=c("red", "blue"), lty = 1:2, pch = c(19,18), cex=0.8)



# Computing the average number of stages of the true data-generating trees:

# Computes the number of stages in a given stage tree
numberStages <- function(t) {
  L <- c()
  for (j in names(t$tree)) {
    L<-append(L,length(unique(t$stages[[j]])))
  }
  return(sum(L))
}

# Computes the average number of stages for a given list of random stage trees:
averageStages <- function(L) {
  return((1/length(L))*sum(sapply(L,numberStages)))
}

# Given D, the output of randomCStreeDataMult or randomStreeDataMult run with parameters m and n,
# returns the average number of stages in the generated trees for each p = i/10.
stageTestData <- function(D,m,n) {
  L <- c()
  for (i in 1:10) {
    Li <- c()
    for (j in 1:m) {
      Li[[j]] <- D[[i]][[j]][[1]]
    }
    L[[i]] <- averageStages(Li)
  }
  return(L)
}

### Example: 
avgNoStages2_4_1 <- stageTestData(randCS4mult,2,4)
avgNoStages2_4_1
#Plotting the results:
x <- c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1)
plot(x,avgNoStages2_4_1,type = "b",frame=FALSE,pch = 19, col = "red", 
     xlab = "Probability of Merging (p)", ylab = " Average Number of Stages", 
     xlim = c(0,1),ylim = c(0,16))


# Computing the average structural Hamming Distance:

# Average Structural Hamming Distances: Return five lists based on the output of 
# averagePredAccuracies, the average structural Hamming distance between the data-generating tree 
# and the learned staged tree, the same values for the data-generating tree and the learned CStree, 
# the same values for the learned staged tree and learned CStree, the proportion of learned CStrees with
# structural Hamming distance 0 from the true tree, and the structural Hamming distances of each learned
# CStree from the true tree.  The output is a list, and each of the above is the i-th element of that list
# in the given order.
averageStructuralHamming <- function(r,m,n,q) {
  ham_true_s <- c()
  for (i in 1:10) {
    tsi <- c()
    for (j in 1:m) {
      tsi[[j]] <- hamming_stages(r[[1]][[i]][[j]][[1]],r[[2]][[i]][[j]])
    }
    ham_true_s[[i]] <- tsi
  }
  ham_true_s_avg <- c()
  for (i in 1:10) {
    ham_true_s_avg[[i]] <- sum(ham_true_s[[i]])/m
  }
  ham_true_cs <- c()
  for (i in 1:10) {
    tsi <- c()
    for (j in 1:m) {
      tsi[[j]] <- hamming_stages(r[[1]][[i]][[j]][[1]],r[[3]][[i]][[j]])
    }
    ham_true_cs[[i]] <- tsi
  }
  ham_true_cs_avg <- c()
  for (i in 1:10) {
    ham_true_cs_avg[[i]] <- sum(ham_true_cs[[i]])/m
  }
  ham_true_cs_prop <- c()
  for (i in 1:10) {
    si <- 0
    for (j in 1:m) {
      if (ham_true_cs[[i]][[j]]==0) {
        si <- si+1
      }
    }
    ham_true_cs_prop[[i]] <- si
  }
  ham_s_cs <- c()
  for (i in 1:10) {
    tsi <- c()
    for (j in 1:m) {
      tsi[[j]] <- hamming_stages(r[[2]][[i]][[j]],r[[3]][[i]][[j]])
    }
    ham_s_cs[[i]] <- tsi
  }
  ham_s_cs_avg <- c()
  for (i in 1:10) {
    ham_s_cs_avg[[i]] <- sum(ham_s_cs[[i]])/m
  }
  return(list(ham_true_s_avg,ham_true_cs_avg,ham_s_cs_avg,ham_true_cs_prop,ham_true_cs))
}

### Example:
sHD2_4_1 <- averageStructuralHamming(aPA_CS_2_4_1_10000,2,4,1)
sHD2_4_1[[1]] #average structural Hamming distances from the learned staged trees to the true trees for p = i/10 for i from 1 to 10.
sHD2_4_1[[2]] #average structural Hamming distances from the CStrees to the true trees for p = i/10 for i from 1 to 10.
#We can again plot the results to compare:
x <- c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1)
plot(x,sHD2_4_1[[1]],type = "b",frame=FALSE,pch = 19, col = "red", 
     xlab = "Probability of Merging (p)", ylab = "Average Structural Hamming Distance", 
     xlim = c(0,1),ylim = c(0,10))
lines(x, sHD2_4_1[[2]], pch = 18, col = "blue", type = "b", lty = 2,ylim=c(1,10))
legend("topright", legend=c("True vs. Staged Trees", "True vs. CStrees"),
       col=c("red", "blue"), lty = 1:2, pch =c(19,18), cex=0.8)


# Using the last entry in the output of averageStructuralHamming, we can also compare the number 
# of stages of the true tree to the structural Hamming distance of the learned tree from the true tree 
# for each learned CStree.  The input r is the output of averagePredAccuraciesMult with parameters m,n,q.
NoStageVsHD <- function(r,m,n,q) {
  sHD <- averageStructuralHamming(r,m,n,q)[[5]]
  NSvSHD <- c()
  HD <- c()
  NS <- c()
  for (i in 1:10) {
    for (j in 1:m) {
      NS[[length(NS)+1]] <- numberStages(r[[1]][[i]][[j]][[1]])
      HD[[length(HD)+1]] <- sHD[[i]][[j]]
    }
    
  }
  NSvSHD[[1]] <- NS
  NSvSHD[[2]] <- HD
  return(NSvSHD)
}

### Example: 
NSvHD_CS_2_4_1_10000 <- NoStageVsHD(aPA_CS_2_4_1_10000,2,4,1)
plot(NSvHD_CS_2_4_1_10000[[1]], NSvHD_CS_2_4_1_10000[[2]], col = "red",
     xlab="Number of Stages", ylab="Structural Hamming Distance", pch=19)
legend("topleft", legend=c("10 000 samples"),
       col=c("red"), pch = c(19), cex=0.8)


##########################################################################################################
# The code run on n=7 variables to produce the average predictive accuracy comparisons of learned staged 
# trees with learned CStrees.  Running the first line 'Data_S_10_7_1 <- randomStreeDataMult(10,7,1)' will
# produce a new random data set with the same parameters.  To work with the data set generated for the
# paper uncomment the second line, enter in the directory to which you have downloaded the file
# simulations_S_10_7_1.R and run it.
Data_S_10_7_1 <- randomStreeDataMult(10,7,1)
# load("/PATH_TO_WORKING_DIRECTORY/simulations_S_10_7_1.R") #loaded R object should already have name 'Data_S_10_7_1'

# Producing the average predictive accuracy plots:
aPA10_7_1_c <- averagePredAccuraciesMult(Data_S_10_7_1,1,10,7,1)
x <- c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1)
plot(x,aPA10_7_1_c[[4]],type = "b",frame=FALSE,pch = 19, col = "red", 
     xlab = "Probability of Merging (p)", ylab = "Average Predictive Accuracy", 
     xlim = c(0,1),ylim = c(0,1))
lines(x, aPA10_7_1_c[[5]], pch = 18, col = "blue", type = "b", lty = 2,ylim=c(1,100))
legend("bottomright", legend=c("Staged Trees", "CStrees"),
       col=c("red", "blue"), lty = 1:2, pch = c(19,18), cex=0.8)

# Producing the average structural Hamming distance plot:
sHD10_7_1_c <- averageStructuralHamming(aPA10_7_1_c,10,7,1)
plot(x,sHD10_7_1_c[[1]],type = "b",frame=FALSE,pch = 19, col = "red", 
     xlab = "Probability of Merging (p)", ylab = "Average Structural Hamming Distance", 
     xlim = c(0,1),ylim = c(0,100))
lines(x, sHD10_7_1_c[[2]], pch = 18, col = "blue", type = "b", lty = 2,ylim=c(1,100))
legend("topright", legend=c("True vs. Staged Trees", "True vs. CStrees"),
       col=c("red", "blue"), lty = 1:2, pch =c(19,18), cex=0.8)

# Producing the average number of stages for the true trees plot:
avgNoStages10_7_1_c <- stageTestData(Data_S_10_7_1,10,7)
plot(x,avgNoStages10_7_1_c,type = "b",frame=FALSE,pch = 19, col = "red", 
     xlab = "Probability of Merging (p)", ylab = " Average Number of Stages", 
     xlim = c(0,1),ylim = c(0,130))
##########################################################################################################


##########################################################################################################
# The code run on n=6 variables to produce the analysis of the BHC-CS CStree learning algorithm given 
# different sample sizes.  Running the first line 'CStestData10_6_1_c <- randomCStreeDataMult(10,6,1)' will
# produce a new random data set with the same parameters.  To work with the data set generated for the
# paper uncomment the second line, enter in the directory to which you have downloaded the file
# simulations_CS_10_6_1.R and run it.

#CS test data constructed with polynomially adjusted number of bernoulli trials (corrected)
CStestData10_6_1_c <- randomCStreeDataMult(10,6,1)
# load("/PATH_TO_WORKING_DIRECTORY/simulations_CS_10_6_1.R") #loaded R object should already have name 'CStestData10_6_1_c'

# Producing the average number of stages of the true trees plot:
avgNoStagesCS10_6_1_c <- stageTestData(CStestData10_6_1_c,10,6)
x <- c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1)
plot(x,avgNoStagesCS10_6_1_c,type = "b",frame=FALSE,pch = 19, col = "red", 
     xlab = "Probability of Merging (p)", ylab = "Average Number of Stages", 
     xlim = c(0,1),ylim = c(0,65))

# Learning CStrees for the data sets of sizes 1000, 10000, and 100000, respectively:
aPACS10_6_1_1000_c <- averagePredAccuraciesMult(CStestData10_6_1_c,3,10,6,1)
aPACS10_6_1_10000_c <- averagePredAccuraciesMult(CStestData10_6_1_c,1,10,6,1)
aPACS10_6_1_100000_c <- averagePredAccuraciesMult(CStestData10_6_1_c,4,10,6,1)

# Producing the average structural Hamming distance plot:
sHDCS10_6_1_1000_c <- averageStructuralHamming(aPACS10_6_1_1000_c,10,6,1)
sHDCS10_6_1_10000_c <- averageStructuralHamming(aPACS10_6_1_10000_c,10,6,1)
sHDCS10_6_1_100000_c <- averageStructuralHamming(aPACS10_6_1_100000_c,10,6,1)
plot(x,sHDCS10_6_1_1000_c[[2]],type = "b",frame=FALSE,pch = 19, col = "red", 
     xlab = "Probability of Merging (p)", ylab = "Average Structural Hamming Distance", 
     xlim = c(0,1),ylim = c(0,60))
lines(x, sHDCS10_6_1_10000_c[[2]], pch = 18, col = "blue", type = "b", lty = 2,ylim=c(1,100))
lines(x, sHDCS10_6_1_100000_c[[2]], pch = 17, col = "green", type = "b", lty = 2,ylim=c(1,100))
legend("topright", legend=c("1 000 samples", "10 000 samples", "100 000 samples"),
       col=c("red", "blue", "green"), lty = 1:3, pch = c(19,18,17), cex=0.8)

# Producing the plot comparing the true number of stages with the structural Hamming distance of the 
# learned tree to the true tree:
NSvHDCS10_6_1_1000_c <- NoStageVsHD(aPACS10_6_1_1000_c,10,6,1)
NSvHDCS10_6_1_10000_c <- NoStageVsHD(aPACS10_6_1_10000_c,10,6,1)
NSvHDCS10_6_1_100000_c <- NoStageVsHD(aPACS10_6_1_100000_c,10,6,1)
plot(NSvHDCS10_6_1_1000_c[[1]], NSvHDCS10_6_1_1000_c[[2]], col = "red",
     xlab="Number of Stages", ylab="Structural Hamming Distance", pch=19)
points(NSvHDCS10_6_1_10000_c[[1]], NSvHDCS10_6_1_10000_c[[2]], col = "blue", pch=18)
points(NSvHDCS10_6_1_100000_c[[1]], NSvHDCS10_6_1_100000_c[[2]], col = "green", pch=17)
legend("topleft", legend=c("1 000 samples", "10 000 samples", "100 000 samples"),
       col=c("red", "blue", "green"), pch = c(19,18,17), cex=0.8)
##########################################################################################################


# The previous learning function only produces a BIC-optimal tree with a fixed causal ordering determined
# by the order of the variables in the dataframe.  To learn a BIC-optimal tree with unknown causal ordering
# we write a function that considers all possible causal orderings (BHC-CS-perm):

# A function that takes in a dataframe D and returns the optimal tree given by 
# CS_bhc for each possible causal ordering of the variables.  The output is a list in which the
# first element is the causal ordering of the variables in the dataframe and the second is the 
# optimal CStree learned by CStree_bhc for that causal ordering.  
CStree_bhc_perm <- function(D) {
  orignames <- names(D)
  n <- length(names(D))
  Vars<-list()
  for (j in 1:n) {
    Vars[[j]] <- c(1:2)
  }
  t <- sevt(Vars,full=TRUE)
  P <- Permn(names(t$tree))
  p <- nrow(P)
  names(D) <- P[1,]
  Trees <- c()
  for (i in 1:p) {
    d <- D[,P[i,]]
    t <- CStree_bhc(d)
    renames <- c()
    for (k in 1:n) {
      renames[[k]] <- orignames[[strtoi(substr(P[i,k],2,2))]]
    }
    Treesi <- list(renames,t)
    Trees[[i]] <- Treesi
  }
  learnedT <- Trees[[1]]
  b <- BIC(Trees[[1]][[2]])
  for (L in Trees) {
    bb <- BIC(L[[2]])
    if (b < bb) {
      learnedT <- L
      b <- bb
    }
  }
  return(learnedT)
}

# Example with the lizards data set from bnlearn:
library("bnlearn")
data(lizards)
optLiz <- CStree_bhc_perm(lizards) # the output is a list of two elements: the first the causal ordering of the variables, the second the CStree.
plot(optLiz[[2]]) # plot the learned CStree
optLiz[[1]] # see which variables correspond to which level of the tree.

#Example with the coronary data set from bnlearn:
library("bnlearn")
data(coronary)
optCoronary <- CStree_bhc_perm(coronary)
plot(optCoronary[[2]])
optCoronary[[1]]

############