myCortexdf <- read.xls("Data_Cortex_Nuclear.xls", sheet = 1, header = TRUE)

# Mice data
#--------------------------------------------------------------
summary(myCortexdf)
names(myCortexdf)
summary(myCortexdf$class)
names(myObsCortex)
myObsCortex<- myCortexdf[which(myCortexdf$class== "c-SC-m"),]
myObsCortex<- myCortexdf[which(myCortexdf$class== "c-SC-s"),]
summary(myObsCortex)
# We select the variables according to column two of Table 3
myObsCortex<- subset(myObsCortex,select = c(AKT_N,pNUMB_N,pNR1_N,pCAMKII_N))
# Then we discretize by quantile
myObsCortex$AKT_Ncat <- cut(myObsCortex$AKT_N,quantile(myObsCortex$AKT_N, probs = c(0,0.5,1)),include.lowest =TRUE)
myObsCortex$pNUMB_Ncat <- cut(myObsCortex$pNUMB_N,quantile(myObsCortex$pNUMB_N, probs = c(0,0.5,1)),include.lowest =TRUE)
myObsCortex$pNR1_Ncat  <- cut(myObsCortex$pNR1_N,quantile(myObsCortex$pNR1_N, probs = c(0,0.5,1)),include.lowest =TRUE)
myObsCortex$pCAMKII_Ncat  <- cut(myObsCortex$pCAMKII_N,quantile(myObsCortex$pCAMKII_N, probs = c(0,0.5,1)),include.lowest =TRUE)
# After discretize
myObsdata <- data.frame(myObsCortex$pNR1_Ncat,myObsCortex$pNUMB_Ncat,myObsCortex$pNR1_Ncat,myObsCortex$pCAMKII_Ncat)                    
names(myObsdata) <- c("V1","V2","V3","V4") #This relabels the variables in the dataset with the variable names produced by CStrees(p,2).
levels(myObsdata$V1) <- c(1:2) #This relabels the outcomes to match the outcome names used by CStrees(p,2).
levels(myObsdata$V2) <- c(1:2)
levels(myObsdata$V3) <- c(1:2)
levels(myObsdata$V4) <- c(1:2)
CStreesList <- CStrees(4,2)
CStrees(4,2)
M <- CStreesList[[1]]
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
plot(MM)
summary(MM)
t
summary(myObsdata)

?quantile
?cut
x<- cut(myObsCortex$pBRAF_N,quantile(myObsCortex$pBRAF_N, probs = c(0,0.5,1)),include.lowest =TRUE)
x <- quantile(myObsCortex$pBRAF_N, probs = c(0,0.5,1))
x
y<-cut(myObsCortex$pBRAF_N,quantile(myObsCortex$pBRAF_N, probs = c(0,0.5,1)))
summary(x)
nrow(y)
myObsCortex$pBRAF_N[c(100:104)]
typeof(myObsCortex)

nrow(myObsCortex)
#---- preprocessing for interventions.
myIntCortex<- myCortexdf[which(myCortexdf$class== "c-CS-s"),]
myIntCortex<- subset(myIntCortex,select = c(pBRAF_N,pNUMB_N,pNR1_N,pCAMKII_N))
summary(myIntCortex)
nrow(myIntCortex)
#--------------------------------------------------------------