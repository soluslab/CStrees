##

## elements in different lines get assigned a different element
# TO DO later: Turn the code below into a function that reads
# partitions for any dimension d.
### This part creates the 3Dpartitions
z<-readLines(con = "threeDpartitions.txt")
listOfchars <- lapply(z, h)
ThreeDpartitions <-lapply(listOfchars,f)
ThreeDpartitions <-lapply(ThreeDpartitions,k)
length(ThreeDpartitions)
# Takes a string of characters and breaks it down into
# a list using the delimiter ,
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
#######################3
three3ps

as.list(listOfchars[[2]][[1]])


