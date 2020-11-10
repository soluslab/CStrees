##

## elements in different lines get assigned a different element
# TO DO later: Turn the code below into a function that reads
# partitions for any dimension d.
### This part creates the 3Dpartitions
z<-readLines(con = "threeDpartitions.txt")
listOfchars <- lapply(z, h)
ThreeDpartitions <-lapply(listOfchars,f)
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
    pieces <-strsplit(mystring," ")
    lapply(pieces, strtoi)
}
#######################3
three3ps

