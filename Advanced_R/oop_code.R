library(reshape2)
## Load any other packages that you may need to execute your code
data <- read_csv("data/MIE.csv")

#The variables in the dataset are:
#*id*: the subject identification number
#*visit*: the visit number which can be 0, 1, or 2
#*room*: the room in which the monitor was placed
#*value*: the level of pollution in micrograms per cubic meter
#*timepoint*: the time point of the monitor value for a given visit/room

#You will need to design a class called “LongitudinalData” that characterizes the structure of this longitudinal dataset. You will also need to design classes to represent the concept of a “subject”, a “visit”, and a “room”.

#creating classes
setClass("LongitudinalData", "data.frame")
setClass("subject","data.frame")
setClass("visit", "data.frame")
setClass("room", "data.frame")

#In addition you will need to implement the following functions
#make_LD: a function that converts a data frame into a “LongitudinalData” object
#subject: a generic function for extracting subject-specific information
#visit: a generic function for extracting visit-specific information
#room: a generic function for extracting room-specific information

#set genertic methods
setGeneric("make_LD", function(x){standardGeneric("make_LD")})
setGeneric("subject", function(x,...){standardGeneric("subject")})
setGeneric("summary", function(x){standardGeneric("summary")})
setGeneric("print", function(x){standardGeneric("print")})
setGeneric("visit", function(x,...){standardGeneric("visit")})
setGeneric("room", function(x,...){standardGeneric("room")})

#function that converts a data frame into a “LongitudinalData” object, make "tbl_df" for tibble
setMethod("make_LD", "data.frame", function(x){
	#df <- tibble("ID"=x$id, "Visit"=x$visit, "Room"=x$room, "Value"=x$value, "Timepoint"=x$timepoint)
	df <- data.frame("ID"=x@.Data[[1]], "Visit"=x@.Data[[2]], "Room"=x@.Data[[3]], "Value"=x@.Data[[4]], "Timepoint"=x@.Data[[5]])
	df <- new("LongitudinalData", df)
		return(df)})

#function that converts a data frame into a "subject" object
setMethod("subject", "LongitudinalData", function(x,num){
	#x <- as.data.frame(x)
	df <- new("subject", x[x$ID==num,])
	return(df)})

#function that converts a subject object into a "summary" object
setMethod("summary", "subject", function(x){
	df <- new("summary", x)
	return(df)})

#function that converts a subject object into a "summary" object
setMethod("summary", "room", function(x){
	#x <- base::summary(x)
	df <- new("summary", x)
	return(df)})

#function that converts a data frame into a "vist" object
setMethod("visit", "subject", function(x,num){
	df <- new("visit", x[x$Visit==num,])
	return(df)})

#function that converts a data frame into a "vist" object
setMethod("room", "visit", function(x,st){
	df <- new("room", x[x$Room==st,])
	return(df)})

#make print statement for class
#setMethod("print", class, function(x){class(x)[1]})

#make method to print for LD object
setMethod("print","LongitudinalData", function(x){paste("Longitudinal data with", length(unique(as.data.frame(x)$ID)), "subjects")})

#print function for subject
setMethod("print","subject", function(x){
	if(length(unique(x$ID))==0){
		return(NULL)
	}
	else{
		paste("Subject ID:", unique(x$ID))}})

#print function for summary
setMethod("print", "summary", function(x){
	cat("ID:", unique(x$ID),"\n")
	x <- subset(x, select=-c(ID, Timepoint))
	if (length(unique(x$Room))==1){
		base::summary(x$Value)
	}
	else{
		dcast(data = x, formula = Visit ~ Room, fun.aggregate = mean, value.var = "Value")}})

#print function for summary
setMethod("print", "visit", function(x){
	cat("ID:", unique(x$ID),"\n")
	cat("Visit", unique(x$Visit) ,"\n")})

#print function for summary
setMethod("print", "room", function(x){
	cat("ID:", unique(x$ID),"\n")
	cat("Visit:", unique(x$Visit) ,"\n")
	cat("Room:", unique(x$Room) ,"\n")})

#For each generic/class combination you will need to implement a method, although not all combinations are necessary (see below). You will also need to write print and summary methods for some classes (again, see below).
#To complete this Part, you can use either the S3 system, the S4 system, or the reference class system to implement the necessary functions. 
#For this assessment, you will need to implement the necessary functions to be able to execute the code in the following script file:

x <- make_LD(data)
print(class(x))
print(x)

## Subject 10 doesn't exist
out <- subject(x, 10)
print(out)

out <- subject(x, 14)
print(out)

out <- subject(x, 54) %>% summary
print(out)

out <- subject(x, 14) %>% summary
print(out)

out <- subject(x, 44) %>% visit(0) %>% room("bedroom")
print(out)

## Show a summary of the pollutant values
out <- subject(x, 44) %>% visit(0) %>% room("bedroom") %>% summary
print(out)

out <- subject(x, 44) %>% visit(1) %>% room("living room") %>% summary
print(out)
