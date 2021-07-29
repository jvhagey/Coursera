library(purrr)
library(microbenchmark)

#Part 1: Factorial Function
#The objective of Part 1 is to write a function that computes the factorial of an integer greater than or equal to 0. Recall that the factorial of a number n is n * (n-1) * (n-2) * … * 1. The factorial of 0 is defined to be 1. 
#For this Part you will need to write four different versions of the Factorial function:

#1. Factorial_loop: a version that computes the factorial of an integer using looping (such as a for loop)

factorial_loop <- function(x){
	n <- 1
	if(x==0){
		return(1)
	}
	for(i in 1:x) {
		n <-n*((1:x)[i])
		}
	return(n)
}

#checking value
factorial_loop(5)


#2. Factorial_reduce: a version that computes the factorial using the reduce() function in the purrr package. Alternatively, you can use the Reduce() function in the base package.

factorial_reduce <- function(n){
	if(n==0){
		return(1)
	}
	else{
		reduce(c(1:n), function(x, y){
			y*x
		})
	}
}

#checking value
factorial_reduce(5)


#3. Factorial_func: a version that uses recursion to compute the factorial.

factorial_rec <- function(n){
	if(n==0){
		return(1)
	}
	else{
		n * factorial_rec(n-1)
	}
}

#checking value
factorial_rec(5)


#4. Factorial_mem: a version that uses memoization to compute the factorial.

memory <- c(0, 1, rep(NA, 100))

factorial_mem <- function(n){
	if(!is.na(memory[n])){ 
		return(memory[n])
	}
	else{
		fac <- factorial(n)
		memory[n] <<- fac
	}
	memory[5]
}

#checking value
factorial_mem(5)


#After writing your four versions of the Factorial function, use the microbenchmark package to time the operation of these functions and provide a summary of their performance. In addition to timing your functions for specific inputs, make sure to show a range of inputs in order to demonstrate the timing of each function for larger inputs.

for(i in seq(1, 20, by =10)){
	print(map(i, ~ (microbenchmark(factorial_loop(i), factorial_reduce(i), factorial_rec(i), factorial_mem(i)))))
}
