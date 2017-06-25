MUTATION <- 
function(str,.op.mutation ="sim", .pr.mut = 0.1) 
{
	len 	<- length(str)
	size 	<-  max(2,round(len * .pr.mut))
	
	if(.op.mutation == "ivm") 
		MUT <- (mutation.ivm(str, size, len))
	if(.op.mutation == "ism") 
		MUT <- (mutation.ism(str, size, len))
	if(.op.mutation == "em") 
		MUT <- (mutation.em(str,size, len))
	if(.op.mutation == "sim") 
		MUT <- (mutation.sim(str, size, len))

	return(MUT)
}

mutation.ivm <- 
function(str,size, len) 
{
	index1  <- sample(1:(len-size),1)
	index2  <- sample((index1+1):len,1)
	subtour <- str[index2:index1]
	str 	<- str[!str %in% subtour]
	len 	<- length(str)
	ran 	<- sample(1:len,1)
	append(str[0:(ran-1)],append(subtour,str[ran:len]))
}

mutation.ism <- 
function(str,size, len) 
{
	for (i in 1:size) {
		random_points 	<- sample(2:(len-1), 2, replace=F)
		index1 			<- random_points[1]
		index2 			<- random_points[2]
		number 			<- str[index1]
		str 			<- str[-index1]
		str 			<- append(append(str[1:(index2-1)],number),str[index2:(len-1)])
	}

	str
}

mutation.em <- 
function(str,size, len) 
{
	for (i in 1:size) {
		random_points 	<- sample(1:len, 2, replace=F)
		index1 			<- random_points[1]
		index2			<- random_points[2]
		number1 		<- str[index1]
		number2 		<- str[index2]
		str[index1]		<- number2
		str[index2] 	<- number1
	}

	str
}

mutation.sim <- 
function(str,size, len)
{
	index1 				<- sample(1:(len-size),1)
	index2 				<- sample((index1+1):len,1)
	str[index1:index2] 	<- str[index2:index1]
	str
}