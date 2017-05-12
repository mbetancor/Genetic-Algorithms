SOLVE_TSP <- function(data, number_offspring) {
	init_population <- replicate(number_offspring,list(sample(nrow(data))))
	
}

## Calculates the upper bound triangle of a matrix

UP_TRIANGLE_SUM <- function(matrix) {
	ans <- 0

	while(length(matrix)>1){

		for (i in matrix[1,]){
			ans <- ans + i
		}

		matrix <- matrix[-1,-1]
	}

	ans + matrix
}


BERTIN_SUM <- function(matrix) {

	ans <- 0
	n <- ncol(matrix)
	reverse_iterator <- n

	while(reverse_iterator>1){

		for (m in 2:reverse_iterator){
			for (i in 1:n){
				ans <- ans + (matrix[1,n] + matrix[m,n]) %% 2
			}
		}

		matrix <- matrix[-1,]
		reverse_iterator <- reverse_iterator -1
	}

	ans
}

GRAPH_TRIANGULATION <- function(matrix) {


}

DRAW_GRAPH <- function(matrix, states) {

}


AUX_GET_DISTANCE <- function(perm,data){

	nrows <- nrow(data)
	len <- length(perm)
	m <- matrix(0,len,2)
	d <- 0
	
	for (i in 1:len){
		num <- perm[i]
		m[i,1] <- data[num,1]
		m[i,2] <- data[num,2]

	}

	for (k in 1:(nrows-1)){
		k2 <- k+1
		x1 <- m[k,1]
		x2 <- m[k2,1]
		y1 <- m[k,2]
		y2 <- m[k2,2]
		d <- d + EUC_DIST(x1,y1,x2,y2)
	}

 	d
}



EUC_DIST <- function(x1,y1,x2,y2){

	x <- (x1-x2)**2
	y <- (y1-y2)**2
	(x+y)**(0.5)
}


REMOVE_DUPLICATES <- function(elements, sorted_list) {

	i <- length(sorted_list)

	while(i>0) {
		new_pos <- which(sorted_list==elements[i])

		if (length(new_pos)) {
			elements[i] <- elements[new_pos]
			elements <- elements[-new_pos]
			sorted_list <- sorted_list[-new_pos]
		}

		else {
			i <- i-1
		}
	}

	print(elements)
	print(sorted_list)

}

