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

PMX <- function(str1, str2) {
	len <- length(str1)
	random_cut_points <- sample(1:len,2, replace=F)
	rc1 <- min(random_cut_points)
	rc2 <- max(random_cut_points)

	ms1 <- str1[rc1:rc2]
	ms2 <- str2[rc1:rc2]

	off1[rc1:rc2] <- ms2
	off2[rc1:rc2] <- ms1

	order <- order(ms1)

	index <- 1

	for i in order:
		map_section_1[index] <- ms1[i]
		map_section_2[index] <- ms2[i]
		index <- index+1

	map_sections_arr <- REMOVE_DUPLICATES(map_section_1, map_section_2)

	ms1 <- map_sections_arr[1]
	ms2 <- map_sections_arr[2]

	for i in (1:len)[-(rc1:rc2)]:
		item_list_1 <- which(ms1==str1[i])
		item_list_2 <- which(ms2==str2[i])

		if (length(item_list_1)){
			off1[i] <- ms2[item_list_1]
		}

		else {
			off1[i] <- str1[i]
		}

		if (length(item_list_2)){
			off2[i] <- ms1[item_list_2]
		}

		else {
			off2[i] <- str2[i]
		}

	print(off1)
	print(off2)

}


# Cycle CrossOver (CX)

CX <- function(str1, str2) {
	len <- length(str1)
	off <-  integer(len)
	seeked <- str2[1]
	elem <- str1[1]
	off[1] <- elem
	p <- 0

	while(length(which[off==0])>0){

		i <- which[off==0][1]

		if (seeked==0){
			if (p==0){
				seeked <-str2[i]
				elem <- str1[i]
				off[i] <- elem
			}
			else {
				seeked <-str1[i]
				elem <- str2[i]
				off[i] <- elem
			}
		}
		if (p==0){ #parent2

			while(str2[i]!=elem){
				i <- i+1
			}

			elem <- str1[i]
		} 
		else {

			while(str1[i]!=elem){
				i <- i+1
			}

			elem <- str2[i]

		}

		off[i] <- elem

		if (elem == seeked) {
			 p <- (p+1) %% 2
			 seeked <- 0

		}
	
	}

	off
}
