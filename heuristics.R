source('problems.R')

heuristic.tsp <- 
function(data, stress) {
	
	len 			<- nrow(data)
	adn_sequences 	<- c()
	order_x 		<- order(data[,1]) 
	order_y 		<- order(data[,2])
	cut_points_x 	<- sort(sample(1:len,stress,replace=F))
	cut_points_y 	<- sort(sample(1:len,stress,replace=F))

	for (i in 1:(stress-1)) 
	{
		px_start 	<- cut_points_x[i]
		px_end 		<- cut_points_x[i+1]
		sequence_x	<- order_x[px_start:px_end]

		for(j in 1:(stress-1))
		{
			py_start 		<- cut_points_y[j]
			py_end 			<- cut_points_y[j+1]
			sequence_y		<- order_y[py_start:py_end]
			adn_intersect 	<- intersect(sequence_x, sequence_y) 
			adn_sequences 	<- c(adn_sequences, adn_intersect)
		}
	}

	adn_sequences 	<- unique(adn_sequences)
	pool_remaining 	<- sample(setdiff(c(1:len), adn_sequences))

	return(list(c(adn_sequences,pool_remaining)))
}

heuristic.lop <- 
function(data,stress) 
{
	mapply ( function(x) list(heuristic.lop_deep(data,x)), 1:stress)
}

heuristic.lop_deep <- 
function(data,n) {

	len 		<- nrow(data)
	perm_row 	<- rev(order(mapply(function(x) sum(data[x,][n:len]), 	1:len)))
	aux_matrix 	<- get_matrix(perm_row,c(1:len), data)
	perm_col 	<- order(mapply(function(x) sum(aux_matrix[,x][n:len]), 1:len))

	return(cbind(perm_row,perm_col)) 
}

# heuristic.tab <- 
# function(data, perc_row = 1, perc_col = 1)
# {
# 	nrows 		<- nrow(data)
# 	ncols 		<- ncol(data)

# 	sample_row	<- sample(1:nrows, round(nrows*perc_row), replace=F)
# 	sample_col	<- sample(1:ncols, round(ncols*perc_col), replace=F)

# 	perm_row 	<- order(mapply(function(x) sum(data[x,]), sample_row))
# 	aux_matrix 	<- get_matrix(perm_row,c(1:ncols), data)
# 	perm_col 	<- order(mapply(function(x) sum(aux_matrix[,x]),sample_col))

# 	return(list(perm_row,perm_col)) 

# }