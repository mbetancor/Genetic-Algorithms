source('generations.R')

SOLVE <- function(data, .problem="TSP", pop.size=20, pr.elit=0.15, pr.cross=0.6, pr.mut=0.1, op.cross="pmx", op.mut="sim", .time.limit=60)
{
	t 					<- Sys.time()
	elitists 			<- max(round(pop.size*pr.elit),3)
	init_population_row <- replicate(pop.size,list(sample(nrow(data)))) 
		
	if (.problem=="LOP" || .problem=="BER")
	{
		init_population_col <- replicate(pop.size,list(sample(ncol(data))))
		pop_results 		<- fun_cost(data,init_population_row, init_population_col, .problem)
	}
	
	if (.problem=="TSP")
	{
		print("pasa")
		pop_results 	<- fun_cost(data,init_population_row)
	}

	order_results 			<- order(pop_results)
	min_results 			<- c()
	best_permutations_row 	<- c()
	best_permutations_col 	<- c() 
	n_its 					<- 0
	
	if (.problem=="TSP")
	{
		while(as.numeric(Sys.time()-t, units="secs") < .time.limit) 
		{
			init_population_row 	<- flat(init_population_row,order_results[1:elitists])
			init_population_row 	<- generate_offspring(init_population_row, op.cross, op.mut, pop.size, pr.cross, pr.mut)
			pop_results 			<- fun_cost(data, init_population_row)
			order_results 			<- order(pop_results)
			min_results 			<- c(min_results,min(pop_results))
			best_permutations_row 	<- append(best_permutations_row, init_population_row[order_results[1]])
			n_its 					<- n_its +1
		}

		ANSWER <- list(rbind(best_permutations_row,min_results),n_its)		
	}

	if (.problem=="LOP" || .problem=="BER")
	{
		while(as.numeric(Sys.time()-t, units="secs") < .time.limit) 
		{
			init_population_row 	<- flat(init_population_row,order_results[1:elitists])
			init_population_col 	<- flat(init_population_col,order_results[1:elitists])
			init_population_row 	<- generate_offspring(init_population_row,op.cross, op.mut, pop.size, pr.cross, pr.mut)
			init_population_col 	<- generate_offspring(init_population_col,op.cross, op.mut, pop.size, pr.cross, pr.mut)
			pop_results 			<- fun_cost(data, init_population_row, init_population_col, .problem)
			order_results 			<- order(pop_results)
			min_results 			<- c(min_results,min(pop_results))
			best_permutations_row 	<- append(best_permutations_row, init_population_row[order_results[1]])
			best_permutations_col 	<- append(best_permutations_col, init_population_col[order_results[1]])
			n_its 					<- n_its +1
		}

		ANSWER <- list(rbind(rbind(best_permutations_row, best_permutations_col),min_results),n_its)	
	}

	return(ANSWER)
		
}	

fun_cost <- 
function(data, init_population_row, init_population_col = 0, .problem="TSP")
{
	if (.problem=="TSP")
		COST <- mapply(function(x) fun_cost.tsp(x,data), init_population_row)
	if (.problem=="LOP") 
		COST <- mapply(function(x,y) fun_cost.lop(x,y,data), init_population_row, init_population_col)
	if (.problem=="BER") 
		COST <- mapply(function(x,y) fun_cost.bertin(x, y, data), init_population_row, init_population_col)

	return(COST)
}

fun_cost.tsp <- 
function(perm,data)
{
	nrows 	<- nrow(data)
	len 	<- length(perm)
	m 		<- matrix(0,len,2)
	d 		<- 0
	
	for (i in 1:len)
	{
		num 	<- perm[i]
		m[i,1] 	<- data[num,1]
		m[i,2] 	<- data[num,2]
	}

	for (k in 1:(nrows-1))
	{
		k2 	<- k+1
		x1 	<- m[k,1]
		x2 	<- m[k2,1]
		y1 	<- m[k,2]
		y2 	<- m[k2,2]
		d 	<- d + euc_dist(x1,y1,x2,y2)
	}

 	d
}

euc_dist <- 
function(x1,y1,x2,y2)
{
	x <- (x1-x2)**2
	y <- (y1-y2)**2
	(x+y)**(0.5)
}


fun_cost.lop <- 
function(list_row, list_col, data) 
{	
	len 		<- length(list_row)
	new_matrix 	<- get_matrix(list_row, list_col ,data)
	ans 		<- new_matrix[len,len] # elemento del final

	while(len>1) 
	{
		for (j in 1:len)
		{
			ans <- ans + new_matrix[j,1]
		}
		
		new_matrix 	<- new_matrix[-1,-1]
		len 		<- len -1	
	}

	ans
}

fun_cost.bertin <- 
function(list_row, list_col, data)
{
	ans 		<- 0
	ncols 		<- length(list_col)
	nrows 		<- length(list_row)
	new_matrix 	<- get_matrix(list_row, list_col ,data)

	for (i in 1:(ncols-1)){
		c1 <- new_matrix[i,]
		c2 <- new_matrix[i+1,]
		ans <- sum((c1+c2) %%2) + ans
	}

	ans
}


get_matrix <- 
function(list_row, list_col, data) 
{
	ncols 		<- length(list_col)
	nrows 		<- length(list_row)
	new_matrix 	<- matrix(c(integer(nrows*ncols)), nrow=nrows, ncol=ncols)
	
	for (i in 1:nrows)
	{
		pos 			<- list_row[i]
		new_matrix[i,] 	<- unlist(data[pos,])
	}

	aux_matrix 	<- new_matrix

	for (i in 1:ncols)
	{
		pos 			<- list_col[i]
		new_matrix[,i] 	<- unlist(aux_matrix[,pos])
	}

	new_matrix

}

flat <- 
function(vector, order)
{
	ans <- c()
	for (i in order){
		ans <- c(ans,vector[i])
	}
	ans
}