SOLVE <- function(data, problem="tsp", pop.size=20, pr.elit=0.15, pr.cross=0.6, pr.mut=0.1, op.cross="pmx", op.mut="sim", .time.limit=60)
{
	t <- Sys.time()
	elitists <- max(round(pop.size*pr.elit),3)
	init_population <- replicate(pop.size,list(sample(nrow(data)))) 
	pop_results <- fun_cost(problem,init_population,data)
	order_results <- order(pop_results)
	min_results <- c()
	best_permutations <- c() 
	n_its <- 0
	
	while(as.numeric(Sys.time()-t, units="secs") < .time.limit) {
		init_population <- flat(init_population,order_results[1:elitists])
		init_population <- generate_offspring(init_population, op.cross , op.mut , pop.size, pr.cross, pr.mut)
		pop_results <- fun_cost(problem,init_population,data)
		order_results <- order(pop_results)
		min_results <- c(min_results,min(pop_results))
		best_permutations <- append(best_permutations, init_population[order_results[1]])
		n_its <- n_its +1
	}
	
	list(rbind(best_permutations,min_results),n_its)
}	



fun_cost <- 
function(.problem="tsp", init_population, data)
{
	if (.problem=="tsp") 
		COST <- mapply(function(x) fun_cost.tsp(x,data), init_population)
	if (.problem=="lop") 
		COST <- mapply(function(x) f_cost.lop(init_population[,x],data), 1:length(init_population))
	if (.problem=="bertin") #todo
		COST <- fun_cost.bertin(init_population, data)

	return(COST)
}

fun_cost.tsp <- 
function(perm,data)
{
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
		d <- d + euc_dist(x1,y1,x2,y2)
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
function(matrix, data) 
{	
	m1 <- unlist(matrix[1])
	# rearrange original matrix
	len <- length(m1)

	new_matrix <- get_matrix(matrix,data)

	ans <- new_matrix[len,len] # elemento del final

	while(len>1) {
		for (j in 1:len){
			ans <- ans + new_matrix[j,1]
		}
		new_matrix <- new_matrix[-1,-1]
		len <- len -1	
	}

	ans
}


get_matrix <- 
function(matrix,data) 
{
	m1 <- unlist(matrix[1])
	m2 <- unlist(matrix[2])
	# rearrange original matrix
	new_matrix <- matrix(c(integer(len*len)), nrow=len, ncol=len)
	
	for (i in 1:len){
		pos <- m1[i]
		new_matrix[i,] <- unlist(data[pos,])
	}

	aux_matrix <- new_matrix

	for (i in 1:len){
		pos <- m2[i]
		new_matrix[,i] <- unlist(aux_matrix[,pos])
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