library(ggplot2)
library(grid)
library(gridExtra)
library(foreach)
library(parallel)
library(doParallel)
library(colorspace)
library(RColorBrewer)
library(bertin)

source('problems.R')
source('heuristics.R')

ISLANDS <-
	function(data, PDFfilename="ga.pdf", .problem="TSP",
		.pop.size=20, .pr.elit=0.15, pr.cross=0.6, pr.mut=0.1, op.cross="pmx", op.mut="sim",
		init.individuals = 200, .per.sampling=0.1,
		.emigration.policy = TRUE, .per.limit=0.05,
		.number.its = 40, .time.it= 15,
		.enable.graphics = TRUE )
{

	vplayout 	<- function(x, y) viewport(layout.pos.row = x, layout.pos.col = y)
	pdf(file=PDFfilename)
	if (.problem=="TSP")
	{
		R_PERM	<- sample(1:nrow(data)) # RANDOM_PERMUTATION
		INIT 	<- fun_cost.tsp(R_PERM, data)
		ANS 	<- ISLANDS.TSP(data, PDFfilename, .pop.size, .pr.elit, pr.cross, pr.mut, op.cross, op.mut,
			init.individuals, .per.sampling,
			.emigration.policy, .per.limit,
			.number.its, .time.it
			)
		ANS_PERM<- ANS[1]$best_permutations
		ANS_RES	<- ANS[2]$min_results

		M_START <- mapply( function(x) cbind(data[x,1],data[x,2]), R_PERM)
		M_END	<- mapply( function(x) cbind(data[x,1],data[x,2]), ANS_PERM)

		hamiltonian_graph 	<- qplot(xlab="X-axe", ylab="Y-axe",
			M_START[1,],
			M_START[2,],
			main=paste("with result", INIT) ) +
			geom_point(color="palevioletred1") +
			geom_line(color="grey")

		hamiltonian_graph_2 <- qplot(xlab="X-axe", ylab="Y-axe",
			M_END[1,],
			M_END[2,],
			main=paste("with result", ANS_RES) ) +
			geom_point(color="palevioletred1") +
			geom_line(color="grey")

		grid.arrange(hamiltonian_graph)
		grid.newpage()
		grid.arrange(hamiltonian_graph_2)
	}

	if (.problem=="LOP")
	{
		data <- data - min(data)
		data <- data / max(data)

		R_PERM1	<- sample(1:nrow(data)) # RANDOM_PERMUTATION
		R_PERM2 <- sample(1:ncol(data))
		INIT 	<- fun_cost.lop(R_PERM1, R_PERM2, data)
		ANS 	<- ISLANDS.LOP(data, PDFfilename, .pop.size, .pr.elit, pr.cross, pr.mut, op.cross, op.mut,
			init.individuals, .per.sampling,
			.emigration.policy, .per.limit,
			.number.its, .time.it
			)
		ANS_ROW	<- ANS[1]$best_permutations_row
		ANS_COL <- ANS[2]$best_permutations_col
		ANS_RES	<- ANS[3]$min_results

		M_START	<- get_matrix(R_PERM1, R_PERM2, data)
		M_END 	<- get_matrix(ANS_ROW, ANS_COL, data)

		PLOT_START 	<- myImagePlot(M_START, xlab="X-axe", ylab="Y-axe", zlim=c(0,1))
		PLOT_END 	<- myImagePlot(M_END, xlab="X-axe", ylab="Y-axe", zlim=c(0,1))

		print(PLOT_START)
		print(PLOT_END)
	}

	if (.problem=="TAB")
	{

		R_PERM1	<- sample(1:nrow(data)) # RANDOM_PERMUTATION
		R_PERM2 <- sample(1:ncol(data))
		INIT 	<- fun_cost.tab(R_PERM1, R_PERM2, data)
		ANS 	<- ISLANDS.TAB(data, PDFfilename, .pop.size, .pr.elit, pr.cross, pr.mut, op.cross, op.mut,
			init.individuals, .per.sampling,
			.emigration.policy, .per.limit,
			.number.its, .time.it
			)
		ANS_ROW	<- ANS[1]$best_permutations_row
		ANS_COL <- ANS[2]$best_permutations_col
		ANS_RES	<- ANS[3]$min_results

		M_START	<- get_matrix(R_PERM1, R_PERM2, data)
		M_END 	<- get_matrix(ANS_ROW, ANS_COL, data)

		PLOT_START 	<- myImagePlot(M_START, xlab="X-axe", ylab="Y-axe", zlim=c(0,1))
		PLOT_END 	<- myImagePlot(M_END, xlab="X-axe", ylab="Y-axe", zlim=c(0,1))

		print(PLOT_START)
		print(PLOT_END)
	}

	dev.off()

	return(c(INIT, ANS_RES))

}

ISLANDS.TSP <- 
	function(data, PDFfilename="ga.pdf", pop.size=20, pr.elit=0.15, pr.cross=0.6, pr.mut=0.1, op.cross="pmx", op.mut="sim",
		init.individuals = 200, .per.sampling = 0.1,
		.emigration.policy = TRUE, .per.limit=0.05,
		.number.its=10 ,.time.it=15)
{

	# Calculate the number of cores
	n_cores <- detectCores() -1 # usually 3

	init_population <- replicate(n_cores, replicate(200, heuristic.tsp(data,40)))
	len_population 	<- max(round(pr.elit*pop.size),3)

	if (.emigration.policy)
	{
		for (i in 1:.number.its)
		{
			# Initiate cluster
			cluster	<- makeCluster(n_cores)
			registerDoParallel(cluster)

			results <- (foreach(init_populations=init_population, .combine=cbind, .export=ls(.GlobalEnv)) %dopar% 
					SOLVE.TSP(data, init_population, op.cross, op.mut, pr.cross, pr.mut, pop.size, pr.elit, .time.it))

			migrations <- sample(1:n_cores)

			init_population <- 
			mapply ( 
				function(x,y) 
					list( c( results[1,][ (1:(len_population-1)) + x*len_population],
					results[1,][len_population + y*len_population])),
					0:(n_cores-1),
					migrations-1
			)
		}

		stopCluster(cluster)
		closeAllConnections()
	}


	else {

		for (i in 1:.number.its)
		{
			# Initiate cluster
			cluster	<- makeCluster(n_cores)
			registerDoParallel(cluster)

			results <- (foreach(init_populations=init_population, .combine=cbind, .export=ls(.GlobalEnv)) %dopar% 
					SOLVE.TSP(data, init_population, op.cross, op.mut, pr.cross, pr.mut, pop.size, pr.elit, .time.it))

			nums	<- results[2,]
			per 	<- nums / min(nums)
			worse	<- which (per >= (1+ .limit.perf))
			best 	<- order(per)[1:len_population]
			it 		<- 1

			for (w in worse)
			{
				fittest <- best[it]
				results[1,w]$best_permutations <- results[1,fittest]$best_permutations
				it 		<- ( it ) %% (len_population-1) + 1 
			}
			
			init_population <- mapply ( function(x)  list(results[1,][(1:len_population)*x]), 1:n_cores )

			stopCluster(cluster)
			closeAllConnections()
		}
	}
	
	flat_res	<- unlist(results[2,])
	ANS_INDEX 	<- which(flat_res==min(flat_res))

	return(results[,ANS_INDEX])
}

SOLVE.TSP <- 
	function(data, init_population, op.cross="pmx", op.mut="sim", pr.cross=0.6, pr.mut=0.1, pop.size=20, pr.elit=0.15, .time.it=10)
{
	t 					<- Sys.time()
	elitists 			<- max(round(pop.size*pr.elit),3)
	pop_results 		<- fun_cost(data,init_population)
	order_results		<- order(pop_results)

	while(as.numeric(Sys.time()-t, units="secs") < .time.it) 
	{
		init_population 		<- flat(init_population,order_results[1:elitists])
		init_population			<- generate_offspring(init_population, op.cross, op.mut, pop.size, pr.cross, pr.mut)
		pop_results 			<- fun_cost(data, init_population)
		order_results 			<- order(pop_results)
		min_results 			<- pop_results[order_results[1:elitists]]
		best_permutations 	 	<- init_population[order_results[1:elitists]]
	}

	ANSWER <- rbind(best_permutations,min_results)	

	return(ANSWER)
}

ISLANDS.LOP <- 
	function(data, PDFfilename="ga.pdf", pop.size=20, pr.elit=0.15, pr.cross=0.6, pr.mut=0.1, op.cross="pmx", op.mut="sim",
		init.individuals = 40, .per.sampling = 0.1,
		.emigration.policy = TRUE, .per.limit=0.05,
		.number.its=10 ,.time.it=15)
{
	# Calculate the number of cores
	n_cores <- detectCores() -1 # usually 3

	init_population <- heuristic.lop(data, init.individuals)
	init_population_row <- mapply(function(x) list(init_population[[x]][1,]), 1:init.individuals)
	init_population_row <- matrix(rep(init_population_row,n_cores), nrow=init.individuals, ncol=n_cores)
	init_population_col <- mapply(function(x) list(init_population[[x]][2,]), 1:init.individuals)
	init_population_col <- matrix(rep(init_population_col,n_cores), nrow=init.individuals, ncol=n_cores)
	len_population 	<- max(round(pr.elit*pop.size),3)


	if (.emigration.policy)
	{
		for (i in 1:.number.its)
		{
			# Initiate cluster
			cluster	<- makeCluster(n_cores)
			registerDoParallel(cluster)

			results <- (foreach(x=1:n_cores, .combine='cbind', .export=ls(.GlobalEnv)) %dopar% 
				SOLVE.LOP(data, init_population_row[,x], init_population_col[,x], pop.size, pr.elit, pr.cross, pr.mut, op.cross, op.mut, .time.it))

			migrations <- sample(1:n_cores)

			init_population_row <- 
			mapply ( 
				function(x,y) 
					cbind( c( results[1,][ (1:(len_population-1)) + x*len_population],
					results[1,][len_population + y*len_population])),
					0:(n_cores-1),
					migrations-1
			)

			init_population_col <- 
			mapply ( 
				function(x,y) 
					cbind( c( results[2,][ (1:(len_population-1)) + x*len_population],
					results[2,][len_population + y*len_population])),
					0:(n_cores-1),
					migrations-1
			)

		}

		stopCluster(cluster)
		closeAllConnections()
	}


	else {

		for (i in 1:.number.its)
		{
			# Initiate cluster
			cluster	<- makeCluster(n_cores)
			registerDoParallel(cluster)

			results <- (foreach(x=1:n_cores, .combine='c', .export=ls(.GlobalEnv)) %dopar% 
					SOLVE.LOP(data, init_population_row[,x], init_population_col[,x], pop.size, pr.elit, pr.cross, pr.mut, op.cross, op.mut, .time.it)
					)

			nums	<- results[3,]
			per 	<- nums / min(nums)
			worse	<- which (per >= (1+ .limit.perf))
			best 	<- order(per)[1:len_population]
			it 		<- 1

			for (w in worse)
			{
				fittest <- best[it]
				results[1,w]$best_permutations <- results[1,fittest]$best_permutations
				it 		<- ( it ) %% (len_population-1) + 1 
				results[2,w]$best_permutations <- results[2,fittest]$best_permutations
				it 		<- ( it ) %% (len_population-1) + 1

			}
			
			init_population <- mapply ( function(x)  list(results[1,][(1:len_population)*x]), 1:n_cores )

			stopCluster(cluster)
			closeAllConnections()
		}
	}
	
	flat_res	<- unlist(results[3,])
	ANS_INDEX 	<- which(flat_res==min(flat_res))

	return(results[,ANS_INDEX])
}

SOLVE.LOP <- 
function(data, init_population_row, init_population_col, pop.size=20, pr.elit=0.15, pr.cross=0.6, pr.mut=0.1, op.cross="pmx", op.mut="sim", .time.it=60)
{
	t 						<- Sys.time()
	elitists 				<- max(round(pop.size*pr.elit),3)
	pop_results 			<- fun_cost(data,init_population_row, init_population_col, .problem="LOP")
	order_results 			<- order(pop_results)
	min_results 			<- c(0)
	best_permutations_row 	<- c()
	best_permutations_col 	<- c() 

	while(as.numeric(Sys.time()-t, units="secs") < .time.it) 
	{
		init_population_row 	<- flat(init_population_row,order_results[1:elitists])
		init_population_col 	<- flat(init_population_col,order_results[1:elitists])
		init_population_row 	<- generate_offspring(init_population_row,op.cross, op.mut, pop.size, pr.cross, pr.mut)
		init_population_col 	<- generate_offspring(init_population_col,op.cross, op.mut, pop.size, pr.cross, pr.mut)
		pop_results 			<- fun_cost(data, init_population_row, init_population_col, .problem="LOP")
		order_results 			<- order(pop_results)
		min_results 			<- pop_results[order_results[1:elitists]]
		best_permutations_row 	<- init_population_row[order_results[1:elitists]]
		best_permutations_col 	<- init_population_col[order_results[1:elitists]]
	}

	ANSWER <- rbind(best_permutations_row,best_permutations_col,min_results)	

	return(ANSWER)
		
}


ISLANDS.TAB <- 
	function(data, PDFfilename="ga.pdf", pop.size=20, pr.elit=0.15, pr.cross=0.6, pr.mut=0.1, op.cross="pmx", op.mut="sim",
		init.individuals = 20, .per.sampling = 0.1,
		.emigration.policy = TRUE, .per.limit=0.05,
		.number.its=10 ,.time.it=15)
{
	# Calculate the number of cores
	n_cores <- detectCores() -1 # usually 3

	init_population_row <- replicate(n_cores,replicate(init.individuals,list(sample(nrow(data)))))
	init_population_col <- replicate(n_cores,replicate(init.individuals,list(sample(ncol(data)))))
	len_population 	<- max(round(pr.elit*pop.size),3)

	if (.emigration.policy)
	{
		for (i in 1:.number.its)
		{
			# Initiate cluster
			cluster	<- makeCluster(n_cores)
			registerDoParallel(cluster)

			results <- (foreach(x=1:n_cores, .combine='cbind', .export=ls(.GlobalEnv)) %dopar% 
				SOLVE.TAB(data, init_population_row[,x], init_population_col[,x], pop.size, pr.elit, pr.cross, pr.mut, op.cross, op.mut, .time.it))

			migrations <- sample(1:n_cores)

			init_population_row <- 
			mapply ( 
				function(x,y) 
					cbind( c( results[1,][ (1:(len_population-1)) + x*len_population],
					results[1,][len_population + y*len_population])),
					0:(n_cores-1),
					migrations-1
			)

			init_population_col <- 
			mapply ( 
				function(x,y) 
					cbind( c( results[2,][ (1:(len_population-1)) + x*len_population],
					results[2,][len_population + y*len_population])),
					0:(n_cores-1),
					migrations-1
			)

		}

		stopCluster(cluster)
		closeAllConnections()
	}


	else {

		for (i in 1:.number.its)
		{
			# Initiate cluster
			cluster	<- makeCluster(n_cores)
			registerDoParallel(cluster)

			results <- (foreach(x=1:n_cores, .combine='cbind', .export=ls(.GlobalEnv)) %dopar% 
					SOLVE.TAB(data, init_population_row[,x], init_population_col[,x], pop.size, pr.elit, pr.cross, pr.mut, op.cross, op.mut, .time.it)
					)

			nums	<- results[3,]
			per 	<- nums / min(nums)
			worse	<- which (per >= (1+ .limit.perf))
			best 	<- order(per)[1:len_population]
			it 		<- 1

			for (w in worse)
			{
				fittest <- best[it]
				results[1,w]$best_permutations <- results[1,fittest]$best_permutations
				it 		<- ( it ) %% (len_population-1) + 1 
				results[2,w]$best_permutations <- results[2,fittest]$best_permutations
				it 		<- ( it ) %% (len_population-1) + 1

			}
			
			init_population <- mapply ( function(x)  list(results[1,][(1:len_population)*x]), 1:n_cores )

			stopCluster(cluster)
			closeAllConnections()
		}
	}
	
	flat_res	<- unlist(results[3,])
	ANS_INDEX 	<- which(flat_res==min(flat_res))

	return(results[,ANS_INDEX])

}

SOLVE.TAB <- 
function(data, init_population_row, init_population_col, pop.size=20, pr.elit=0.15, pr.cross=0.6, pr.mut=0.1, op.cross="pmx", op.mut="sim", .time.it=60)
{
	t 						<- Sys.time()
	elitists 				<- max(round(pop.size*pr.elit),3)
	pop_results 			<- fun_cost(data,init_population_row, init_population_col, .problem="TAB")
	order_results 			<- order(pop_results)
	min_results 			<- c(0)
	best_permutations_row 	<- c()
	best_permutations_col 	<- c() 

	while(as.numeric(Sys.time()-t, units="secs") < .time.it) 
	{
		init_population_row 	<- flat(init_population_row,order_results[1:elitists])
		init_population_col 	<- flat(init_population_col,order_results[1:elitists])
		init_population_row 	<- generate_offspring(init_population_row,op.cross, op.mut, pop.size, pr.cross, pr.mut)
		init_population_col 	<- generate_offspring(init_population_col,op.cross, op.mut, pop.size, pr.cross, pr.mut)
		pop_results 			<- fun_cost(data, init_population_row, init_population_col, .problem="TAB")
		order_results 			<- order(pop_results)
		min_results 			<- pop_results[order_results[1:elitists]]
		best_permutations_row 	<- init_population_row[order_results[1:elitists]]
		best_permutations_col 	<- init_population_col[order_results[1:elitists]]
	}

	ANSWER <- rbind(best_permutations_row,best_permutations_col,min_results)	
	return(ANSWER)
}



# source: http://phaget4.org/R/image_matrix.html
# ----- Define a function for plotting a matrix ----- #
myImagePlot <- function(x, ...){
     min <- min(x)
     max <- max(x)
     yLabels <- rownames(x)
     xLabels <- colnames(x)
     title <-c()
  # check for additional function arguments
  if( length(list(...)) ){
    Lst <- list(...)
    if( !is.null(Lst$zlim) ){
       min <- Lst$zlim[1]
       max <- Lst$zlim[2]
    }
    if( !is.null(Lst$yLabels) ){
       yLabels <- c(Lst$yLabels)
    }
    if( !is.null(Lst$xLabels) ){
       xLabels <- c(Lst$xLabels)
    }
    if( !is.null(Lst$title) ){
       title <- Lst$title
    }
  }
# check for null values
if( is.null(xLabels) ){
   xLabels <- c(1:ncol(x))
}
if( is.null(yLabels) ){
   yLabels <- c(1:nrow(x))
}

layout(matrix(data=c(1,2), nrow=1, ncol=2), widths=c(4,1), heights=c(1,1))

 # Red and green range from 0 to 1 while Blue ranges from 1 to 0
 ColorRamp <- rgb( seq(1,1,length=256),  # Red
                   seq(1,0,length=256),  # Green
                   seq(1,1,length=256))  # Blue
 ColorLevels <- seq(min, max, length=length(ColorRamp))

 # Reverse Y axis
 # reverse <- nrow(x) : 1
 # yLabels <- yLabels[reverse]
 # x <- x[reverse,]

 # Data Map
 par(mar = c(3,5,2.5,2))
 image(1:length(xLabels), 1:length(yLabels), t(x), col=ColorRamp, xlab="",
 ylab="", axes=FALSE, zlim=c(min,max))
 if( !is.null(title) ){
    title(main=title)
 }
axis(BELOW<-1, at=1:length(xLabels), labels=xLabels, cex.axis=0.7)
 axis(LEFT <-2, at=1:length(yLabels), labels=yLabels, las= HORIZONTAL<-1,
 cex.axis=0.7)

 # Color Scale
 par(mar = c(3.5,3,2.5,2))
 image(1, ColorLevels,
      matrix(data=ColorLevels, ncol=length(ColorLevels),nrow=1),
      col=ColorRamp,
      xlab="",ylab="",
      xaxt="n")

 layout(1)
}
# ----- END plot function ----- #