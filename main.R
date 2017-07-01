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

MAIN <- function(data, PDFfilename="ga.pdf", .problem="TSP" ,pop.size=20, pr.elit=0.15, pr.cross=0.6, pr.mut=0.1, op.cross="pmx", op.mut="sim", .time.limit=10)
{

	# Calculate the number of cores
	n_cores <- detectCores() -1 # usually 3

	# Initiate cluster
	cluster	<- makeCluster(n_cores)
	registerDoParallel(cluster)

	if (.problem=="LOP") {
		data <- data - min(data)
		data <- data / max(data)
	}

	if (length(pop.size)>1) 
	{
		legend 	<- "population size"
		array 	<- pop.size 
		results <- (foreach(pops.size=pop.size, .combine=rbind, .export=ls(.GlobalEnv)) %dopar% 
			SOLVE(data, .problem, pops.size, pr.elit, pr.cross , pr.mut, op.cross, op.mut, .time.limit))
	}

	else if (length(pr.elit)>1) 
	{
		legend 	<- "elitism percentage"
		array 	<- pr.elit 
		results <- (foreach(prs.elit=pr.elit, .combine=rbind, .export=ls(.GlobalEnv)) %dopar% 
			SOLVE(data, .problem, pop.size, prs.elit, pr.cross , pr.mut, op.cross, op.mut, .time.limit))	
	}
		
	else if (length(pr.cross)>1) 
	{
		legend 	<- "probability crossover"
		array 	<- pr.cross
		results <- (foreach(prs.cross=pr.cross, .combine=rbind, .export=ls(.GlobalEnv)) %dopar% 
			SOLVE(data, .problem, pop.size, pr.elit, prs.cross , pr.mut, op.cross, op.mut, .time.limit))	
	}

	else if (length(pr.mut)>1) 
	{
		legend 	<- "probability mutation"
		array 	<- pr.mut
		results <- (foreach(prs.mut=pr.mut, .combine=rbind, .export=ls(.GlobalEnv)) %dopar% 
			SOLVE(data, .problem, pop.size, pr.elit, pr.cross , prs.mut, op.cross, op.mut, .time.limit))
	}
		
	else if (length(op.cross)>1)
	{
		legend 	<- "operator crossover"
		array 	<- op.cross
		results <- (foreach(ops.cross=op.cross, .combine=rbind, .export=ls(.GlobalEnv)) %dopar% 
			SOLVE(data, .problem, pop.size, pr.elit, pr.cross , pr.mut, ops.cross, op.mut, .time.limit))
	}
		
	else if (length(op.mut)>1) 
	{
		legend 	<- "operator mutation"
		array 	<- op.mut 
		results <- (foreach(ops.mut=op.mut, .combine=rbind, .export=ls(.GlobalEnv)) %dopar% 
			SOLVE(data, .problem, pop.size, pr.elit, pr.cross , pr.mut, op.cross, ops.mut, .time.limit))
	}
	
	else 
	{	
		array 	<- c("")
		legend 	<- paste("problem", .problem, "pop size", pop.size, "pr.elit", pr.elit, "pr.cross", pr.cross, "pr.mut", pr.mut, "op.cross", op.cross, "op.mut", op.mut)
		results <- SOLVE(data, .problem, pop.size, pr.elit, pr.cross , pr.mut, op.cross, op.mut, .time.limit)
	}

	aux_results 	<- list()
	aux_qu_results 	<- list()
	matrix_start 	<- list()
	matrix_end 		<- list()
	n_its 			<- c()
	quarters 		<- c()
	quarter_results <- c()
	init_results	<- c()
	end_results 	<- c()
	counter			<- 1

	for (its in results[,2])
	{
		n_its 	 <- c(n_its,its)
		quarters <- c(quarters, round(its*0.75))
	}

	if (.problem=="TSP") 
	{
		for (result in results[,1]) 
		{
			total_its    	<- n_its[counter]
			curr_quarter 	<- quarters[counter]
			matrix_start 	<- append(matrix_start,  	list(plot_data(data,result[,1]$best_permutations_row)))
			matrix_end   	<- append(matrix_end,    	list(plot_data(data,result[,total_its]$best_permutations_row)))
			aux_results  	<- append(aux_results,   	list(mapply(function(x) result[,x]$min_results, 1:total_its)))
			aux_qu_results 	<- append(aux_qu_results, 	list(mapply(function(x) result[,x]$min_results, curr_quarter:total_its)))
			init_results 	<- c(init_results,			result[,1]$min_results)
			end_results 	<- c(end_results, 		 	result[,total_its]$min_results)
			quarter_results <- c(quarter_results, 	 	result[,curr_quarter]$min_results)
			counter			<- counter+1
		}
	}

	if (.problem=="LOP" || .problem=="BER")
	{
		for (result in results[,1]) 
		{
			total_its    	<- n_its[counter]
			curr_quarter 	<- quarters[counter]
			matrix_start 	<- append(matrix_start,  	list(get_matrix(result[,1]$best_permutations_row,result[,1]$best_permutations_col, data)))
			matrix_end   	<- append(matrix_end,    	list(get_matrix(result[,total_its]$best_permutations_row,result[,total_its]$best_permutations_col, data)))
			aux_results  	<- append(aux_results,   	list(mapply(function(x) result[,x]$min_results, 1:total_its)))
			aux_qu_results 	<- append(aux_qu_results, 	list(mapply(function(x) result[,x]$min_results, curr_quarter:total_its)))
			init_results 	<- c(init_results,			result[,1]$min_results)
			end_results 	<- c(end_results, 		 	result[,total_its]$min_results)
			quarter_results <- c(quarter_results, 	 	result[,curr_quarter]$min_results)
			counter			<- counter+1
		}
	}

	# plots

	vplayout 		<- function(x, y) viewport(layout.pos.row = x, layout.pos.col = y)

	max_init_result <- max(init_results)
	max_qu_result 	<- max(quarter_results) 
	min_end_result 	<- min(end_results)
	counter 		<- max(1,counter-1)

	pdf(file=PDFfilename)
	for (iterator in 1:counter)
	{
		iterations 		<- n_its[iterator]
		curr_quarter	<- quarters[iterator]
		arr_results 	<- unlist(aux_results[iterator])
		arr_qu_results 	<- unlist(aux_qu_results[iterator])
		main_title 		<- paste(legend, array[iterator])

		plot1 			<- ggplot(data.frame(x=1:iterations, y=arr_results), aes(1:iterations, arr_results)) +  
			labs(x="number of iterations reached",y="results",title=main_title) + 
			ylim(c(min_end_result, max_init_result)) + 
			geom_point(size=2) +
 	 		scale_colour_gradientn(colours = terrain.colors(7))
		
		plot2 			<- qplot(x=curr_quarter:iterations, y=arr_qu_results, 
			xlab="number of iterations reached", ylab="results", 
			main=main_title,
			ylim=c(min_end_result, max_qu_result)) +
			geom_point(size=2) +
			scale_colour_distiller(palette = "RdPu")

		grid.arrange(plot1,plot2)

		if (.problem=="TSP") 
		{
			repr_start_x		<- mapply(function(x) 	matrix_start[iterator][[1]][x,1],1:nrow(data))
			repr_start_y		<- mapply(function(x) 	matrix_start[iterator][[1]][x,2],1:nrow(data))
			repr_end_x			<- mapply(function(x)	matrix_end[iterator][[1]][x,1],1:nrow(data))
			repr_end_y			<- mapply(function(x)	matrix_end[iterator][[1]][x,2],1:nrow(data))

			hamiltonian_graph 	<- qplot(xlab="X-axe", ylab="Y-axe",
			repr_start_x,
			repr_start_y,
			main=paste(main_title,paste("with result",init_results[iterator]))) +
			geom_point(color="palevioletred1") +
			geom_line(color="grey")

			hamiltonian_graph_end <- qplot(xlab="X-axe", ylab="Y-axe",
			repr_end_x,
			repr_end_y,
			main=paste(main_title,paste("with result",end_results[iterator]))) +
			geom_point(color="palevioletred1") +
			geom_line(color="grey")


		
			grid.arrange(hamiltonian_graph)
			grid.newpage()
			grid.arrange(hamiltonian_graph_end)
		}
		
		if (.problem=="LOP")
		{
			levelplot <- myImagePlot(matrix_end[iterator][[1]], xlab="X-axe", ylab="Y-axe", zlim=c(0,1))
			print(levelplot)

		}	

		if (.problem=="BER")
		{
			levelplot <- plot.bertin(matrix_end[iterator][[1]])
			print(levelplot)
		} 	

	}

	dev.off()

	stopCluster(cluster)

}

plot_data <- 
function(data, perm){
	
	x <- c()
	y <- c()
	
	for (i in perm)
	{
		x <- c(x,data[i,1]) # x
		y <- c(y,data[i,2]) # y
	}	
	
	cbind(x,y)
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
 reverse <- nrow(x) : 1
 yLabels <- yLabels[reverse]
 x <- x[reverse,]

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