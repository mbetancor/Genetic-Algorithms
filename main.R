library(ggplot2)
library(grid)
library(gridExtra)
library(foreach)
library(doParallel)
library(bertin)

source('generations.R')

SOLVE_PARALLEL <- function(data,f_cost, pop_size, p_elitism, p_crossover, p_mutation, fun_crossover, fun_mutation, time_limit){

	# Calculate the number of cores

	no_cores <- detectCores() -1 # usually 3


	# Initiate cluster

	cl <- makeCluster(no_cores)

	registerDoParallel(cl)

	results <- (foreach(ps_crossover=c(0.7,0.8,0.9)) %do% SOLVE(data,f_cost, pop_size, p_elitism, ps_crossover, p_mutation, fun_crossover, fun_mutation, time_limit))	

	stopCluster(cl)

	res1 <- results[[1]][[1]]
	res2 <- results[[2]][[1]]
	res3 <- results[[3]][[1]]

	n_its_1 <- results[[1]][[2]]
	n_its_2 <- results[[2]][[2]]
	n_its_3 <- results[[3]][[2]]

	matrix1_start <- unlist(AUX_PLOT_DATA(data,res1[1,1]))
	matrix2_start <- unlist(AUX_PLOT_DATA(data,res2[1,1]))
	matrix3_start <- unlist(AUX_PLOT_DATA(data,res3[1,1]))

	matrix1_end <- unlist(AUX_PLOT_DATA(data,res1[1,n_its_1]))
	matrix2_end <- unlist(AUX_PLOT_DATA(data,res2[1,n_its_2]))
	matrix3_end <- unlist(AUX_PLOT_DATA(data,res3[1,n_its_3]))


	vplayout <- function(x, y) viewport(layout.pos.row = x, layout.pos.col = y)

	plot1 <- qplot(xlab="number of iterations reached", ylab="results", 1:n_its_1,unlist(res1[2,]), main="0.7 crossover")
	plot2 <- qplot(xlab="number of iterations reached", ylab="results", 1:n_its_2,unlist(res2[2,]), main="0.8 crossover")
	plot3 <- qplot(xlab="number of iterations reached", ylab="results", 1:n_its_3,unlist(res3[2,]), main="0.9 crossover")

	q1 <- round(n_its_1*0.75)
	q2 <- round(n_its_2*0.75)
	q3 <- round(n_its_3*0.75)

	plot10 <- qplot(xlab="number of iterations reached", ylab="results", q1:n_its_1,unlist(res1[2,][q1:n_its_1]), main="0.7 crossover")
	plot11 <- qplot(xlab="number of iterations reached", ylab="results", q2:n_its_2,unlist(res2[2,][q2:n_its_2]), main="0.8 crossover")
	plot12 <- qplot(xlab="number of iterations reached", ylab="results", q3:n_its_3,unlist(res3[2,][q3:n_its_3]), main="0.9 crossover")

	plot4 <- qplot(xlab="X-axe", ylab="Y-axe",matrix1_start[1,], matrix1_start[2,], main="TSP[start] 0.7 crossover") + geom_line(color="steelblue")
	plot5 <- qplot(xlab="X-axe", ylab="Y-axe",matrix2_start[1,], matrix2_start[2,], main="TSP[start] 0.8 crossover") + geom_line(color="steelblue")
	plot6 <- qplot(xlab="X-axe", ylab="Y-axe",matrix3_start[1,], matrix3_start[2,], main="TSP[start] 0.9 crossover") + geom_line(color="steelblue")

	plot7 <- qplot(xlab="X-axe", ylab="Y-axe",matrix1_end[1,], matrix1_end[2,], main="TSP[end] 0.7 crossover") + geom_line(color="darkblue")
	plot8 <- qplot(xlab="X-axe", ylab="Y-axe",matrix2_end[1,], matrix2_end[2,], main="TSP[end] 0.8 crossover") + geom_line(color="darkblue")
	plot9 <- qplot(xlab="X-axe", ylab="Y-axe",matrix3_end[1,], matrix3_end[2,], main="TSP[end] 0.9 crossover") + geom_line(color="darkblue")


	# 9 figures arranged in 3 rows and 3 columns

	grid.arrange(plot1,plot2, plot3, plot10, plot11, plot12, ncol=3, nrow=2)
	#qplot(1:n_its,results[[chosen]])
	grid.newpage()
	grid.arrange(plot4, plot7, plot5, plot8, plot6, plot9, ncol=2, nrow=3)

	print(res1[,n_its_1])
	print(res2[,n_its_2])
	print(res3[,n_its_3])
}



MAIN <- function(data, PDFfilename,problem="tsp" ,pop.size=20, pr.elit=0.15, pr.cross=c(0.4,0.6,0.8), pr.mut=0.1, op.cross="pmx", op.mut="sim", .time.limit=10)
{

	# Calculate the number of cores

	no_cores <- detectCores() -1 # usually 3

	# Initiate cluster

	cl <- makeCluster(no_cores)

	registerDoParallel(cl)

	popsize_len <- length(pop.size)
	prelit_len <- length(pr.elit)
	prcross_len <- length(pr.cross)
	prmut_len <- length(pr.mut)
	opcross_len <- length(op.cross)
	opmut_len <- length(op.mut)

	long <- max(c(popsize_len,prelit_len,prcross_len,prmut_len,opcross_len,opmut_len))
	print(long)

	if (popsize_len>1) 
		print("test0")
		l <- popsize_len
		legend <- "population size "
		results <- (foreach(pops.size=pop.size) %do% SOLVE(data, problem, pops.size, pr.elit, pr.cross , pr.mut, op.cross, op.mut, .time.limit))

	if (prelit_len>1) 
		print("test1")
		l <- prelit_len
		legend <- "elitism percentage "
		results <- (foreach(prs.elit=unlist(pr.elit)) %do% SOLVE(data, problem, pop.size, prs.elit, pr.cross , pr.mut, op.cross, op.mut, .time.limit))	

	if (prcross_len>1) 
		l <- prcross_len
		legend <- "probability crossover "
		results <- (foreach(prs.cross=c(0.4,0.6,0.8)) %do% SOLVE(data, problem, pop.size, pr.elit, prs.cross , pr.mut, op.cross, op.mut, .time.limit))	

	if (prmut_len>1) 
		print("test3")
		l <- prmut_len
		legend <- "probability mutation "
		results <- (foreach(prs.mut=pr.mut) %do% SOLVE(data, problem, pop.size, pr.elit, pr.cross , prs.mut, op.cross, op.mut, .time.limit))

	if (opcross_len>1) 
		print("test4")
		l <- opcross_len
		legend <- "operator crossover "
		results <- (foreach(ops.cross=op.cross) %do% SOLVE(data, problem, pop.size, pr.elit, pr.cross , pr.mut, ops.cross, op.mut, .time.limit))

	if (opmut_len>1) 
		print("test5")
		l <- opmut_len
		legend <- "operator mutation "
		results <- (foreach(ops.mut=op.mut) %do% SOLVE(data, problem, pop.size, pr.elit, pr.cross , pr.mut, op.cross, ops.mut, .time.limit))


	# desglose de resultados

	print("hola")

	aux_results <- list()
	matrix_start <- list()
	matrix_end <- list()
	n_its <- c()
	quarters <- c()
	quarter_results <- c()

	print("tests fuerza bruta")
	print(results)
	# print(results[[1]][[1]])
	# print(results[[2]][[1]])
	# print(results[[3]][[1]])


	for (i in 1:long) {
		result <- results[[i]][[1]]
		its <- results[[i]][[2]]

		matrix_start <- append(matrix_start, unlist(AUX_PLOT_DATA(data, result[1,1])))
		matrix_end <- append(matrix_end, unlist(AUX_PLOT_DATA(data, result[1,its])))

		aux_results <- append(aux_results, result)
		n_its <- c(n_its, its)
		quarters <- c(quarters, round(its*0.75))
		quarter_results <- c(quarter_results, result[2,round(its*0.75)])

	}

	print("resultados")
	print(long)
	print("adios")


	# vplayout <- function(x, y) viewport(layout.pos.row = x, layout.pos.col = y)

	# plot1 <- qplot(xlab="number of iterations reached", ylab="results", 1:n_its_1,unlist(res1[2,]), main="0.7 crossover")
	# plot2 <- qplot(xlab="number of iterations reached", ylab="results", 1:n_its_2,unlist(res2[2,]), main="0.8 crossover")
	# plot3 <- qplot(xlab="number of iterations reached", ylab="results", 1:n_its_3,unlist(res3[2,]), main="0.9 crossover")


	# plot10 <- qplot(xlab="number of iterations reached", ylab="results", q1:n_its_1,unlist(res1[2,][q1:n_its_1]), main="0.7 crossover")
	# plot11 <- qplot(xlab="number of iterations reached", ylab="results", q2:n_its_2,unlist(res2[2,][q2:n_its_2]), main="0.8 crossover")
	# plot12 <- qplot(xlab="number of iterations reached", ylab="results", q3:n_its_3,unlist(res3[2,][q3:n_its_3]), main="0.9 crossover")

	# plot4 <- qplot(xlab="X-axe", ylab="Y-axe",matrix1_start[1,], matrix1_start[2,], main="TSP[start] 0.7 crossover") + geom_line(color="steelblue")
	# plot5 <- qplot(xlab="X-axe", ylab="Y-axe",matrix2_start[1,], matrix2_start[2,], main="TSP[start] 0.8 crossover") + geom_line(color="steelblue")
	# plot6 <- qplot(xlab="X-axe", ylab="Y-axe",matrix3_start[1,], matrix3_start[2,], main="TSP[start] 0.9 crossover") + geom_line(color="steelblue")

	# plot7 <- qplot(xlab="X-axe", ylab="Y-axe",matrix1_end[1,], matrix1_end[2,], main="TSP[end] 0.7 crossover") + geom_line(color="darkblue")
	# plot8 <- qplot(xlab="X-axe", ylab="Y-axe",matrix2_end[1,], matrix2_end[2,], main="TSP[end] 0.8 crossover") + geom_line(color="darkblue")
	# plot9 <- qplot(xlab="X-axe", ylab="Y-axe",matrix3_end[1,], matrix3_end[2,], main="TSP[end] 0.9 crossover") + geom_line(color="darkblue")


	# # 9 figures arranged in 3 rows and 3 columns

	# grid.arrange(plot1,plot2, plot3, plot10, plot11, plot12, ncol=3, nrow=2)
	# #qplot(1:n_its,results[[chosen]])
	# grid.newpage()
	# grid.arrange(plot4, plot7, plot5, plot8, plot6, plot9, ncol=2, nrow=3)

	# print(res1[,n_its_1])
	# print(res2[,n_its_2])
	# print(res3[,n_its_3])



	stopCluster(cl)
	
}	


AUX_PLOT_DATA <- function(data, perm){

	x <- c()
	y <- c()
		for (i in perm) {
			x <- c(x,data[i,1]) # x
			y <- c(y,data[i,2]) # y
		}
	
	rbind(x,y)

}
