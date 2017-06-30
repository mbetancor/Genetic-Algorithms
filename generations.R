source('crossovers.R')
source('mutations.R')

generate_offspring <- 
function(population, .op.cross="pmx", .op.mut="sim", .times=20, .pr.cross = 0.6 , .pr.mut = 0.1) {
	next_gen <- population
	len <- length(population)
	.times <- round(.times/2)

	if (.op.cross=="ap")
		for (i in 1:(2*.times)){
			indexes 	 <- sample(1:len,2,replace=F)
			individual_1 <- population[[indexes[1]]]
			individual_2 <- population[[indexes[2]]]

			# crossover + mutation
			individual 	<- MUTATION(CROSSOVER(individual_1,individual_2, .op.cross, .pr.cross), .op.mut, .pr.mut)	
			next_gen 	<- c(next_gen,individual)
		}
	else 
		for (i in 1:.times){
			indexes		 <- sample(1:len,2,replace=F)
			individual_1 <- population[[indexes[1]]]
			individual_2 <- population[[indexes[2]]]

			# crossover
			ans <- CROSSOVER(individual_1, individual_2, .op.cross, .pr.cross)
			individual_1 <- ans[1,]
			individual_2 <- ans[2,]
			
			# mutation
			individual_1 <- MUTATION(individual_1, .op.mut, .pr.mut)
			individual_2 <- MUTATION(individual_2, .op.mut, .pr.mut)
				
			next_gen 	 <- c(next_gen,list(individual_1),list(individual_2))
		}

	next_gen
}