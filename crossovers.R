CROSSOVER <- 
function(str1, str2,.op.cross ="pmx", .pr.cross=0.6) 
{
	len <- length(str1)
	size <-  len * .pr.cross 
	off1 <- integer(len)
	off2 <- integer(len)

	if(.op.cross == "pmx") 
		CROSS <- (crossover.pmx(str1,str2, off1, off2, size, len))
	if(.op.cross == "ox1") 
		CROSS <- (crossover.ox1(str1,str2, off1, off2, size, len))
	if(.op.cross == "ox2") 
		CROSS <- (crossover.ox2(str1,str2, off1, off2, size, len))
	if(.op.cross == "pos") 
		CROSS <- (crossover.pos(str1,str2, off1, off2, size, len))
	if(.op.cross == "ap2") 
		CROSS <- (rbind(crossover.ap(str1,str2, off1, off2, len),(crossover.ap(str2,str1, off1,off2, len))))
	if(.op.cross == "ap") 
		CROSS <- (rbind(crossover.ap(str1,str2, off1, off2, len)))

	return(CROSS)
}

crossover.pmx <- 
function(str1, str2, off1, off2, size, len) 
{
	rc1 <- sample(1:(len-size),1)
	rc2 <- sample((rc1+1):len,1)

	map_section_1 <- integer(rc2-rc1)
	map_section_2 <- integer(rc2-rc1)

	ms1 <- str1[rc1:rc2]
	ms2 <- str2[rc1:rc2]

	off1[rc1:rc2] <- ms2
	off2[rc1:rc2] <- ms1

	order <- order(ms1)
	index <- 1

	for (i in order){
		map_section_1[index] <- ms1[i]
		map_section_2[index] <- ms2[i]
		index <- index+1
	}

	map_sections_arr <- remove_duplicates(map_section_1, map_section_2)

	ms1 <- map_sections_arr[,1]
	ms2 <- map_sections_arr[,2]

	for (i in (1:len)[-(rc1:rc2)]) {
		item_list_1 <- which(ms2==str1[i])
		item_list_2 <- which(ms1==str2[i])

		if (length(item_list_1)){
			off1[i] <- ms1[item_list_1]
		}

		else {
			off1[i] <- str1[i]
		}

		if (length(item_list_2)){
			off2[i] <- ms2[item_list_2]
		}

		else {
			off2[i] <- str2[i]
		}
	}

	return(rbind(off1,off2))
}


crossover.ox1 <-
function(str1, str2, off1, off2,size, len) 
{

	size <- len - size
	rc1 <- sample(1:(len-size),1)
	rc2 <- sample((rc1+1):len,1)

	off1[rc1:rc2] <- str1[rc1:rc2]
	off2[rc1:rc2] <- str2[rc1:rc2]

	path_sequence <- append((rc2:len)[(rc2:len)>rc2], (1:rc1)[(1:rc1)<rc1])

	path1 <- append(str2[(rc2:len)[(rc2:len)>rc2]],str2[(1:rc2)])
	path2 <- append(str1[(rc2:len)[(rc2:len)>rc2]],str1[(1:rc2)])

	path1 <- path1[!path1 %in% str1[rc1:rc2]]
	path2 <- path2[!path2 %in% str2[rc1:rc2]]

	index <- 1
	for (i in path_sequence){
		off1[i] <- path1[index]
		off2[i] <- path2[index]
		index <- index+1	
	}

	return(rbind(off1,off2))

}

crossover.ox2 <- 
function(str1, str2, off1, off2,size, len) 
{
	rc1 <- sample(1:(len-size),1)
	rc2 <- sample((rc1+2):len,1)

	off1[1:rc1] <- str1[1:rc1]
	off1[rc2:len] <- str1[rc2:len]
	
	off2[1:rc1] <- str2[1:rc1]
	off2[rc2:len] <- str2[rc2:len]

	path_sequence <- c( (rc1+1): (rc2-1))

	path1 <- str2
	path2 <- str1

	path1 <- path1[!path1 %in% off1]
	path2 <- path2[!path2 %in% off2]

	index <- 1
	for (i in path_sequence){
		off1[i] <- path1[index]
		off2[i] <- path2[index]
		index <- index+1	
	}

	return(rbind(off1,off2))

}

crossover.pos <- 
function(str1, str2, off1, off2, size, len)
{

	n_random_points <- sample(2:size,1)
	random_points <- sample(1:len, n_random_points, replace=F)
	path_sequence <- (1:len)[!(1:len) %in% random_points]
	acc1 <- c()
	acc2 <- c()
	index <- 1

	for (i in random_points){
		off1[i] <- str2[i]
		off2[i] <- str1[i]
		acc1 <- c(acc1, str1[i])
		acc2 <- c(acc2, str2[i])
	}

	seq1 <- str1[!str1 %in% acc2]
	seq2 <- str2[!str2 %in% acc1]

	for (i in path_sequence){
		off1[i] <- seq1[index]
		off2[i] <- seq2[index]
		index <- index+1 
	}
	
	return(rbind(off1,off2))

}

crossover.ap <- 
function(str1, str2, .pr.cross, len) 
{

	off <- c()
	turn <- 0
	parents <- rbind(str1,str2)
	while(length(off)<len){

		p <- parents[turn+1,]
		while(p[1] %in% off) {
			p <- p[-1]
		}
		
		off <- c(off,p[1])
		turn <- ((turn+1) %% 2)	
	}

	return(off)

}

remove_duplicates <- 
function(elements, sorted_list) 
{

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

	cbind(elements, sorted_list)

}