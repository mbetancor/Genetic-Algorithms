## ----installPackage,eval=FALSE-------------------------------------------
## install.packages('rscala')

## ----installScala,eval=FALSE---------------------------------------------
## rscala::scalaInstall()

## ----scalaInfo,eval=FALSE------------------------------------------------
## rscala::scalaInfo(verbose=TRUE)

## ----library-------------------------------------------------------------
library('rscala')

## TODO: LOCAL SEARCH, PARALELLISM


BRUTE_FORCE <- function (data, nperms){
     
    t <- proc.time()

    nrows <- nrow(data)
    perm <- RANDOM_POPULATION(nperms,nrows)
    distances <- AUX_GET_DISTANCE(perm,data)
    
    
    min_d <- distances[1]
    index <- 1
    for (i in 2:length(distances)){
        if (distances[i]<min_d){
            min_d <- distances[i]
            index <- i
        }
    }
    
    print(proc.time()-t)
    append(list(min_d),perm[i])
    
}


SOLVE <- function(data, n_off, type, method, limit, sa_per){

	t <- proc.time()

	nrows <- nrow(data)
	perm <- RANDOM_POPULATION(n_off,nrows)
	distances <- AUX_GET_DISTANCE(perm,data)

	# calcula fitness
	min <- min(distances)
	f <- (min(distances)+max(distances))/2

	arr <- character(0)

	for (i in 1:n_off){
		if (distances[i]<=f){
			arr <- append(arr, perm[i])
		}
	}

	off <- AUX_GENERATE_OFFSPRING(type, method,limit,arr,data,min, sa_per)
    
    d <- AUX_GET_DISTANCE(off,data)
    
    min_d <- d[1]
    index <- 1
    for (i in 2:length(d)){
        if (d[i]<min_d){
            min_d <- d[i]
            index <- i
        }
    }
    
    print(proc.time()-t)
    append(list(min_d),off[i])
    
	
}


PLOT_TSP <- function(data,perm){

	x <- character(0)
	y <- character(0)

	for (i in perm){
		x <- append(x,data[i,1])
		y <- append(y,data[i,2])
	}

	plot.new()
	plot.window(range(x, na.rm = TRUE), range(y, na.rm = TRUE))
	polypath(x,y, col="grey", rule="winding")

}

AUX_GENERATE_OFFSPRING <- function(type, method, limit, parents, data, f, sa_per){
    
    i <-1

	while(length(parents)>1 && i<limit){
        p <- sample(1:length(parents),2,replace=F)
        n1 <- min(p)
        n2 <- max(p)
            
        p1 <- unlist(parents[n1])
        p2 <- unlist(parents[n2])

 		parents <- parents[-n2][-n1]

        if (type=="crossover"){
        	
            children <- method(p1,p2)

			c1 <- children[1]
			c2 <- children[2]
			
			d1 <- AUX_GET_DISTANCE(c1,data)
			d2 <- AUX_GET_DISTANCE(c2,data)
            
            distance <- min(d1,d2)
            if (distance==d1){
                child <- c1
            } else {
                child <- c2
            }
            
        } else {
            if (type=="crossover1") {
            	child <- list(method(p1,p2))
            } else {
            	child <- list(method(p1))
            }

        	distance <- AUX_GET_DISTANCE(child,data)

    	}
        
        ran <- sample(1:100,1)
            
        if (distance<=f || ran <= sa_per){
            parents <- append(child,parents)
        } else {
            parents <- append(parents,append(list(p1),list(p2)))
        }
            
        if (distance<=f){
            f <- (f +distance)/2
        }

        i <- i +1

	}

	parents

}

AUX_GET_DISTANCE <- function (perm,data){

	n_off <- length(perm)
	
	
	nrows <- nrow(data)
	distances <- integer(n_off)

	for(i in 1:n_off){
		d <- 0
		k <- 1
		while(k<nrows){
			n <- unlist(perm[i])[k]
			m <- unlist(perm[i])[k+1]
			x1 <- data[n,1]
			x2 <- data[m,1]
			y1 <- data[n,2]
			y2 <- data[m,2]
			d <- d + EUC_DIST(x1,y1,x2,y2)
			k <- k+1
		}
        
		distances[i] <- d
	}

	distances
}

RANDOM_POPULATION <- function(n_off,nrows){

	ans <- character(0)

    for (i in 1:n_off){
		ans <- append(ans,list(sample(nrows)))
	}

	ans

}

EUC_DIST <- function(x1,y1,x2,y2){

	x <- (x1-x2)**2
	y <- (y1-y2)**2
	(x+y)**(0.5)
}



PMX <- function(str1, str2) {
	len <- length(str1)
	mapsections <- sample(2:(len-1),2, replace=F)
	index1 <- min(mapsections)
	index2 <- max(mapsections)
	replacements <- cbind(str1[index1:index2],str2[index1:index2])

    
	end <- len-index2
	child1 <- append(integer(index1-1),append(str2[index1:index2],integer(end)))
	child2 <- append(integer(index1-1),append(str1[index1:index2],integer(end)))
    
    path <- append(1:(index1-1),(index2+1):len)
   
    
	replacements <- AUX_REMOVE_DUP(replacements)
    
    l <- length(replacements[,2])
    
    for (i in path){
		if (is.element(str1[i],child1)){
            k <- 1
            while (k<=l){
				if (str1[i]==replacements[k,2]){
					child1[i] <- replacements[k,1]
                    break
				}
                k <- k+1
			}
		} else {
			child1[i] <- str1[i]
		}

        if (is.element(str2[i],child2)){
            k <- 1
            while (k<=l){
                if (str2[i]==replacements[k,1]){
                    child2[i] <- replacements[k,2]
                    break
                }
                k <- k+1
            }
        } else {
            child2[i] <- str2[i]
        }
    }

	append(list(child1), list(child2))
	
}

AUX_REMOVE_DUP <- function(arr){
    
    i <- 1
    while (i<length(arr[,2])){
        
        elem1 <- arr[i,1]
        j <- i+1
        
        while(j<=length(arr[,2])){
            if (elem1==arr[j,2]){
                aux <- arr[j,1]
                arr <- arr[-j,]
                # caso número incorrecto de dimensiones
                if (length(arr)==2){
                    arr <- cbind(arr[1],arr[2])
                }
                arr[i,1] <- aux
                break;
            }
            j<- j+1
        }
        
        elem2 <- arr[i,2]
        
        j <- i+1
        while(j<=length(arr)/2){
            if (elem2==arr[j,1]){
                aux <- arr[j,2]
                arr <- arr[-j,]
                # caso número incorrecto de dimensiones
                if (length(arr)==2){
                    arr <- cbind(arr[1],arr[2])
                }
                arr[i,2] <- aux
                break;
            }
            j <- j+1
        }
        
        i<- i+1
        
    }
    
    arr
}


AUX_REORDER_LIST <- function(start1, end1, start2, end2, str){
	
	f1 <- str[start1:end1]
	f2 <- str[start2:end2]
	f <- append(f1,f2)

	if (end1 <= start1){
		f <- f2
	}
	
	if (start2>=length(str)){
		f <- f1
	}

	f
}

# Cycle CrossOver (CX)

CX <- function(str1, str2) {
	len <- length(str1)
	off <-  integer(len)
	seeked <- str2[1]
	elem <- str1[1]
	off[1] <- elem
	p <- 0

	while(length(off[off==0])>0){

		i <- AUX_GET_FZERO(off)

		if (seeked ==0){
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

AUX_GET_FZERO <- function(str){

	i <- 1
	len <- length(str)
	while(i<len){
		if (str[i] == 0) break;
		i <- i+1
	}

	i
}


# OX1

OX1 <- function(str1, str2) {
	len <- length(str1)
	off1 <- integer(len)
	off2 <- integer(len)

	mapsections <- sample(2:(len-1),2, replace=F)
	index1 <- min(mapsections)
	index2 <- max(mapsections)

	off1[index1:index2] <- str1[index1:index2]
	off2[index1:index2] <- str2[index1:index2]
	path1 <- append(str2[(index2+1):len],str2[1:(index2)])
	path2 <- append(str1[(index2+1):len],str1[1:(index2)])

	i <- index2+1
	n <- 1 #off1
	m <- 1 #off2

	while(i<=len){

		while (is.element(path1[n],off1)){
			n <- n+1
		}
		off1[i] <- path1[n]

		if (is.element(path2[m],off1)){
			m <- m+1
		}
		off2[i] <- path2[n]

		i<- i+1
	}

	for (i in 1:(index1-1)){

		while (is.element(path1[n],off1)){
			n <- n+1
		}
		off1[i] <- path1[n]

		if (is.element(path2[m],off1)){
			m <- m+1
		}
		off2[i] <- path2[n]

	}
	
	
    append(list(off1),list(off2))

}



# DM


DM <- function(str) {
	len <- length(str)
	mapsections <- sample(2:(len-1),2, replace=F)
	index1 <- min(mapsections)
	index2 <- max(mapsections)
	subtour <- str[index1:(index2-1)]
	str_aux <- append(str[0:(index1-1)], str[(index2):len])
	ran <- sample(1:length(str_aux),1)
	append(str_aux[1:(ran-1)],append(subtour,str_aux[ran:length(str_aux)]))

}


# IVM

IVM <- function(str) {
	len <- length(str)
	mapsections <- sample(2:(len-1),2, replace=F)
	index1 <- min(mapsections)
	index2 <- max(mapsections)
	subtour <- str[(index2-1):(index1)]
	str_aux <- append(str[0:(index1-1)], str[(index2):len])
	ran <- sample(1:length(str_aux),1)
	append(str_aux[0:(ran-1)],append(subtour,str_aux[ran:length(str_aux)]))
}

# ISM

ISM <- function(str) {
	len <- length(str)
	mapsections <- sample(2:(len-1),2, replace=F)
	index1 <- min(mapsections)
	index2 <- max(mapsections)
	start <- append(str[0:(index1-1)],str[(index1+1):(index2)])
	end <- append(str[index1],str[(index2+1):len])
	
	append(start,end)
}
	

AUX_GET_FZERO <- function(str){

	i <- 1

	while(i<=length(str)){
		if (str[i] == 0) break;
		i <- i+1
	}

	i
}
