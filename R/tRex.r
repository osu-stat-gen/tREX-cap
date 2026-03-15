

.nz <- function(Contact, bias, n, repn, repb, argv0, argv1, argv2, nHL, vepsilon, thinning, gear, result){
    posterior=.C("nz",
    as.double(Contact),
    as.double(bias),
    as.integer(n),
    as.integer(repn),
    as.integer(repb),
    as.double(argv0),
    as.double(argv1),
    as.double(argv2),
    as.integer(nHL),
    as.double(vepsilon),
    as.integer(thinning),
    as.integer(gear),
    result=double(6*10000+n*30000+n+n*(n-1)*0.5+4),
    PACKAGE="tRexCAP")
}


.nz2 <- function(Contact, n, repn, repb, argv0, argv1, argv2, nHL, vepsilon, thinning, gear, result){
    posterior=.C("nz2",
    as.double(Contact),
    as.integer(n),
    as.integer(repn),
    as.integer(repb),
    as.double(argv0),
    as.double(argv1),
    as.double(argv2),
    as.integer(nHL),
    as.double(vepsilon),
    as.integer(thinning),
    as.integer(gear),
    result=double(6*10000+n*30000+n+n*(n-1)*0.5+4),
    PACKAGE="tRexCAP")
}

.bn <- function(Contact, bias, n, repn, repb, argv0, argv1, argv2, nHL, vepsilon, thinning, gear, result){
    posterior=.C("bn",
    as.double(Contact),
    as.double(bias),
    as.integer(n),
    as.integer(repn),
    as.integer(repb),
    as.double(argv0),
    as.double(argv1),
    as.double(argv2),
    as.integer(nHL),
    as.double(vepsilon),
    as.integer(thinning),
    as.integer(gear),
    result=double(2*10000+12*10000+n*30000+4),
    PACKAGE="tRexCAP")
}


.bn2 <- function(Contact, n, repn, repb, argv0, argv1, argv2, nHL, vepsilon, thinning, gear, result){
    posterior=.C("bn2",
    as.double(Contact),
    as.integer(n),
    as.integer(repn),
    as.integer(repb),
    as.double(argv0),
    as.double(argv1),
    as.double(argv2),
    as.integer(nHL),
    as.double(vepsilon),
    as.integer(thinning),
    as.integer(gear),
    result=double(12*10000+n*30000+2),
    PACKAGE="tRexCAP")
}


.rotate=function(rSi){
    S=t(matrix(nrow=3,unlist(rSi)))
    ns=nrow(S)
	x1=S[1,1]
	y1=S[1,2]
	z1=S[1,3]
	
	S[,1]=S[,1]-x1
	S[,2]=S[,2]-y1
	S[,3]=S[,3]-z1
	
	x=S[ns,1]
	y=S[ns,2]
	z=S[ns,3]
	d=sqrt(x^2+y^2)
	
	Rz=matrix(ncol=3,nrow=3,0)
	Rz[1,1]=x/d #cos
	Rz[1,2]=y/d #sin
	Rz[2,1]=-y/d
	Rz[2,2]=x/d
	Rz[3,3]=1
	SR=S%*%t(Rz)
	
	x=SR[ns,1]
	y=SR[ns,2]
	z=SR[ns,3]
	d=sqrt(x^2+y^2+z^2)
	
	R2x=matrix(ncol=3,nrow=3,0)
	R2x[1,1]=x/d #cos
	R2x[1,3]=z/d #sin
	R2x[3,1]=-z/d
	R2x[3,3]=x/d
	R2x[2,2]=1
	S=SR%*%t(R2x)
	
	x=S[2,1]
	y=S[2,2]
	z=S[2,3]
	d=sqrt(y^2+z^2)	
	R2x=matrix(ncol=3,nrow=3,0)
	R2x[2,2]=z/d
	R2x[2,3]=-y/d 
	R2x[3,2]=y/d
	R2x[3,3]=z/d
	R2x[1,1]=1

	
	SR=S%*%t(R2x)
    
    if(SR[3,2]<0){
        SR[,2]=-SR[,2]
    }
	
	as.vector(t(SR))
}


trex<-
function(contact,bias=NULL,mcmc=10000,burn=10000,jump.beta1,jump.x=0.1,jump.u=0.1,jump.cov1=0.1,jump.cov2=0.1,leapfrog.L=10,leapfrog.e=0.001,method.type=c("tRex","tPAM","bn"))
{
    n = nrow(contact)
    switch(method.type,
    tRex={
        if(!is.null(bias)){
            if(n != nrow(bias)) stop("tRex : The number of loci is different from the number of bias.\n")
            if(ncol(bias) != 3) stop("tRex : The number of columns is not 3.\n")
        }
        if(max(abs(contact-t(contact)))!=0) stop("tRex : The contact matrix is not symmetric.\n")
        if(is.null(jump.x)) stop("tRex : A jumping rule for X is not specified.\n")
        if(is.null(jump.u)) stop("tRex : A jumping rule for U is not specified.\n")
    },
    tPAM={
        if(!is.null(bias)){
            if(n != nrow(bias)) stop("tPAM : The number of loci is different from the number of bias.\n")
            if(ncol(bias) != 3) stop("tPAM : The number of columns is not 3.\n")
        }
        if(max(abs(contact-t(contact)))!=0) stop("tPAM : The contact matrix is not symmetric.\n")
    },
    bn={
        if(!is.null(bias)){
            if(n != nrow(bias)) stop("Bayesian Nonparametric : The number of loci is different from the number of bias.\n")
            if(ncol(bias) != 3) stop("Bayesian Nonparametric : The number of columns is not 3.\n")
        }
        if(max(abs(contact-t(contact)))!=0) stop("Bayesian Nonparametric : The contact matrix is not symmetric.\n")
    }
    )
    
    if(!is.null(bias)){
        bias=log(bias)
        m.bias=colMeans(bias)
        bias2=bias
        bias2[,1]=bias[,1]-m.bias[1]
        bias2[,2]=bias[,2]-m.bias[2]
    }
    
    thinning=mcmc/10000
    if(thinning %%1 != 0) stop("Thinning is not integer.")

    
    Contact=rep(0,n*(n-1)*0.5)
    k=1
    for(i in 1:(n-1)){
        for(j in (i+1):n){
            Contact[k]=contact[i,j]
            k=k+1
        }
    }
    
    if(is.null(bias)){
        if(method.type=="tRex") output=.nz2(Contact, n, mcmc, burn, jump.beta1, jump.x, jump.u, leapfrog.L, leapfrog.e, thinning,1)
        else if(method.type=="tPAM") output=.nz2(Contact, n, mcmc, burn, jump.beta1, jump.cov1, jump.cov2, leapfrog.L, leapfrog.e, thinning,0)
        else if(method.type=="bn") output=.bn2(Contact, n, mcmc, burn, jump.beta1, jump.cov1, jump.cov2, leapfrog.L, leapfrog.e, thinning,0)
        
    }else{
        if(method.type=="tRex") output=.nz(Contact, as.vector(t(bias2)), n, mcmc, burn, jump.beta1, jump.x, jump.u, leapfrog.L, leapfrog.e, thinning,1)
        else if(method.type=="tPAM") output=.nz(Contact, as.vector(t(bias2)), n, mcmc, burn, jump.beta1, jump.cov1, jump.cov2, leapfrog.L, leapfrog.e, thinning,0)
        else if(method.type=="bn") output=.bn(Contact, as.vector(t(bias2)), n, mcmc, burn, jump.beta1, jump.cov1, jump.cov2, leapfrog.L, leapfrog.e, thinning,0)
    }
    
    R=output[["result"]]
    
    S=matrix(nrow=10000,ncol=3*n,0)
    
    k=1
    for(i in 1:10000){
        for(j in 1:(3*n)){
            S[i,j]=R[k]
            k=k+1
        }
    }
    
    
    # loci=apply(S,1,.rotate)
    loci = S
    # loci=t(loci)
    
    X=matrix(ncol=1,nrow=n,0)
    U=matrix(ncol=n,nrow=n,0)

    if(method.type=="tRex"){
        Param=matrix(nrow=10000,ncol=5,0)
        colnames(Param)=c("beta1","cov1","cov2","sigma.x","sigma.u")
        for(i in 1:10000){
            for(j in 1:5){
                Param[i,j]=R[k]
                k=k+1
            }
        }

        for(i in 1:n){
            X[i,1]=R[k];k=k+1
        }
        for(i in 1:(n-1)){
            for(j in (i+1):n){
                U[i,j]=R[k];U[j,i]=R[k];k=k+1
            }
        }
    }else if(method.type=="tPAM"){
        Param=matrix(nrow=10000,ncol=3,0)
        colnames(Param)=c("beta1","cov1","cov2")
        for(i in 1:10000){
            for(j in 1:3){
                Param[i,j]=R[k]
                k=k+1
            }
        }
    }else if(method.type=="bn"){
        if(!is.null(bias)){
            Param=matrix(nrow=10000,ncol=13,0)
            colnames(Param)=c("cov1","cov2","b1","b2","b3","b4","b5","b6","b7","b8","b9","b10","b11")
            for(i in 1:10000){
                for(j in 1:2){
                    Param[i,j]=R[k]
                    k=k+1
                }
            }
            for(i in 1:10000){
                for(j in 3:13){
                    Param[i,j]=R[k]
                    k=k+1
                }
            }
        }else{
            Param=matrix(nrow=10000,ncol=11,0)
            colnames(Param)=c("b1","b2","b3","b4","b5","b6","b7","b8","b9","b10","b11")
            for(i in 1:10000){
                for(j in 1:11){
                    Param[i,j]=R[k]
                    k=k+1
                }
            }
        }
    }
    acceptance=c(R[k],R[(k+1)],R[(k+2)],R[(k+3)])
    k=k+4
    lk=matrix(nrow=10000,ncol=1,0)
    for(i in 1:10000){
        lk[i,1]=R[k]
        k=k+1
    }
    list(loci,Param,X,U,acceptance,lk)
}



#' @title get posterior samples of coordinates of a specific locus.
#' @param psample An output of tRex().
#' @param i The locus to get coordinates. 
#' @return A Matrix of 10000 rows. 1st, 2nd, and 3rd colums corresponds to x,y, and z coordinates. Each rows corresponds to a posterior sample  of coordinates.
#' @export
#' 
get.sample=function(psample,i){
	where=c(3*(i-1)+1,3*(i-1)+2,3*i)
	S=psample[[1]]
	Si=S[,where]
	colnames(Si)=c("X","Y","Z")
	Si
}

#' @title get posterior samples of parameters.
#' @param psample An output of tRex().
#' @return A matrix of 10000 rows and 5 columns of \beta_1, cov_1, cov_2, \sigma_x, and \sigma_u.
#' @export
#' 
get.param=function(psample){
	Param=psample[[2]]
}	
	
#' @title Summarizing posterior samples of loci coordinates.
#' @description Posterior coordinates can be summarized using either posterior mode or posterior mean.
#' @param psample An output of tRex().
#' @param method which method to summarize the coordinates. One of "Mode", "Mean", or "MAP".
#' @return A matrix of three columns and n rows, where n is the number of loci.
#' @export
#' 
summarize.struct=function(psample,method=c("Mode","Mean","MAP")){
	which=match.arg(method)
	S=psample[[1]]
	lk=psample[[6]]
	n=ncol(S)/3
	end=3*n-2
	switch(which,
		"Mode"={	
			mmestloci=posterior.mode(mcmc(S))
			mmestloci[c(1,2,3,5,end+1,end+2)]=0
			tmaxS=matrix(nrow=3,unlist(mmestloci))
			maxS=t(tmaxS)
		},
		"Mean"={
			#mmestloci=colMeans(S)
			#mmestloci[c(1,2,3,5,end+1,end+2)]=0
			#tmaxS=matrix(nrow=3,mmestloci)
			#maxS=t(tmaxS)
			calmat=matrix(nrow=10000,ncol=n*(n-1)*0.5,0)
			for(i in 1:10000){
				ai=matrix(ncol=3,nrow=n,S[i,],byrow=T)
				calmat[i,]=as.vector(dist(ai,upper=T))
			}

			standard=colMeans(calmat)
			delta=vector()
			for(i in 1:nrow(calmat)){
				delta[i]=sum(abs(calmat[i,]-standard))
			}
			id=which.min(delta)
			maxS=matrix(ncol=3,byrow=T,S[id,])
		},
		"MAP"={
			id=which.max(lk)
			maxS=matrix(ncol=3,byrow=T,S[id,])
		}
	)
}

#' @title Draw 3D structure
#' @description Using a rgl package, this function provides a visualization function of the 3D structure
#' @param coordinates A matrix of three coloumns correspoding to three coordinates.
#' @return The figure.
#' @export
draw.struct=function(coordinates){
	x=coordinates[,1]
	y=coordinates[,2]
	z=coordinates[,3]
	n=nrow(coordinates)

	df=data.frame(x=x,y=y,z=z)
	df$color="yellow"
	df$size=2


	plot3d(df$x, df$y, df$z, col=df$color, size=df$size, type='s',lwd=1,xlab="",ylab="",zlab="")

	for(i in 1:(n-1)){
		segments3d(df$x[i:(i+1)],df$y[i:(i+1)],df$z[i:(i+1)],col="brown",lwd=3)
	}
}


mctrex <- function(k, bias = NULL, contact, cutlist, save_mcmc = FALSE) {
  # partition the matrix
  
  start <- cutlist[[k]][1]
  end <- cutlist[[k]][2]
  
  mat <- contact[start:end, start:end]
  if(!is.null(bias)){
      bias_term = bias[start:end,]
  } else{
      bias_term = NULL
  }
  # run tRex mcmc=100000, burn=100000
  psample <- trex(contact = mat, bias = bias_term, mcmc = 100000, burn = 100000, jump.beta1 = 0.05, jump.x = 0.2, jump.u = 0.2, leapfrog.L = 10, leapfrog.e = 0.001, method.type = "tRex")
  param <- get.param(psample)
  
  str1 <- summarize.struct(psample, method = "Mean")
  # X is the random effect
  # same as X in data_generation.r
  
  if (save_mcmc){
    return(list(coords = str1, X = psample[[3]], params = colMeans(param), mcmc = psample))
  } else {
    return(list(coords = str1, X = psample[[3]], params = colMeans(param)))      
  }
}


Rx_apply <- function(S, theta) {
  R <- matrix(0, ncol = 3, nrow = 3)
  R[1, 1] <- 1
  R[2, 2] <- cos(theta)
  R[2, 3] <- -sin(theta)
  R[3, 2] <- sin(theta)
  R[3, 3] <- cos(theta)
  a <- S %*% t(R)
  # a1 <- t(R %*% t(S))
}

Ry_apply <- function(S, theta) {
  R <- matrix(0, ncol = 3, nrow = 3)
  R[1, 1] <- cos(theta)
  R[1, 3] <- sin(theta)
  R[2, 2] <- 1
  R[3, 1] <- -sin(theta)
  R[3, 3] <- cos(theta)
  a <- S %*% t(R)
  # a <- t(R %*% t(S))
}

Rz_apply <- function(S, theta) {
  R <- matrix(0, ncol = 3, nrow = 3)
  R[1, 1] <- cos(theta)
  R[1, 2] <- -sin(theta)
  R[2, 1] <- sin(theta)
  R[2, 2] <- cos(theta)
  R[3, 3] <- 1
  a <- S %*% t(R)
  # a <- t(R %*% t(S))
}

# fast way to calculate vector distance
mat_dist <- function(x, y) {
  sqrt(outer(rowSums(x^2), rowSums(y^2), "+") - tcrossprod(x, 2 * y))
}

# beta_min <- function(beta, llambdax, y12, ddd) {
#   y <- sum(log(ddd) * ddd^beta * exp(llambdax)) - sum(y12 * log(ddd))
#   return(y)
# }

# input a vector of break points and the dimension of the contact matrix
# output a list of the start and end breakpoints 

get.breakpoints <- function(breaks=NULL, n, block_size = 40, noverlap = 1){
  numloci = n
  
  # if breaks is null, then we assume that the block_size is equal
  if(is.null(breaks)){
    nblock <- ceiling(numloci / block_size)
    cutlist <- list()
    for(k in 1:nblock){
      if (k == 1) {
        start <- 1
        end <- block_size + noverlap
      } else if (k == nblock) {
        start <- (k - 1) * block_size - noverlap + 1
        end <- n
      } else {
        start <- (k - 1) * block_size - noverlap + 1
        end <- k * block_size + noverlap
      }
      cutlist[[k]] <- c(start, end)
    }
  } else {
    lenbreaks <- length(breaks)
    #if(breaks[lenbreaks]>n){}  to check if we have a valid break, breaks[lenbreaks]=n, also breaks[1]=1
    cutlist <- list()
    for(k in 1:(lenbreaks-1)){
      if (k == 1) {
        start <- breaks[k]
        end <- breaks[k+1] -1 + noverlap
      } else if (k == lenbreaks-1) {
        start <- breaks[k] - noverlap 
        end <- breaks[k+1]
      } else {
        start <- breaks[k] - noverlap 
        end <- breaks[k+1] -1 + noverlap
      }
      cutlist[[k]] <- c(start, end)
    }
  }
  return(cutlist)
}


# Run the cut part of cut and paste
Cut <- function(contact, bias = NULL, breaks = NULL, block_size = 40, noverlap = 1, CPU, save_mcmc = FALSE){
  n = ncol(contact)
  cutlist = get.breakpoints(breaks=breaks, n=n, block_size=block_size, noverlap=noverlap)
  nblock = length(cutlist)
  
  
  # Check if the operating system is Windows
  if (.Platform$OS.type == "windows") {
    # Load doParallel for Windows
    if (requireNamespace("doParallel", quietly = TRUE)) {
      library(doParallel)
      cl <- makeCluster(CPU)
      registerDoParallel(cl)
      on.exit(stopCluster(cl))  # Ensure the cluster stops after use
    } else {
      warning("doParallel is required for parallel processing on Windows, but it is not installed.")
    }
  } else {
    # Load doMC for Unix-like systems (Linux, macOS)
    if (requireNamespace("doMC", quietly = TRUE)) {
      library(doMC)
      registerDoMC(cores = CPU)
    } else {
      warning("doMC is required for parallel processing on Unix-like systems, but it is not installed.")
    }
  }
  
  
  result <- foreach(k = 1:nblock) %dopar% {
    mctrex(k = k, bias = bias, contact = contact, cutlist = cutlist, save_mcmc = save_mcmc) 
  }
  
  return(list(cutlist, result, noverlap))  # double check
}

# run the paste part of cut and paste
Paste_orig <- function(contact, cutresult, CPU){
  cutlist <- cutresult[[1]] # these are the breakpoints
  result <- cutresult[[2]] # these are the cuts
  noverlap <- cutresult[[3]] # this is noverlap used for the cuts
  
  y2 <- contact 
  N <- nrow(y2)
  
  sr <- 1
  # initial block
  start <- 1    #shouldn't need this for the initial block:  max(1, (sr - 1) * block_size - noverlap + 1)
  likelihood <- NULL
  y11 <- y2[start:nrow(y2), start:nrow(y2)]
  block1_end <- floor(length(cutlist)/2)-1  #floor(N / block_size / 2) - 1
  k <- 1
  for (k in sr:block1_end) {
    cat("k =", k, "\n")
    if (k == sr) {
      S1 <- result[[k]]$coords ## L1
      X1 <- result[[k]]$X[, 1]
    }
    
    # remove M1 and M2 margins, fix the last of M1 and first of M2 both to be (0,0,0)
    n1 <- nrow(S1)
    S1 <- t(t(S1) - S1[(n1 - noverlap + 1), ])
    S2 <- result[[(k + 1)]]$coords ## L2
    S2 <- S2[(noverlap + 1):nrow(S2), ]
    S2 <- t(t(S2) - S2[1, ])
    n2 <- nrow(S2)
    X2 <- result[[(k + 1)]]$X[, 1]
    
    fixed_point <- n1 - noverlap + 1
    if (fixed_point > N) break
    cat("Fixed Point", n1 - noverlap + 1, "\n")
    y12 <- y11[1:(fixed_point - 1), fixed_point:(fixed_point - 1 + n2)]
    fixed_point <- n1 - noverlap + 1
    llambdax <- exp(outer(X1[1:(n1 - noverlap)], X2[-(1:noverlap)], FUN = "+"))
    n1 <- nrow(S1)
    S1 <- S1[1:(n1 - noverlap), ]
    # running APG in C++
    res <- minimizer(50, y12,S1,S2,llambdax, threads=CPU)
    # rotate wrt the minimizer
    nS2 <- t(Rz(res[3]) %*% Ry(res[2]) %*% Rx(res[1]) %*% t(S2))
    # paste results to current matrices
    X1 <- c(X1[1:(n1 - noverlap)], X2[-(1:noverlap)])
    S1 <- rbind(S1[1:(fixed_point - 1), ], nS2)
    likelihood <- rbind(likelihood, c(k, res))
  }
  
  SS1 <- S1
  X11 <- X1
  
  sr <- block1_end + 2
  block2_end <- length(cutlist)-1   #ceiling(N / block_size) - 1
  
  start <- max(1, cutlist[[sr]][1])  #max(1, (sr - 1) * block_size - noverlap + 1)
  y22 <- y2[start:nrow(y2), start:nrow(y2)] 
  for (k in sr:block2_end) {
    cat("k =", k, "\n")
    if (k == sr) {
      S1 <- result[[k]]$coords ## L1
      X1 <- result[[k]]$X[, 1]
    }
    
    n1 <- nrow(S1)
    S1 <- t(t(S1) - S1[(n1 - noverlap + 1), ])
    S2 <- result[[(k + 1)]]$coords ## L2
    S2 <- S2[(noverlap + 1):nrow(S2), ]
    S2 <- t(t(S2) - S2[1, ])
    n2 <- nrow(S2)
    X2 <- result[[(k + 1)]]$X[, 1]
    
    fixed_point <- n1 - noverlap + 1
    if (fixed_point > N) break
    cat("Fixed Point", n1 - noverlap + 1, "\n")
    y12 <- y22[1:(fixed_point - 1), fixed_point:(fixed_point - 1 + n2)]
    fixed_point <- n1 - noverlap + 1
    llambdax <- exp(outer(X1[1:(n1 - noverlap)], X2[-(1:noverlap)], FUN = "+"))
    n1 <- nrow(S1)
    S1 <- S1[1:(n1 - noverlap), ]
    #min first arugment defualt 50
    res <- minimizer(50, y12,S1,S2,llambdax, threads=CPU)
    nS2 <- t(Rz(res[3]) %*% Ry(res[2]) %*% Rx(res[1]) %*% t(S2))
    X1 <- c(X1[1:(n1 - noverlap)], X2[-(1:noverlap)])
    S1 <- rbind(S1[1:(fixed_point - 1), ], nS2)
    likelihood <- rbind(likelihood, c(k, res))
  }
  SS2 <- S1
  X22 <- X1
  
  
  # final mix
  
  S1 <- SS1
  S2 <- SS2
  S1 <- t(t(S1) - S1[(nrow(S1) - noverlap + 1), ])
  n1 <- nrow(S1)
  S2 <- S2[(noverlap + 1):nrow(S2), ]
  S2 <- t(t(S2) - S2[1, ])
  n2 <- nrow(S2)
  fixed_point <- n1 - noverlap + 1
  cat("Final: Fixed Point", n1 - noverlap + 1, "\n")
  y2_new <- y2[1:(fixed_point - 1), fixed_point:(fixed_point - 1 + n2)]
  
  llambdax <- exp(outer(X11[1:(n1 - noverlap)], X22[-(1:noverlap)], FUN = "+"))
  n1 <- nrow(S1)
  S1 <- S1[1:(n1 - noverlap), ]
  res <- minimizer(50, y2_new,S1,S2,llambdax, threads=CPU)
  nS2 <- t(Rz(res[3]) %*% Ry(res[2]) %*% Rx(res[1]) %*% t(S2))
  FS12 <- rbind(S1[1:(fixed_point - 1), ], nS2)
  likelihood <- rbind(likelihood, c(k+1, res))
  
  return(list(FS12, likelihood))
}

Paste <- function(contact, cutresult, CPU){
    cutlist <- cutresult[[1]] # these are the breakpoints
    result <- cutresult[[2]] # these are the cuts
    noverlap <- cutresult[[3]] # this is noverlap used for the cuts
    
    y2 <- contact 
    N <- nrow(y2)


    sr <- 1
    # initial block
    start <- max(1, (sr - 1) * block_size - noverlap + 1)
    likelihood <- list()
    y11 <- y2[start:nrow(y2), start:nrow(y2)]
    block1_end <- floor(N / block_size / 2) - 1
    k <- 1
    for (k in sr:block1_end) {
      cat("k =", k, "\n")
      if (k == sr) {
        S1 <- result[[k]]$coords ## L1
        X1 <- result[[k]]$X[, 1]
      }
    
      # remove M1 and M2 margins, fix the last of M1 and first of M2 both to be (0,0,0)
      n1 <- nrow(S1)
      S1 <- t(t(S1) - S1[(n1 - noverlap + 1), ])
      S2 <- result[[(k + 1)]]$coords ## L2
      S2 <- S2[(noverlap + 1):nrow(S2), ]
      S2 <- t(t(S2) - S2[1, ])
      n2 <- nrow(S2)
      X2 <- result[[(k + 1)]]$X[, 1]
    
      fixed_point <- n1 - noverlap + 1
      if (fixed_point > N) break
      cat("Fixed Point", n1 - noverlap + 1, "\n")
      y12 <- y11[1:(fixed_point - 1), fixed_point:(fixed_point - 1 + n2)]
      
       ## Attempt APG method first
      apg_result = tryCatch({
          llambdax <- exp(outer(X1[1:(n1 - noverlap)], X2[-(1:noverlap)], FUN = "+"))
          n1 <- nrow(S1)
          # running APG in C++
          S1_trim <- S1[1:(n1 - noverlap), ]
          #min first arugment defualt 50
          res <- minimizer(50, y12,S1_trim,S2,llambdax, threads=CPU)
          
              
          # rotate wrt the minimizer
          nS2 <- t(Rz(res[3]) %*% Ry(res[2]) %*% Rx(res[1]) %*% t(S2))
          # paste results to current matrices
          X1_new <- c(X1[1:(n1 - noverlap)], X2[-(1:noverlap)])
          S1_new <- rbind(S1_trim[1:(fixed_point - 1), ], nS2)
          likelihood_new <- c(k, res)
          
          list(S1 = S1_new, X1 = X1_new, lk = likelihood_new)
          
          }, error = function(e) {
        #   message("APG method failed, falling back to lk_iso3_loglin: ", e$message)
          cat("APG method failed, falling back to lk_iso3_loglin: ", e$message, "\n")
          NULL
        })
      
        # fallback to do the MC
        
        if (is.null(apg_result)) {
          best <- lk_iso3_loglin(S1, S2, noverlap, y12, CPU, X1, X2)
          nS2 <- best$result[[1]]
          X1  <- best$result[[2]]
          S1  <- rbind(S1[1:(fixed_point - 1), ], nS2)
          likelihood[[length(likelihood) + 1]] = best$lkv[which.max(best$lkv[,7]), ]
        } else {
          S1 <- apg_result$S1
          X1 <- apg_result$X1
          
          likelihood[[length(likelihood) + 1]] <- apg_result$lk
        }
      
    }
    
    SS1 <- S1
    X11 <- X1
    
    sr <- block1_end + 2
    block2_end <- ceiling(N / block_size) - 1
    
    start <- max(1, (sr - 1) * block_size - noverlap + 1)
    y22 <- y2[start:nrow(y2), start:nrow(y2)]
    for (k in sr:block2_end) {
      cat("k =", k, "\n")
      if (k == sr) {
        S1 <- result[[k]]$coords ## L1
        X1 <- result[[k]]$X[, 1]
      }
    
      n1 <- nrow(S1)
      S1 <- t(t(S1) - S1[(n1 - noverlap + 1), ])
      S2 <- result[[(k + 1)]]$coords ## L2
      S2 <- S2[(noverlap + 1):nrow(S2), ]
      S2 <- t(t(S2) - S2[1, ])
      n2 <- nrow(S2)
      X2 <- result[[(k + 1)]]$X[, 1]
      
      fixed_point <- n1 - noverlap + 1
      if (fixed_point > N) break
      cat("Fixed Point", n1 - noverlap + 1, "\n")
      y12 <- y22[1:(fixed_point - 1), fixed_point:(fixed_point - 1 + n2)]
      
      ## Attempt APG method first
      apg_result = tryCatch({
          llambdax <- exp(outer(X1[1:(n1 - noverlap)], X2[-(1:noverlap)], FUN = "+"))
          n1 <- nrow(S1)
          S1_trim <- S1[1:(n1 - noverlap), ]
          #min first arugment defualt 50
          res <- minimizer(50, y12,S1_trim,S2,llambdax, threads=CPU)
          
          nS2 <- t(Rz(res[3]) %*% Ry(res[2]) %*% Rx(res[1]) %*% t(S2))
          X1_new <- c(X1[1:(n1 - noverlap)], X2[-(1:noverlap)])
          S1_new <- rbind(S1_trim[1:(fixed_point - 1), ], nS2)
          likelihood_new <- c(k, res)
    
          list(S1 = S1_new, X1 = X1_new, lk = likelihood_new)
          
          }, error = function(e) {
        #   message("APG method failed, falling back to lk_iso3_loglin: ", e$message)
          cat("APG method failed, falling back to lk_iso3_loglin: ", e$message, "\n")
          NULL
        })
      
        # fallback to do the MC
        
        if (is.null(apg_result)) {
          best <- lk_iso3_loglin(S1, S2, noverlap, y12, CPU, X1, X2)
          nS2 <- best$result[[1]]
          X1  <- best$result[[2]]
          S1  <- rbind(S1[1:(fixed_point - 1), ], nS2)
          likelihood[[length(likelihood) + 1]] = best$lkv[which.max(best$lkv[,7]), ]
        } else {
          S1 <- apg_result$S1
          X1 <- apg_result$X1
          
          likelihood[[length(likelihood) + 1]] <- apg_result$lk
        }
      
      
    }
    SS2 <- S1
    X22 <- X1
    
    
    # final mix
    
    S1 <- SS1
    S2 <- SS2
    S1 <- t(t(S1) - S1[(nrow(S1) - noverlap + 1), ])
    n1 <- nrow(S1)
    S2 <- S2[(noverlap + 1):nrow(S2), ]
    S2 <- t(t(S2) - S2[1, ])
    n2 <- nrow(S2)
    fixed_point <- n1 - noverlap + 1
    cat("Final: Fixed Point", n1 - noverlap + 1, "\n")
    y2_new <- y2[1:(fixed_point - 1), fixed_point:(fixed_point - 1 + n2)]
    
    
     ## Attempt APG method first
    apg_result = tryCatch({
      llambdax <- exp(outer(X11[1:(n1 - noverlap)], X22[-(1:noverlap)], FUN = "+"))
      n1 <- nrow(S1)
      # running APG in C++
      S1_trim <- S1[1:(n1 - noverlap), ]
      #min first arugment defualt 50
      res <- minimizer(50, y2_new,S1_trim,S2,llambdax, threads=CPU)
      
      
      # rotate wrt the minimizer
      nS2 <- t(Rz(res[3]) %*% Ry(res[2]) %*% Rx(res[1]) %*% t(S2))
      # paste results to current matrices
      X1_new <- c(X1[1:(n1 - noverlap)], X2[-(1:noverlap)])
      S1_new <- rbind(S1_trim[1:(fixed_point - 1), ], nS2)
      likelihood_new <- c(-1, res)
      
      list(S1 = S1_new, X1 = X1_new, lk = likelihood_new)
      
      }, error = function(e) {
    #   message("APG method failed, falling back to lk_iso3_loglin: ", e$message)
      cat("APG method failed, falling back to lk_iso3_loglin: ", e$message, "\n")
      NULL
    })
    
    # fallback to do the MC
    
    if (is.null(apg_result)) {
      best <- lk_iso3_loglin(S1, S2, noverlap, y2_new, CPU, X11, X22)
      nS2 <- best$result[[1]]
      X1  <- best$result[[2]]
      FS12  <- rbind(S1[1:(fixed_point - 1), ], nS2)
      likelihood[[length(likelihood) + 1]] = best$lkv[which.max(best$lkv[,7]), ]
    } else {
      FS12 <- apg_result$S1
      X1 <- apg_result$X1
      
      likelihood[[length(likelihood) + 1]] <- apg_result$lk
    }
    
    return(list(FS12, likelihood))

}



# run both cut and paste
CutAndPaste <- function(contact, breaks = NULL, block_size = 40, noverlap = 1, CPU){
  
  cuts = Cut(contact=contact, breaks = breaks, block_size = block_size, noverlap = noverlap, CPU=CPU)
  paste = Paste(contact = contact, cutresult=cuts, CPU=CPU)
  
  return(list(cuts, paste))
  
}


Rx_R <- function(S, theta) {
    R <- matrix(0, ncol = 3, nrow = 3)
    R[1, 1] <- 1
    R[2, 2] <- cos(theta)
    R[2, 3] <- -sin(theta)
    R[3, 2] <- sin(theta)
    R[3, 3] <- cos(theta)
    a <- S %*% t(R)
    # a1 <- t(R %*% t(S))
}

Ry_R <- function(S, theta) {
    R <- matrix(0, ncol = 3, nrow = 3)
    R[1, 1] <- cos(theta)
    R[1, 3] <- sin(theta)
    R[2, 2] <- 1
    R[3, 1] <- -sin(theta)
    R[3, 3] <- cos(theta)
    a <- S %*% t(R)
    # a <- t(R %*% t(S))
}

Rz_R <- function(S, theta) {
    R <- matrix(0, ncol = 3, nrow = 3)
    R[1, 1] <- cos(theta)
    R[1, 2] <- -sin(theta)
    R[2, 1] <- sin(theta)
    R[2, 2] <- cos(theta)
    R[3, 3] <- 1
    a <- S %*% t(R)
    # a <- t(R %*% t(S))
}

# fast way to calculate vector distance
#this is sqrt
mat_dist <- function(x, y) {
    sqrt(outer(rowSums(x^2), rowSums(y^2), "+") - tcrossprod(x, 2 * y))
    # outer(rowSums(x^2), rowSums(y^2), "+") - tcrossprod(x, 2 * y)
}

#square ddd and exp log lambda
beta_min <- function(beta, llambdax, y12, ddd) {
    ddd <- ddd^2
    y <- sum(log(ddd) * ddd^beta * exp(llambdax)) - sum(y12 * log(ddd))
    return(y)
}



lk_iso3_loglin=function(S1,S2,noverlap,y12,CPU,X1,X2){
    n1 <- nrow(S1)
    n2 <- nrow(S2)
    fixed_point <- n1 - noverlap + 1

    X <- c(X1[1:(n1 - noverlap)], X2[-(1:noverlap)])
    llambdax <- outer(X1[1:(n1 - noverlap)], X2[-(1:noverlap)], FUN = "+")
   
    cal_lk=function(x){
        ddd=matrix(nrow=(fixed_point-1),ncol=ncol(y12),0)
        
        a1=runif(1,0,2*pi)
        a2=runif(1,0,2*pi)
        a3=runif(1,0,2*pi)
        
        nS2=Rx_R(S2,a1)
        nS2=Ry_R(nS2,a2)
        nS2=Rz_R(nS2,a3)
        
        scale=runif(1)*0.2+0.9
        nS2=scale*nS2
        
        rrefx=runif(1);rrefy=runif(1);rrefz=runif(1)
        refx=ifelse(rrefx<0.5,-1,1)
        refy=ifelse(rrefy<0.5,-1,1)
        refz=ifelse(rrefz<0.5,-1,1)
        nS2[,1]=refx*nS2[,1];nS2[,2]=refy*nS2[,2];nS2[,3]=refz*nS2[,3]



		for(i in 1:(n1-noverlap)){
		    for(j in fixed_point:(fixed_point-1+n2)){
		        jstar=j-fixed_point+1
		        ddd[i,jstar]=sqrt(sum((S1[i,]-nS2[jstar,])^2))
		        #cat("i,j,jstar",i,j,jstar,"\n")
		    }
		}
        
        beta1=runif(1,min=-3,max=0)
        log.lambda=beta1*log(ddd)+llambdax
        lambda=exp(log.lambda)
	
        lv=-lambda+y12*log.lambda
        v=sum(lv)
        pick=c(a1,a2,a3,refx,refy,refz,v,scale,beta1)
        return(pick)
    }
    

    NN=50000
    
    value=mclapply(1:NN,cal_lk,mc.cores = getOption("mc.cores", CPU))
    lkv=do.call(rbind,value)
    pick=lkv[which.max(lkv[,7]),]

    cat(pick,"\n")

    nS2=Rx_R(S2,pick[1])
    nS2=Ry_R(nS2,pick[2])
    nS2=Rz_R(nS2,pick[3])
    nS2[,1]=pick[4]*nS2[,1];nS2[,2]=pick[5]*nS2[,2];nS2[,3]=pick[6]*nS2[,3]
    nS2=pick[8]*nS2
    # list(nS2,X,pick[7])
    
    return(list(result = list(nS2, X, pick[9]), lkv = lkv))
}

