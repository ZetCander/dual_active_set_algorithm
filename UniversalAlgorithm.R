reductDispersion <- function(y) {
	max <- max(y)
	mult <- 1

	if (max > 1e+04) {
		mult <- 2e+04 / max / length(y)
		y <- mult * y
	}

	return(list(y = y, mult = mult))
}

reverseConvert <- function(y, mult) {
		y <- y / mult

		return(y)
}

changePools <- function(finite_diff, k, pools = vector("numeric")) {

	n <- length(finite_diff)

	new_pools <- unique(c(pools, (1 : n)[finite_diff < 0]))
	new_pools <- sort(new_pools)

	return(new_pools)
}

finiteDiff <- function(z, k) {

	n <- length(z)
	
	finite_diff <- z[-c(1 : k)]
	for (i in 1 : (k - 1)) {
		diff <- k - i
		finite_diff <- finite_diff + ((-1) ^ i) * choose(k, i) * z[-c(c(1 : diff), c((n - i + 1) : n))]
	}
	finite_diff <- finite_diff + ((-1) ^ k) * z[-c((n - k + 1) : n)]

	return(finite_diff)
}

universalAlgorithm <- function(y, k, iter = 100, digits = 10, isNeedOuts = FALSE, isNeedPools = FALSE) {
  print(Sys.time())

	# length of input vector
	n <- length(y)

	# vector of intermediate values
	z <- vector("numeric", n)

	# vector of pools
	pools <- vector("numeric")
  
	temp_val <- reductDispersion(y)
	# first approximation
	z <- temp_val$y
	mult <- temp_val$mult
  
	# indicator of monotonicity
	conv_ind <- FALSE

	# iteration
	i <- 1

	if (isNeedPools && file.exists("pools.txt")) {
		file.remove("pools.txt")
	}
	if (isNeedOuts && file.exists("outs.txt")) {
		file.remove("outs.txt")
	}

	#for first iteration
	finite_diff <- finiteDiff(z, k)
	pools <- changePools(finite_diff, k)

	if (length(pools) > 0) {

		while(!conv_ind && i < iter) {

			if (isNeedPools) {
				cat(paste("pools_", i, ": ", sep=""), file="pools.txt", append = TRUE)
				cat(paste(pools, sep=";"), file="pools.txt", append = TRUE) 
				cat("\n\n", file="pools.txt", append = TRUE)
			}
			if (isNeedOuts) {
				cat(paste("outs_", i, ": ", sep=""), file="z.txt", append = TRUE)
				cat(paste(outs, sep=";"), file="z.txt", append = TRUE) 
				cat("\n\n", file="z.txt", append = TRUE)
			}

			z <- quadProg(z, k, pools)
			finite_diff <- finiteDiff(z, k)

			if (all(round(finite_diff, 10) >= 0)) {
				print(Sys.time())
				conv_ind <- TRUE
				break;
			}
			else {
				pools <- changePools(finite_diff, k, pools)
				i <- i + 1
			}
		}
	}
	else {
		conv_ind <- TRUE
	}

	print(Sys.time())
	
	if (i == iter && !conv_ind) {
		cat(paste('Message: Limit of iterations has been reached', "\n"))
	}
	else {
		cat(paste('Message: The solution was found on iteration ', i, "\n"))
	}

	if (mult != 1) {
		z <- reverseConvert(z, mult)
	}

	err <- sum((y - z)^2) / n 

	return(list(out = z, iter = i, err = err))
}

quadProg <- function(z, k, pools) {
	# We need it to solve quadratic optimization problem
	require(quadprog)
  
	n <- length(z)
  
	binom_coeffs <- sapply(k : 0, function(i) ((-1) ^ i) * choose(k, i))

	dmat <- diag(2, n, n)
	dvec <- 2 * z
	bvec <- rep(0, n)
  
	amat <- matrix(0, n, n)
	
	prev <- pools[1]
	for(i in pools) {
		amat[i : (i + k), i] <- binom_coeffs

		if ((i - prev) > 1 && (i - prev) < k) {
			for(j in (prev + 1) : (i - 1)) {
				amat[j : (j + k), j] <- binom_coeffs
			}
		}
		prev <- i
	}
  
	solution <- solve.QP(Dmat = dmat, dvec = dvec, Amat = amat, bvec = bvec, meq = 0)$solution

	return(solution)
}