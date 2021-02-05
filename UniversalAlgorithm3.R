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

	return(new_pools)
}

finiteDiff <- function(z, x_diff) {

	n <- length(z)

	finite_diff <- vector("numeric", n - 3)
	for (i in 1 : (n - 3)) {
		temp <- x_diff[i] + x_diff[i + 1] + x_diff[i + 2]
		finite_diff[i] <- x_diff[i + 1] * (z[i + 3] - z[i])
		finite_diff[i] <- finite_diff[i] + temp * (z[i + 1] - z[i + 2]) 
	}

	return(round(finite_diff))
}

universalAlgorithm <- function(y, iter = 100, x = NULL, digits = 10, isNeedOuts = FALSE, isNeedPools = FALSE) {
  print(Sys.time())

	k <- 3 # Only for 3-monotone
	# length of input vector
	n <- length(y)

	# vector of intermediate values
	z <- vector("numeric", n)
	if (is.null(x)) {
		x <- 1:n
	}
	x_diff <- x[-1] - x[-n]

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

	if (isNeedPools && file.exists("pools1.txt")) {
		file.remove("pools1.txt")
	}
	if (isNeedOuts && file.exists("outs1.txt")) {
		file.remove("outs1.txt")
	}

	#for first iteration
	finite_diff <- finiteDiff(z, x_diff)
	pools <- changePools(finite_diff, k)

	if (length(pools) > 0) {

		while(!conv_ind && i < iter) {

			if (isNeedPools) {
				cat(paste("pools_", i, ": ", sep=""), file="pools1.txt", append = TRUE)
				cat(paste(pools, sep=";"), file="pools1.txt", append = TRUE) 
				cat("\n\n", file="pools1.txt", append = TRUE)
			}
			if (isNeedOuts) {
				cat(paste("outs_", i, ": ", sep=""), file="z1.txt", append = TRUE)
				cat(paste(z, sep=";"), file="z1.txt", append = TRUE) 
				cat("\n\n", file="z1.txt", append = TRUE)
			}

			z <- quadProg(z, k, pools, x_diff)
			finite_diff <- finiteDiff(z, x_diff)

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

	return(list(out = z, iter = i, err = err, pool = pools))
}

quadProg <- function(z, k, pools, x_diff) {
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
	  temp <- x_diff[i] + x_diff[i + 1] + x_diff[i + 2]  
	  coefs <- c(x_diff[i + 1], -temp, temp, -x_diff[i + 1])
		amat[i : (i + k), i] <- coefs

		if ((i - prev) > 1 && (i - prev) < k) {
			for(j in (prev + 1) : (i - 1)) {
				amat[j : (j + k), j] <- coefs
			}
		}
		prev <- i
	}
  
	solution <- solve.QP(Dmat = dmat, dvec = dvec, Amat = amat, bvec = bvec, meq = 0)$solution

	return(solution)
}