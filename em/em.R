

pe_random <- function()
{
	k=3
	n=200
	u=137
	s=21
	
	# random selection of breakpoints B
	B = matrix(round(runif(k*2,0,1000)), nrow=k)
	
	# mixing weights W
	while (T)
	{
		W = runif(k,0,1); 
		W = W / sum(W);
		if (all(W > 1/(4*k)))
		{
			break;
		}
	}
	
	# create n simulated fragments
	Z = rmultinom(n, 1, W)
	
	# generate fragment lengths
	lengths = as.integer(rnorm(n, u, s))
	
	# generate breakpoint offsets
	offsets = as.integer(runif(n, 0, 1) * lengths)
	
	# generate positions
	X = colSums(B[,1] * Z) - offsets
	Y = colSums(B[,2] * Z) - (lengths - offsets)
	plot(X,Y)
}

tol=0.0001
l=0.1

max0 <- function(X)
{
	X[X<0] = 0
	return(X)
}

pe_prob <- function(X, Y, U, S, b, l)
{
	return(exp(-0.5 * ((b[1]+b[2] - X - Y - U) / S)^2 - l * max0(X - b[1]) - l * max0(Y - b[2])))
}

pe_resp_log_likelihood <- function(X, Y, U, S, r, b, l)
{
	return(sum(r * (-0.5 * ((b[1]+b[2] - X - Y - U) / S)^2 - l * max0(X - b[1]) - l * max0(Y - b[2]))))
}

pe_exp <- function(X, Y, U, S, b, l)
{
	return(-1/2 * ((b[1]+b[2]-X-Y-U) / S)^2 - l*max0(X-b[1]) - l*max0(Y-b[2]) )
}

pe_resp <- function(X, Y, U, S, B, W, l)
{
	exponents = matrix(0, length(X), dim(B)[1])
	for (i in 1:dim(B)[1])
	{
		exponents[,i] = -0.5 * ((B[i,1]+B[i,2]-X-Y-U) / S)^2 - l*max0(X-B[i,1]) - l*max0(Y-B[i,2])
	}
	maxexp = apply(exponents, 1, function(x) max(x) )
	resp = matrix(0, length(X), dim(B)[1])
	for (i in 1:dim(B)[1])
	{
		resp[,i] = W[i] * exp(exponents[,i] - maxexp) / rowSums(rep(W,each=length(X)) * exp(exponents - maxexp))
	}
	return(resp)
}

pe_likelihood <- function(X, Y, U, S, B, W, l)
{
	ll = 0
	ls = c()
	for (n in 1:length(X))
	{
		sum = 0
		for (k in 1:dim(B)[1])
		{
			sum = sum + W[k] * pe_prob(X[n], Y[n], U[n], S[n], B[k,], l)
		}
		ll = ll + log(sum)
		ls = c(ls,sum)
	}
	
	if (is.infinite(ll))
	{
		ll = -1e-100
	}
	
	return(ll)
}

pe_max_likelihood <- function(X, Y, U, S, r, l)
{
	Xord = cbind(X,r)[order(X,decreasing=T),]
	Yord = cbind(Y,r)[order(Y,decreasing=T),]
	
	Xsum = cbind(Xord,cumsum(Xord[,2]));
	Ysum = cbind(Yord,cumsum(Yord[,2]));
	
	i = 1
	j = 1
	Xcorner = c(Xsum[1,1])
	Ycorner = c(Ysum[1,1])
	Sum = c(0)
	while (i <= dim(Xord)[1] && j <= dim(Xord)[1])
	{
		while (i < dim(Xord)[1] && Xsum[i,1] == Xsum[i+1,1])
		{
			i = i+1
		}
		
		while (j < dim(Yord)[1] && Ysum[j,1] == Ysum[j+1,1])
		{
			j = j+1
		}
		
		if (abs(Xsum[i,3] - Ysum[j,3]) < 1e-150)
		{
			SAVG = 0.5*(Xsum[i,3]+Ysum[j,3]);

			Xcorner = c(Xcorner, Xsum[i,1])
			Ycorner = c(Ycorner, Ysum[j,1])
			Sum = c(Sum, SAVG)
			
			if (i < dim(Xord)[1] && j < dim(Yord)[1]) {
				Xcorner = c(Xcorner, Xsum[i+1,1])
				Ycorner = c(Ycorner, Ysum[j+1,1])
				Sum = c(Sum, SAVG)
			}
			
			i = i+1;
			j = j+1;
		} 
		else if (Xsum[i,3] < Ysum[j,3]) {
			Xcorner = c(Xcorner, Xsum[i,1])
			Ycorner = c(Ycorner, Ysum[j,1])
			Sum = c(Sum, Xsum[i,3])
			if (i < dim(Xord)[1]) {
				Xcorner = c(Xcorner, Xsum[i+1,1])
				Ycorner = c(Ycorner, Ysum[j,1])
				Sum = c(Sum, Xsum[i,3])
			}
			i = i + 1
		} else {
			Xcorner = c(Xcorner, Xsum[i,1])
			Ycorner = c(Ycorner, Ysum[j,1])
			Sum = c(Sum, Ysum[j,3])
			if (j < dim(Xord)[1]) {
				Xcorner = c(Xcorner, Xsum[i,1])
				Ycorner = c(Ycorner, Ysum[j+1,1])
				Sum = c(Sum, Ysum[j,3])
			}
			j = j + 1
		}
	}
	
	rxyu = sum(r * (X + Y + U) / S^2)
	nks = sum(r / S^2)
	
	partial_ab = rxyu - nks * (Xcorner + Ycorner) + l*Sum
	
	minindex = which(partial_ab > 0)[1]
	
	if (is.na(minindex)) {
		return(c(0,0))
	}
	
	aplusb = (rxyu + l * Sum[minindex]) / nks;
	
	if (minindex == 1) {
		min_a = Xcorner[minindex]
		max_a = aplusb - Ycorner[minindex]
		a = 0.5 * (min_a + max_a)
		b = aplusb - a
	} else if (Sum[minindex] != Sum[minindex-1]) {
		a = Xcorner[minindex]
		b = Ycorner[minindex]
	} else {
		min_a = max(Xcorner[minindex], aplusb - Ycorner[minindex-1])
		max_a = min(Xcorner[minindex-1], aplusb - Ycorner[minindex])
		a = 0.5 * (min_a + max_a)
		b = aplusb - a
	}
	
#	if (minindex > 1 && Sum[minindex] != Sum[minindex-1]) {
#		print("error")
#		return(c())
#	}
	
	
#	pac = sum(r * (X + Y)) / s^2 - sum(r) * (a + b - u) / s^2 + l*sum(r[X>a])
#	pbc = sum(r * (X + Y)) / s^2 - sum(r) * (a + b - u) / s^2 + l*sum(r[Y>b])
	
#	print(pac)
#	print(pbc)
	
	return(c(a,b))
}

pe_plot <- function(X, Y, B) {
	plot(X,Y,xlim=c(min(X)-100,max(X)+100),ylim=c(min(Y)-100,max(Y)+100))
	for (i in 1:dim(B)[1])
	{
		abline(v=B[i,1],col=i)
		abline(h=B[i,2],col=i)
		points(B[i,1],B[i,2],pch=2,col=i)
	}
}

pe_plot_cycle <- function(X, Y, B) {
	for (i in 1:dim(B)[1])
	{
		plot(X,Y,xlim=c(min(X)-100,max(X)+100),ylim=c(min(Y)-100,max(Y)+100))
		abline(v=B[i,1],col=i)
		abline(h=B[i,2],col=i)
		Sys.sleep(1)
	}
	
	pe_plot(X, Y, B)
}

pe_em <- function(X, Y, U, S, k, l) {
	if (k == 1) {
		R0 = matrix(1, length(X), 1)
	} else if (k == length(X)) {
		R0 = diag(k)
	} else {
		# Use KKZ to select k points
		D = cbind(X,Y)
		centers = matrix(D[which.max(X*Y),],1);
		for (i in 2:k) {
			sumsqs = c()
			for (j in 1:dim(centers)[1]) {
				diff = D - matrix(rep(centers[j,],each=length(X)),length(X))
				sumsq = rowSums(diff * diff)
				sumsqs = cbind(sumsqs,sumsq)
			}
			minOfRows = apply(sumsqs, 1, function(x) min(x))
			centers = rbind(centers,D[which.max(minOfRows),]);
		}
	
		# Guess breaks with kmeans
		km = kmeans(cbind(X,Y), centers, k)
		
		R0 = matrix(0, length(X), k)
		
		for (i in 1:length(km$cluster)) {
			R0[i,km$cluster[i]] = 1
		}
	}
	
	R1 = R0
	
	likelihoods = c()
	while (T) {
		# Find max likelihood breakpoints (M Step)
		B2 = c()
		for (i in 1:k)
		{	
			B2 = c(B2, pe_max_likelihood(X, Y, U, S, R1[,i], l))
		}
		B1 = matrix(B2,k,2,byrow=T)
	
		# Find max likelihood mix weights (M Step)
		W1 = colSums(R1) / length(X)
		
		pe_plot(X,Y,B1)
		
		likelihood = pe_likelihood(X, Y, U, S, B1, W1, l)
		
		print(likelihood)
		
		# Check for converged clusters, else check for convergence
		if (length(likelihoods) > 0 && abs(likelihood - likelihoods[length(likelihoods)]) < tol) {
			break
		}
		
		if (length(likelihoods) > 1 && likelihood - likelihoods[length(likelihoods)] < -0.000001) {
			print("error")
			return(c())
		}
				
		# Output likelihood
		likelihoods = c(likelihoods,likelihood)
		
		# Calculate responsibilities (E Step)
		R1 = pe_resp(X, Y, U, S, B1, W1, l)
	}
	
	return(list(likelihood = likelihoods[length(likelihoods)], breaks = B1)) 
}



for (k in 1:min(30,length(X))) {
	estimate = pe_em(X, Y, U, S, k, l)
	print(-2 * estimate$likelihood + k * 2 * log(length(X)))
}




X=c(91638652, 91649816, 91646039, 91650117, 91660697, 91665466, 91662395, 91665466, 91662395, 91665462, 91662391, 91657659, 91659367, 91667207, 91671546, 91672568, 91673924, 91677149, 91661316, 91670866, 91682072, 91670186, 91654846, 91656254, 91665467, 91662396)
Y=c(48766941, 48766922, 48766934, 48766939, 48768520, 48766935, 48766935, 48761786, 48761786, 48766940, 48766940, 48766922, 48766922, 48766922, 48772958, 48772958, 48772958, 48772958, 48772958, 48772958, 48772958, 48772958, 48772958, 48766939, 48761783, 48761783)
U=c(33, 38, 40, 33, 34, 42, 42, 33, 33, 33, 33, 38, 38, 38, 36, 36, 36, 36, 36, 36, 36, 36, 36, 33, 37, 37)
S=c(37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37)

