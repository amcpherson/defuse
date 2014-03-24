

pe_random <- function()
{
	k=5
	n=500
	u=400
	s=50
	
	set.seed(3)
	
	# random selection of breakpoints B
	B = matrix(round(runif(k*2,0,10000)), nrow=k)
	
	# mixing weights W
	W = runif(k,0,1); 
	W = W / sum(W);
	
	# create n simulated fragments
	Z = rmultinom(n, 1, W)
	
	# generate fragment lengths
	lengths = as.integer(rnorm(n, u, s))
	
	# generate breakpoint offsets
	offsets = as.integer(runif(n, 0, 1) * lengths)
	
	# generate positions
	X1 = colSums(B[,1] * Z) - offsets
	Y1 = colSums(B[,2] * Z) - (lengths - offsets)
	
	# create n simulated fragments
	Z = rmultinom(n, 1, W)
	
	u1 = 200
	s1 = 30
	
	# generate fragment lengths
	lengths = as.integer(rnorm(n, u1, s1))
	
	# generate breakpoint offsets
	offsets = as.integer(runif(n, 0, 1) * lengths)
	
	# generate positions
	X2 = colSums(B[,1] * Z) - offsets
	Y2 = colSums(B[,2] * Z) - (lengths - offsets)
	
	plot(X1,Y1,col="red",xlim=c(min(c(X1,X2)),max(c(X1,X2))),ylim=c(min(c(Y1,Y2)),max(c(Y1,Y2))))
	points(X2,Y2,col="blue")
	
	X1a = X1 + Y1 + u
	Y1a = X1 - Y1
	
	X2a = X2 + Y2 + u1
	Y2a = X2 - Y2
	
	plot(X1a,Y1a,col="red",xlim=c(min(c(X1a,X2a)),max(c(X1a,X2a))),ylim=c(min(c(Y1a,Y2a)),max(c(Y1a,Y2a))))
	points(X2a,Y2a,col="blue")
	
	X = c(X1,X2)
	Y = c(Y1,Y2)
	U = c(rep(u,n),rep(u1,n))
	S = c(rep(s,n),rep(s1,n))
	
	noise = 20
	X = c(X1,X2,runif(noise,0,10000))
	Y = c(Y1,Y2,runif(noise,0,10000))
	U = c(rep(u,n),rep(u1,n),rep(u,noise))
	S = c(rep(s,n),rep(s1,n),rep(u,noise))
	
	sprintf("x = {%s};", paste(X,sep="",collapse=","))
	sprintf("y = {%s};", paste(Y,sep="",collapse=","))
	sprintf("u = {%s};", paste(U,sep="",collapse=","))
	sprintf("s = {%s};", paste(S,sep="",collapse=","))
}

lines=read.table("test1")
for (i in 1:length(lines$V1)) {
  lines(c(lines$V1[i],lines$V3[i]),c(lines$V2[i],lines$V4[i]))
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
		ll = -1e100
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
		
		print(B1)
	
		# Find max likelihood mix weights (M Step)
		W1 = colSums(R1) / length(X)
		
		pe_plot(X,Y,B1)
		
		likelihood = pe_likelihood(X, Y, U, S, B1, W1, l)
		
		# Check for converged clusters, else check for convergence
		if (length(likelihoods) > 0 && abs(likelihood - likelihoods[length(likelihoods)]) < tol) {
			break
		}
		
		if (length(likelihoods) > 1 && likelihood - likelihoods[length(likelihoods)] < -0.000001) {
			print("error")
			return(c())
		}
		
		print(likelihood)
				
		# Output likelihood
		likelihoods = c(likelihoods,likelihood)
		
		# Calculate responsibilities (E Step)
		R1 = pe_resp(X, Y, U, S, B1, W1, l)
	}
	
	print("done")
	
	return(list(likelihood = likelihoods[length(likelihoods)], breaks = B1)) 
}

pe_em <- function(X, Y, U, S, l, mindist) {
	k = 1
	R1 = matrix(1, length(X), 1)
	
	while (T) {
		likelihoods = c()
		while (T) {
			# Find max likelihood breakpoints (M Step)
			B2 = c()
			for (i in 1:k)
			{	
				B2 = c(B2, pe_max_likelihood(X, Y, U, S, R1[,i], l))
			}
			B1 = matrix(B2,k,2,byrow=T)
			
			if (k > 1) {
				for (i1 in 1:(k-1)) {	
					for (i2 in (i1+1):k) {
						if (abs(B1[i1,1]-B1[i2,1]) + abs(B1[i1,2]-B1[i2,2]) < mindist) {
							return(k-1)
						}
					}
				}
			}
			
			print(B1)
			
			# Find max likelihood mix weights (M Step)
			W1 = colSums(R1) / length(X)
			
			pe_plot(X,Y,B1)
			
			likelihood = pe_likelihood(X, Y, U, S, B1, W1, l)
			
			# Check for converged clusters, else check for convergence
			if (length(likelihoods) > 0 && abs(likelihood - likelihoods[length(likelihoods)]) < tol) {
				break
			}
			
			if (length(likelihoods) > 1 && likelihood - likelihoods[length(likelihoods)] < -0.000001) {
				print("error")
				return(c())
			}
			
			print(likelihood)
					
			# Output likelihood
			likelihoods = c(likelihoods,likelihood)
			
			# Calculate responsibilities (E Step)
			R1 = pe_resp(X, Y, U, S, B1, W1, l)
		}
				
		pep = matrix(0, length(X), 1)
		for (i in 1:k) {
			pep = pep + W1[i] * pe_prob(X, Y, U, S, B1[i,], l)
		}
		uncovered = which.min(pep)
		R1 = cbind(R1,matrix(0, length(X), 1))
		k = k + 1
		R1[uncovered,] = 0
		R1[uncovered,k] = 1
		
		# Find max likelihood breakpoints (M Step)
		B2 = c()
		for (i in 1:k)
		{	
			B2 = c(B2, pe_max_likelihood(X, Y, U, S, R1[,i], l))
		}
		B1 = matrix(B2,k,2,byrow=T)
		
		# Find max likelihood mix weights (M Step)
		W1 = c(W1 * (1 - 1/k), 1/k)

		# Calculate responsibilities (E Step)
		R1 = pe_resp(X, Y, U, S, B1, W1, l)
	}
	
	print("done")
	
	return(list(likelihood = likelihoods[length(likelihoods)], breaks = B1)) 
}

pe_em(X, Y, U, S, l, 20)

X=c(-43940146, -43940250, -43941593, -43940269, -43940146, -43940250, -43941570, -43940269, -43940250, -43940269, -43943886, -43943990, -43943465, -43943113, -43943886, -43942981, -43943990, -43942983, -43942981, -43943990, -43942983)
Y=c(43941997, 43940231, 43941641, 43940243, 43943869, 43943971, 43943430, 43943983, 43945838, 43945850, 43941997, 43940231, 43941641, 43941234, 43943869, 43942948, 43943971, 43942928, 43944820, 43945838, 43944800)
U=c(33, 33, 36, 33, 33, 33, 33, 33, 33, 33, 33, 33, 36, 33, 33, 33, 33, 33, 33, 33, 33)
S=c(37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37)

X=c(-49024134, -49024164, -49024050, -49024185, -49024339, -49024167, -49024146, -49024138, -49024135, -49024179, -49024143, -49024170, -49024150, -49024170, -49024163, -49024164, -49024173, -49024144, -49024139, -49024145, -49024134, -49024173, -49024130, -49024134, -49024168, -49025753, -49024567, -49025254, -49025703, -49025275, -49025268, -49024639, -49024584, -49024643, -49024565, -49024823, -49024586, -49025723, -49025282, -49024588, -49025268, -49024637, -49024588, -49024573, -49025710, -49025768, -49025963, -49024609, -49024588, -49025739, -49024824, -49024805, -49025808, -49024637, -49025729, -49024627, -49025270, -49025270, -49025172, -49025190, -49025243, -49024583, -49024829, -49024640, -49025263, -49024570, -49025284, -49024572, -49024592, -49025280, -49024576, -49025709, -49024654, -49024680, -49025172, -49025768, -49024661, -49024669, -49025270, -49025256, -49024659, -49024641, -49025970, -49024584, -49025754, -49025190, -49024645, -49025748, -49024673, -49024672, -49024586, -49024677, -49025725, -49025739, -49024815, -49024674, -49025685, -49025283, -49024571, -49024651, -49025274, -49025734, -49025721, -49025257, -49025280, -49024609, -49025962, -49025276, -49024643, -49025246, -49025267, -49025735, -49025280, -49026078, -49026108, -49026088, -49026102, -49026071)
Y=c(142123802, 142123795, 142123767, 142123857, 142123971, 142123802, 142123803, 142123789, 142123775, 142123853, 142123830, 142123849, 142123852, 142123806, 142123776, 142123807, 142123814, 142123793, 142123801, 142123807, 142123784, 142123802, 142123775, 142123776, 142123810, 142125384, 142124204, 142124858, 142125359, 142124951, 142124871, 142124279, 142124214, 142124275, 142124185, 142124492, 142124241, 142125348, 142124926, 142124218, 142124878, 142124297, 142124225, 142124209, 142125369, 142125369, 142125583, 142124279, 142124192, 142125392, 142124497, 142124514, 142125381, 142124270, 142125361, 142124265, 142124902, 142124902, 142124784, 142124850, 142124881, 142124229, 142124519, 142124290, 142124840, 142124220, 142124877, 142124233, 142124215, 142124875, 142124213, 142125340, 142124295, 142124277, 142124787, 142125359, 142124270, 142124302, 142124876, 142124877, 142124341, 142124316, 142125593, 142124216, 142125384, 142124787, 142124290, 142125390, 142124345, 142124338, 142124266, 142124319, 142125392, 142125377, 142124493, 142124330, 142125337, 142124878, 142124175, 142124294, 142124932, 142125389, 142125366, 142124859, 142124911, 142124291, 142125549, 142124934, 142124320, 142124853, 142124862, 142125327, 142124875, 142125712, 142125719, 142125740, 142125721, 142125707)
U=c(33, 33, 34, 33, 33, 33, 33, 33, 34, 33, 33, 33, 34, 33, 33, 33, 33, 33, 33, 33, 33, 33, 33, 33, 35, 33, 33, 34, 33, 33, 33, 33, 33, 33, 33, 33, 33, 33, 38, 33, 33, 36, 33, 33, 33, 33, 33, 33, 33, 33, 33, 33, 37, 35, 33, 33, 34, 34, 60, 34, 33, 33, 33, 33, 34, 33, 34, 33, 34, 35, 34, 33, 33, 33, 61, 33, 33, 33, 34, 33, 33, 33, 33, 33, 33, 35, 33, 33, 33, 35, 33, 33, 33, 33, 33, 33, 33, 33, 33, 35, 35, 33, 33, 33, 33, 34, 33, 33, 33, 33, 35, 33, 35, 33, 33, 33, 33, 33)
S=c(37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37)

X=c(23975556, 23975573, 23975562, 23975556, 23975571, 23975584, 23975573, 23975570, 23975583, 23975577, 23975562, 23975607, 23975573, 23975533, 23975558, 23975574, 23975564, 23975600, 23975567, 23975559, 23975542, 23975605, 23975533, 23975591, 23975586, 23975584, 23975563, 23975611, 23975563, 23975573, 23975595, 23975618, 23975562, 23975590, 23975545, 23975522, 23975535, 23975563, 23975585, 23975570, 23975573, 23975541, 23975568, 23975584, 23975575, 23975566, 23975584, 23975584, 23975605, 23975606, 23975533, 23975538, 23975533, 23975532, 23975570, 23975567, 23975563, 23975583, 23975570, 23975562, 23975596, 23975577, 23975563, 23975566, 23975607, 23975572, 23975562, 23975533, 23975562, 23975570, 23975541, 23975533, 23975536, 23975533, 23975541, 23975544, 23975556, 23975524, 23975562, 23975566, 23975588, 23975585, 23975571, 23975610, 23975533, 23975563, 23975604, 23975548, 23975606, 23975567, 23975571, 23975564, 23975562, 23975522, 23975560, 23975533, 23975535, 23975532, 23975533, 23975562, 23975533, 23975533, 23975589, 23975574, 23975605, 23975573, 23975573, 23975569, 23975562, 23975584, 23975570, 23975563, 23975564, 23975612, 23975604, 23975589, 23975570, 23975567, 23975574, 23975583, 23975606, 23975570, 23975564, 23975560, 23975545, 23975607, 23975573, 23975566, 23975573, 23975564, 23975617, 23975561, 23975570, 23975558, 23975546, 23975562, 23975540, 23975533, 23975541, 23975607, 23975533, 23975536, 23975532, 23975553, 23975569, 23975591, 23975575, 23975553, 23975556, 23975556, 23975613, 23975573, 23975565, 23975586, 23975612, 23975570, 23975604, 23975571, 23975589, 23975570, 23975533, 23975533, 23975532, 23975570, 23975571, 23975556, 23975580, 23975571, 23975586, 23975570, 23975578, 23975563, 23975573, 23975574, 23975571, 23975544, 23975562, 23975542, 23975563, 23975571, 23975560, 23975568, 23975562, 23975548, 23975541, 23975585, 23975583, 23975573, 23975586, 23975565, 23975611, 23975586, 23975533, 23975533, 23975573, 23975533, 23975601, 23975567, 23975604, 23975543, 23975563, 23975533, 23975584, 23975612, 23975534, 23975583, 23975584, 23975563, 23975574, 23975583, 23975562, 23975553, 23975579, 23975612, 23975589, 23975573, 23975562, 23975606, 23975596, 23975562, 23975606, 23975595, 23975564, 23975605, 23975543, 23975590, 23975562, 23975574, 23975570, 23975532, 23975535)
Y=c(71182647, 71182691, 71182676, 71182641, 71182643, 71182663, 71182640, 71182622, 71182698, 71182629, 71182691, 71182621, 71182652, 71182694, 71182633, 71182640, 71182702, 71182644, 71182670, 71182672, 71182657, 71182649, 71182636, 71182674, 71182627, 71182654, 71182698, 71182590, 71182661, 71182682, 71182644, 71182631, 71182633, 71182634, 71182578, 71182622, 71182637, 71182644, 71182661, 71182683, 71182624, 71182653, 71182647, 71182622, 71182677, 71182693, 71182637, 71182680, 71182614, 71182627, 71182614, 71182580, 71182559, 71182579, 71182691, 71182668, 71182691, 71182642, 71182674, 71182619, 71182661, 71182649, 71182680, 71182674, 71182592, 71182678, 71182627, 71182568, 71182578, 71182575, 71182581, 71182635, 71182603, 71182644, 71182609, 71182624, 71182564, 71182611, 71182673, 71182623, 71182622, 71182641, 71182631, 71182641, 71182698, 71182672, 71182636, 71182659, 71182602, 71182695, 71182697, 71182702, 71182593, 71182663, 71182571, 71182606, 71182593, 71182631, 71182634, 71182571, 71182586, 71182635, 71182629, 71182698, 71182600, 71182652, 71182687, 71182641, 71182626, 71182698, 71182629, 71182624, 71182659, 71182620, 71182630, 71182675, 71182679, 71182621, 71182698, 71182675, 71182655, 71182669, 71182601, 71182611, 71182628, 71182670, 71182692, 71182679, 71182682, 71182612, 71182622, 71182663, 71182678, 71182630, 71182637, 71182649, 71182684, 71182694, 71182644, 71182622, 71182619, 71182610, 71182614, 71182602, 71182634, 71182634, 71182640, 71182677, 71182631, 71182621, 71182643, 71182635, 71182694, 71182652, 71182626, 71182688, 71182660, 71182683, 71182665, 71182652, 71182680, 71182634, 71182635, 71182644, 71182624, 71182650, 71182669, 71182695, 71182634, 71182660, 71182688, 71182683, 71182682, 71182634, 71182659, 71182669, 71182667, 71182620, 71182637, 71182602, 71182663, 71182644, 71182618, 71182691, 71182648, 71182640, 71182643, 71182659, 71182637, 71182698, 71182659, 71182607, 71182637, 71182688, 71182612, 71182620, 71182621, 71182641, 71182603, 71182692, 71182595, 71182614, 71182601, 71182612, 71182574, 71182632, 71182584, 71182624, 71182642, 71182613, 71182628, 71182649, 71182676, 71182660, 71182680, 71182624, 71182695, 71182682, 71182612, 71182677, 71182590, 71182617, 71182691, 71182652, 71182626, 71182666, 71182683, 71182671, 71182626, 71182634, 71182637)
U=c(60, 42, 51, 57, 42, 35, 44, 47, 42, 33, 53, 33, 39, 79, 60, 36, 50, 43, 53, 53, 70, 39, 124, 35, 34, 34, 53, 33, 49, 37, 42, 33, 58, 34, 121, 121, 121, 50, 33, 41, 40, 71, 47, 35, 41, 50, 37, 34, 38, 36, 124, 121, 123, 121, 41, 48, 47, 36, 44, 59, 35, 38, 50, 53, 36, 54, 54, 128, 121, 123, 121, 121, 121, 124, 122, 121, 121, 121, 55, 51, 34, 38, 38, 33, 95, 52, 35, 68, 33, 43, 38, 58, 123, 122, 121, 123, 121, 121, 122, 121, 123, 123, 36, 63, 34, 39, 46, 42, 54, 34, 41, 33, 49, 37, 33, 33, 40, 53, 45, 35, 41, 46, 53, 63, 72, 33, 50, 43, 36, 56, 33, 40, 43, 61, 71, 59, 70, 79, 70, 33, 122, 124, 123, 121, 40, 36, 45, 56, 60, 54, 42, 43, 45, 33, 39, 52, 34, 44, 34, 46, 123, 130, 121, 42, 44, 60, 35, 39, 33, 49, 37, 54, 46, 44, 44, 68, 51, 74, 46, 39, 52, 50, 47, 69, 124, 36, 36, 36, 35, 65, 38, 33, 79, 80, 44, 78, 39, 48, 35, 72, 52, 83, 33, 37, 121, 37, 34, 56, 39, 36, 50, 62, 33, 38, 35, 36, 47, 33, 37, 51, 33, 39, 45, 38, 80, 33, 47, 36, 46, 121, 122)
S=c(37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37)


for (k in 1:min(50,length(X))) {
	estimate = pe_em(X, Y, U, S, k, l)
	print(-2 * estimate$likelihood + k * 2 * log(length(X)))
}




k = 1
R1 = matrix(1, length(X), 1)

while (T) {
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
	
	print(-2 * likelihoods[length(likelihoods)] + k * 2 * log(length(X)))
			
	pep = matrix(0, length(X), 1)
	for (i in 1:k) {
		pep = pep + W1[i] * pe_prob(X, Y, U, S, B1[i,], l)
	}
	uncovered = which.min(pep)
	R1 = cbind(R1,matrix(0, length(X), 1))
	k = k + 1
	R1[uncovered,] = 0
	R1[uncovered,k] = 1
	
	# Find max likelihood breakpoints (M Step)
	B2 = c()
	for (i in 1:k)
	{	
		B2 = c(B2, pe_max_likelihood(X, Y, U, S, R1[,i], l))
	}
	B1 = matrix(B2,k,2,byrow=T)
	
	# Find max likelihood mix weights (M Step)
	W1 = c(W1 * (1 - 1/k), 1/k)

	# Calculate responsibilities (E Step)
	R1 = pe_resp(X, Y, U, S, B1, W1, l)
}


source("/Users/amcphers/Projects/svn/code/dev/comrad/tools/asdf.R")



X=c(48841507, 48841580, 48841569, 48841504)
Y=c(22288060, 22288060, 22288060, 22288060)
U=c(33, 34, 36, 37)
S=c(37, 37, 37, 37)


