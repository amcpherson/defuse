
k=3
n=200
u=137
s=21
l=0.1

minW=0.0001
tol=0.0001

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


max0 <- function(X)
{
	X[X<0] = 0
	return(X)
}

pe_prob <- function(X, Y, b, u, s, l)
{
	return(exp(-0.5 * ((b[1]+b[2] - X - Y - u) / s)^2 - l * max0(X - b[1]) - l * max0(Y - b[2])))
}

pe_resp_log_likelihood <- function(X, Y, r, b, u, s, l)
{
	return(sum(r * (-0.5 * ((b[1]+b[2] - X - Y - u) / s)^2 - l * max0(X - b[1]) - l * max0(Y - b[2]))))
}

pe_exp <- function(X, Y, b, u, s, l)
{
	return(-1/2 * ((b[1]+b[2]-X-Y-u) / s)^2 - l*max0(X-b[1]) - l*max0(Y-b[2]) )
}

pe_resp <- function(X, Y, B, W, u, s, l)
{
	exponents = matrix(0, length(X), dim(B)[1])
	for (i in 1:dim(B)[1])
	{
		exponents[,i] = -0.5 * ((B[i,1]+B[i,2]-X-Y-u) / s)^2 - l*max0(X-B[i,1]) - l*max0(Y-B[i,2])
	}
	maxexp = apply(exponents, 1, function(x) max(x) )
	resp = matrix(0, length(X), dim(B)[1])
	for (i in 1:dim(B)[1])
	{
		resp[,i] = W[i] * exp(exponents[,i] - maxexp) / rowSums(rep(W,each=length(X)) * exp(exponents - maxexp))
	}
	return(resp)
}

pe_likelihood <- function(X, Y, B, W, u, s, l)
{
	ll = 0
	ls = c()
	for (n in 1:length(X))
	{
		sum = 0
		for (k in 1:dim(B)[1])
		{
			sum = sum + W[k] * pe_prob(X[n], Y[n], B[k,], u, s, l)
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

pe_max_likelihood <- function(X, Y, r, u, s, l)
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
	
	partial_ab = sum(r * (X + Y)) / s^2 - sum(r) * (Xcorner + Ycorner - u) / s^2 + l*Sum
		
	minindex = which(partial_ab > 0)[1]
	
	if (is.na(minindex)) {
		return(c(0,0))
	}
	
	aplusb = sum(r * (X + Y))/sum(r) + s^2*l*Sum[minindex]/sum(r) + u	
	
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

pe_em <- function(X, Y, k, u, s, l) {
	if (k == 1) {
		B0 = cbind(X[1],Y[1])
	} else if (k == length(X)) {
		B0 = cbind(X,Y)
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
		B0 = km$centers + u/2
	}
	
	# Guess even mix weights
	W0 = matrix(1, 1, dim(B0)[1]) / dim(B0)[1]
	
	# Assign initial parameters
	B1 = B0
	W1 = W0
	
	pe_plot(X,Y,B1)

	likelihoods = c()
	while (T) {
		likelihood = pe_likelihood(X, Y, B1, W1, u, s, l)
		
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
		R1 = pe_resp(X, Y, B1, W1, u, s, l)
	
		# Find max likelihood breakpoints (M Step)
		B2 = c()
		for (i in 1:dim(B1)[1])
		{	
			B2 = c(B2, pe_max_likelihood(X, Y, R1[,i], u, s, l))
		}
		B1 = matrix(B2,k,2,byrow=T)
	
		# Find max likelihood mix weights (M Step)
		W1 = colSums(R1) / n
		
		pe_plot(X,Y,B1)
	}
	
	print(pe_prob(X, Y, B1[1,], u, s, l))
	
	return(list(likelihood = likelihoods[length(likelihoods)], breaks = B1)) 
}



for (k in 1:min(10,length(X))) {
	estimate = pe_em(X, Y, k, u, s, l)
	print(-2 * estimate$likelihood + k * 2 * log(length(X)))
}




# Overestimate k
k0 = 3

# Use KKZ to select k points
D = cbind(X,Y)
centers = matrix(D[which.max(X*Y),],1);
for (i in 2:k0) {
	sumsqs = c()
	for (j in 1:dim(centers)[1]) {
		diff = D - matrix(rep(centers[j,],each=length(X)),length(X))
		sumsq = rowSums(diff * diff)
		sumsqs = cbind(sumsqs,sumsq)
	}
	minOfRows = apply(sumsqs, 1, function(x) min(x) )
	centers = rbind(centers,D[which.max(minOfRows),]);
}

# Guess breaks with kmeans
km = kmeans(cbind(X,Y), centers, k0)
B0 = km$centers + u/2

# Guess even mix weights
W0 = matrix(1, 1, dim(B0)[1]) / dim(B0)[1]

# Assign as initial parameters
B1 = B0
W1 = W0
k1 = k0

# Plot initial parameters
pe_plot(X,Y,B1)

Sys.sleep(1)

likelihoods = c()
while (T) {
	likelihood = pe_likelihood(X, Y, B1, W1, u, s, l)

	# Check for converged clusters, else check for convergence
 	if (length(likelihoods) > 0 && abs(likelihood - likelihoods[length(likelihoods)]) < tol) {
		break
	}
	
	# Output likelihood
	likelihoods = c(likelihoods,likelihood)
	
	# Calculate responsibilities (E Step)
	R1 = pe_resp(X, Y, B1, W1, u, s, l)

	# Find max likelihood breakpoints (M Step)
	B2 = c()
	for (i in 1:dim(B1)[1])
	{	
		B2 = c(B2, pe_max_likelihood(X, Y, R1[,i], u, s, l))
	}
	B1 = matrix(B2,k1,2,byrow=T)

	# Find max likelihood mix weights (M Step)
	W1 = colSums(R1) / n
	#W1 = max0(colSums(R1) - 1) / sum(max0(colSums(R1) - 1))
	
	# Delete dead clusters
	#B1 = B1[W1 > 0,]
	#W1 = W1[W1 > 0]
	#k1 = length(W1)

	# Plot iteration
	pe_plot(X,Y,B1)
	
	print(W1)

	Sys.sleep(1)
}

comp = 1; plot(X[R1[,comp] > 0.5],Y[R1[,comp] > 0.5])




pe_plot_cycle(X,Y,B1)


# Code for resetting clusters
	if (min(W1) < minW){
		print("fail")
		break
		for (i in which(W1 < minW)) {
			B1[i,] = runif(2,0,1000)
			W1[i] = 1 / length(W1)
		}
		W1 = W1 / sum(W1)
	} else



B1





B1 = array(c(800,300,600,100),c(2,2))
rowSums(R1 * t(Z))


Xcorner[minindex-1]
Ycorner[minindex-1]
Xcorner[minindex]
Ycorner[minindex]
partial_ab[minindex-1]
partial_ab[minindex]


partial_a = sum(r * (X + Y)) / s^2 - sum(r) * (a + b + u) / s^2 - l*sum(Xord[Xord[,1] < a,2])
partial_b = sum(r * (X + Y)) / s^2 - sum(r) * (a + b + u) / s^2 - l*sum(Xord[Yord[,1] < b,2])




sum(r * (X + Y)) / s^2 - sum(r) * (a + b + u) / s^2 - l*sum(Xord[Xord[,1] < a,2])
sum(r * (X + Y)) / s^2 - sum(r) * (a + b + u) / s^2 - l*sum(Xord[Yord[,1] < b,2])


pe_likelihood()



Sys.sleep(1)


pe_print <- function(X, Y) {
	cat(sprintf("int xlist[] = {"))
	for (i in 1:length(X)) {
		cat(sprintf("%d", X[i]))
		if (i != length(X)) {
			cat(sprintf(", ", X[i]))
		}
	}
	cat(sprintf("};\n"))

	cat(sprintf("int ylist[] = {"))
	for (i in 1:length(Y)) {
		cat(sprintf("%d", Y[i]))
		if (i != length(Y)) {
			cat(sprintf(", ", Y[i]))
		}
	}
	cat(sprintf("};\n"))
}


astart=-994002
aend=-993802
bstart=991115
bend=991315
avals=c()
bvals=c()
zvals=c()
for (aval in astart:aend)
{
	for (bval in bstart:bend)
	{
		zval = pe_resp_log_likelihood(X, Y, r, c(aval,bval), u, s, l)
		avals = c(avals, aval)
		bvals = c(bvals, bval)
		zvals = c(zvals, zval)
	}
}
zvals=matrix(zvals,nrow=length(astart:aend))
image(astart:aend,bstart:bend,zvals)


		# Calculate responsibilities (E Step)
		R1 = pe_resp(X, Y, B1, W1, u, s, l)
	
		# Find max likelihood breakpoints (M Step)
		B2 = c()
		for (i in 1:dim(B1)[1])
		{	
			B2 = c(B2, pe_max_likelihood(X, Y, R1[,i], u, s, l))
		}
		B1 = matrix(B2,k,2,byrow=T)
	
		# Find max likelihood mix weights (M Step)
		W1 = colSums(R1) / n
		
		pe_plot(X,Y,B1)
		pe_likelihood(X, Y, B1, W1, u, s, l)
		B1



X=c(40050676, 40050799, 40050799, 40050799, 40050799, 40050799, 40050799, 40050799, 40050799, 40050799, 40050799, 40050799, 40050799, 40050799, 40050799)
Y=c(2245, 2019, 2023, 2027, 2031, 2035, 2039, 2043, 2047, 2051, 2055, 2075, 2079, 2083, 2087)

X=c(-786, -786, -242, -242, -378, -606, -378, -606, -782, -782, -782, -351, -351, -351, -579, -579, -579, -807, -807, -807, -302, -530, -758, -302, -530, -758, -302, -530, -758, -327, -555, -783, -327, -555, -783, -327, -555, -783, -351, -351, -351, -579, -579, -579, -807, -807, -807, -544, -772, -544, -772, -544, -772)
Y=c(457, 685, 367, 595, 268, 268, 496, 496, 209, 437, 665, 239, 467, 695, 239, 467, 695, 239, 467, 695, 192, 192, 192, 420, 420, 420, 648, 648, 648, 242, 242, 242, 470, 470, 470, 698, 698, 698, 239, 467, 695, 239, 467, 695, 239, 467, 695, 206, 206, 434, 434, 662, 662)

<<<<<<< .mine

B1 = matrix(1:k*2,k,2,byrow=T)


B1[,1] =c(-960618, -993902, -972421, -988330, -987217)
B1[,2]=c(947790, 991215.408000000054016709327698, 975330, 985084.883333333302289247512817, 990028)
W1=c(0.10526315789473683626198408092, 0.657894736842105309904127352638, 0.0526315789473684181309920404601, 0.15789473684210525439297612138, 0.0263157894736842090654960202301)

pe_likelihood(X, Y, B1, W1, u, s, l)

# Calculate responsibilities (E Step)
R1 = pe_resp(X, Y, B1, W1, u, s, l)

# Find max likelihood breakpoints (M Step)
B2 = c()
for (i in 1:dim(B1)[1])
{	
	B2 = c(B2, pe_max_likelihood(X, Y, R1[,i], u, s, l))
}
B1 = matrix(B2,k,2,byrow=T)

# Find max likelihood mix weights (M Step)
W1 = colSums(R1) / n

pe_likelihood(X, Y, B1, W1, u, s, l)





R1[,1]=c(0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
R1[,2]=c(1, 0, 0, 1, 1, 1, 1, 1, 1, 0, 1, 1, 0, 0, 0, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 0, 1, 1, 1, 0, 0, 0, 0, 0, 1, 1, 1, 1)
R1[,3]=c(0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0)
R1[,4]=c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 1, 1, 0, 0, 0, 0)
R1[,5]=c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)

# Find max likelihood breakpoints (M Step)
B2 = c()
for (i in 1:dim(B1)[1])
{	
	B2 = c(B2, pe_max_likelihood(X, Y, R1[,i], u, s, l))
}
B1 = matrix(B2,k,2,byrow=T)

# Find max likelihood mix weights (M Step)
W1 = colSums(R1) / n

# Plot iteration
pe_plot(X,Y,B1)

pe_likelihood(X, Y, B1, W1, u, s, l)

# Calculate responsibilities (E Step)
R1 = pe_resp(X, Y, B1, W1, u, s, l)

# Find max likelihood breakpoints (M Step)
B2 = c()
for (i in 1:dim(B1)[1])
{	
	B2 = c(B2, pe_max_likelihood(X, Y, R1[,i], u, s, l))
}
B1 = matrix(B2,k,2,byrow=T)

# Find max likelihood mix weights (M Step)
W1 = colSums(R1) / n

pe_likelihood(X, Y, B1, W1, u, s, l)

# Plot iteration
pe_plot(X,Y,B1)
=======
X=c(-1662600, -1662600, -1622198, -1645763)
Y=c(68971268, 68989300, 68956194, 68955447)

X=c(-1804, -1804, -1811)
Y=c(159480110, 159480110, 159424643)

X=c(-5706, -2261, -2838, -3019)
Y=c(66180773, 66106373, 66181391, 66178043)

X=c(15124780, 15124826, 15108406, 15124811, 15124803, 15124828, 15124804, 15124815, 15105843, 15105843, 15111787, 15124795, 15105936, 15124791, 15105843, 15105843, 15107586, 15105843, 15105843, 15124803, 15124803, 15105843, 15105843, 15124783, 15108363, 15124827, 15124783, 15105927, 15124809, 15105843, 15105843, 15105843, 15105843, 15105843, 15105843, 15105843, 15105843, 15105843, 15105843, 15124811, 15124819, 15111797, 15107558, 15105927, 15124783, 15105843, 15105843, 15105843, 15105843, 15105843, 15105843, 15105843, 15105843, 15105843, 15105843, 15105843, 15105843, 15105843, 15124795, 15124802, 15111809, 15124815, 15111787, 15124804, 15124828, 15107586, 15107605, 15124795, 15124780, 15108363, 15124827, 15124780, 15105843, 15124795, 15124836)
Y=c(153, 117, 410, 130, 132, 110, 136, 129, 926, 983, 429, 153, 718, 157, 926, 983, 340, 926, 983, 132, 132, 929, 986, 151, 475, 103, 151, 719, 117, 920, 977, 1115, 1222, 1474, 1600, 1726, 1852, 1991, 2117, 130, 122, 424, 562, 719, 151, 919, 976, 926, 983, 1121, 1228, 1354, 1480, 1606, 1732, 1858, 1997, 2123, 153, 131, 400, 129, 429, 136, 110, 340, 501, 153, 153, 475, 103, 153, 992, 153, 101)

X=c(-993917, -972421, -960630, -993910, -993918, -994491, -993907, -993915, -993907, -960618, -993906, -993901, -988339, -960618, -987217, -993906, -993918, -993899, -993906, -988317, -993932, -993906, -993910, -993915, -993910, -960618, -993918, -994133, -993925, -991084, -988321, -972421, -991084, -988321, -993920, -993917, -993903, -993906)
Y=c(991037, 975193, 947665, 991035, 991049, 992004, 991035, 991036, 991035, 947653, 991036, 991043, 985018, 947653, 989891, 991036, 991034, 991026, 991036, 985599, 991054, 991036, 991026, 991035, 991035, 947653, 991035, 992071, 991035, 987200, 985012, 975193, 987200, 985012, 991039, 991035, 991039, 991036)

X=c(99743978, 99745364, 99746754, 99748144, 99749537, 99743868, 99745258, 99746644, 99748034, 99749424)
Y=c(102087034, 102087034, 102087034, 102087034, 102087034, 102088230, 102088230, 102088230, 102088230, 102088230)

X=c(-1930, -2062, -2049, -2544, -2104, -2554, -2052, -1928, -2065, -2542, -2050, -2540)
Y=c(-3229, -2424, -2411, -2329, -2339, -1356, -2414, -3227, -2427, -2327, -2412, -2325)

R1 = matrix(1:k*length(X),length(X),k,byrow=T)
R1[,1]=c(0, 0, 0, 0, 0, 0, 0, 0, 0, 1)
R1[,2]=c(1, 0, 0, 0, 0, 1, 0, 0, 0, 0)
R1[,3]=c(0, 0, 0, 0, 0, 0, 0, 1, 0, 0)
R1[,4]=c(0, 0, 0, 1, 0, 0, 0, 0, 1, 0)
R1[,5]=c(0, 1, 0, 0, 0, 0, 1, 0, 0, 0)
R1[,6]=c(0, 0, 0, 0, 1, 0, 0, 0, 0, 0)
R1[,7]=c(0, 0, 1, 0, 0, 0, 0, 0, 0, 0)

		# Find max likelihood breakpoints (M Step)
		B2 = c()
		for (i in 1:dim(B1)[1])
		{	
			B2 = c(B2, pe_max_likelihood(X, Y, R1[,i], u, s, l))
		}
		B1 = matrix(B2,k,2,byrow=T)
	
		# Find max likelihood mix weights (M Step)
		W1 = colSums(R1) / n
		
		pe_plot(X,Y,B1)

		# Calculate responsibilities (E Step)
		R1 = pe_resp(X, Y, B1, W1, u, s, l)
	
		pe_max_likelihood(X, Y, R1[,4], u, s, l)
		r=R1[,4]
		
X=c(99743978, 99745364, 99746754, 99748144, 99749537, 99743868, 99745258, 99746644, 99748034, 99749424)
Y=c(102087034, 102087034, 102087034, 102087034, 102087034, 102088230, 102088230, 102088230, 102088230, 102088230)

X=c(0, 1386, 2776, 4166, 5559, -110, 1280, 2666, 4056, 5446)
Y=c(0,    0,    0,    0,    0, 1196, 1196, 1196, 1196, 1196)


X=c(-188, -235, -242, -138, -243, -242, -238, -243, -209, -229, -198, -240, -183, -256, -188, -190, -190, -197, -235, -190, -210, -138, -196, -198, -176, -226, -235, -243, -157, -208, -238, -154, -235, -195, -188, -265, -213, -226, -205, -195, -212, -190, -212, -246, -181, -198, -212, -195, -217, -196, -205, -208, -247, -169, -225, -234, -198, -198, -145, -212, -156, -148, -143, -212, -202, -196, -221, -144, -222, -301, -197, -182, -208, -188, -197, -160, -208, -158, -200, -231, -191, -238, -139, -240, -243, -247, -208, -192, -188, -213, -235, -189, -220, -324, -366, -195, -226, -212, -189, -229, -202, -198, -195, -256, -197, -194)
Y=c(1128, 1200, 1186, 1075, 1157, 1186, 1186, 1164, 1131, 1135, 1126, 1130, 1130, 1148, 1117, 1143, 1145, 1141, 1170, 1130, 1131, 1075, 1132, 1131, 1124, 1144, 1185, 1177, 1099, 1130, 1181, 1111, 1164, 1146, 1097, 1206, 1141, 1139, 1109, 1121, 1135, 1105, 1143, 1151, 1131, 1153, 1147, 1128, 1123, 1157, 1148, 1138, 1192, 1104, 1130, 1160, 1150, 1118, 1206, 1135, 1065, 1052, 1088, 1110, 1140, 1109, 1193, 1084, 1130, 1145, 1094, 1111, 1155, 1129, 1131, 1118, 1140, 1064, 1109, 1177, 1134, 1192, 1055, 1164, 1186, 1170, 1138, 1124, 1102, 1149, 1143, 1130, 1148, 1190, 1212, 1128, 1167, 1164, 1101, 1150, 1137, 1105, 1139, 1179, 1130, 1130)
u=132.772384937238
s=25.1349843494109
l=0.1

