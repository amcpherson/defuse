args <- commandArgs(TRUE)

read_stats <- read.table(args[1], header=T);
cov_stats <- read.table(args[2], header=T);
readlength_trim <- as.numeric(args[3]);
sample_mean <- as.numeric(args[4]);
numfrags <- as.numeric(args[5]);

fraglength_mean <- read_stats$fraglength_mean;
fraglength_stddev <- read_stats$fraglength_stddev;
readlength_max <- min(read_stats$readlength_max,readlength_trim);

fraglength_variance <- read_stats$fraglength_stddev^2;

fraglength_mean_adjusted <- fraglength_mean + fraglength_variance / (fraglength_mean - 2 * readlength_max);
fraglength_variance_adjusted <- fraglength_variance - fraglength_variance^2 / (fraglength_mean - 2 * readlength_max)^2;

sample_mean_variance <- fraglength_variance / numfrags + (numfrags - 1) * cov_stats$span_covariance / numfrags;
sample_mean_variance_adjusted <- fraglength_variance_adjusted / numfrags + (numfrags - 1) * cov_stats$span_covariance / numfrags;

all_mean <- c(fraglength_mean, fraglength_mean_adjusted);
all_variance <- c(sample_mean_variance, sample_mean_variance_adjusted);

fraglength_prob <- dnorm((sample_mean - all_mean)/sqrt(all_variance), log=T);
fraglength_pval <- 2 * pnorm(-abs(sample_mean - all_mean)/sqrt(all_variance));

max <- max(fraglength_prob);

if (sample_mean >= fraglength_mean && sample_mean <= fraglength_mean_adjusted)
{
	print(1);
} else {
	print(fraglength_pval[fraglength_prob == max][1]);
}

