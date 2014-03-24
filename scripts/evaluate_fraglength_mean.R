args <- commandArgs(TRUE)

read_stats <- read.table(args[1], header=T);
cov_stats <- read.table(args[2], header=T);
readlength_trim <- as.numeric(args[3]);
spanning_data <- read.table(args[4], header=F, col.names=c("id","mean","count"));

fraglength_mean <- read_stats$fraglength_mean;
fraglength_stddev <- read_stats$fraglength_stddev;
readlength_max <- min(read_stats$readlength_max,readlength_trim);

fraglength_variance <- read_stats$fraglength_stddev^2;
sample_mean_variance <- fraglength_variance / spanning_data$count + (spanning_data$count - 1) * cov_stats$covariance / spanning_data$count;
fraglength_prob <- dnorm((spanning_data$mean - fraglength_mean)/sqrt(sample_mean_variance), log=T);
fraglength_pval <- 2 * pnorm(-abs(spanning_data$mean - fraglength_mean)/sqrt(sample_mean_variance));

fraglength_test <- 1 - pnorm((fraglength_mean - 2 * readlength_max) / read_stats$fraglength_stddev)

if (fraglength_test < 0.05) {
	fraglength_mean_adjusted <- fraglength_mean + fraglength_variance / (fraglength_mean - 2 * readlength_max);
	fraglength_variance_adjusted <- fraglength_variance - fraglength_variance^2 / (fraglength_mean - 2 * readlength_max)^2;
	sample_mean_variance_adjusted <- fraglength_variance_adjusted / spanning_data$count + (spanning_data$count - 1) * cov_stats$covariance / spanning_data$count;
	fraglength_prob_adjusted <- dnorm((spanning_data$mean - fraglength_mean_adjusted)/sqrt(sample_mean_variance_adjusted), log=T);
	fraglength_pval_adjusted <- 2 * pnorm(-abs(spanning_data$mean - fraglength_mean_adjusted)/sqrt(sample_mean_variance_adjusted));

	pvalue = fraglength_pval * (fraglength_prob > fraglength_prob_adjusted) + fraglength_pval_adjusted * (fraglength_prob <= fraglength_prob_adjusted);
	pvalue[spanning_data$mean >= fraglength_mean & spanning_data$mean <= fraglength_mean_adjusted] = 1;
} else {
	pvalue = fraglength_pval
}

pvalue_data <- data.frame(id=spanning_data$id, pvalue=pvalue);

write.table(pvalue_data, args[5], quote=F, sep="\t", eol="\n", row.names=F, col.names=F)
