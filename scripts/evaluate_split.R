args <- commandArgs(TRUE)

split_pos_cov_stats <- read.table(args[1], header=T);
split_min_cov_stats <- read.table(args[2], header=T);
split_data <- read.table(args[3], header=F, col.names=c("id","sequence","interlen","split_count","split_pos_average","split_min_average"));

no_prediction <- split_data$split_count == 0
split_data$split_count[no_prediction] <- 0

split_pos_pvalue <- 2 * pnorm(-1.0 * abs(split_data$split_pos_average - 0.5)/(sqrt(split_pos_cov_stats$covariance+(1/(12*split_data$split_count)))))
split_min_pvalue <- pnorm((split_data$split_min_average - 0.5)/(sqrt(split_min_cov_stats$covariance+(1/(12*split_data$split_count)))))

split_pos_pvalue[no_prediction] <- 0
split_min_pvalue[no_prediction] <- 0

pvalue_data <- data.frame(id=split_data$id, split_pos_pvalue=split_pos_pvalue, split_min_pvalue=split_min_pvalue);

write.table(pvalue_data, args[4], quote=F, sep="\t", eol="\n", row.names=F, col.names=F)

