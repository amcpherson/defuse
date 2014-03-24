args <- commandArgs(TRUE)

cov_stats <- read.table(args[1], header=T);
split_data <- read.table(args[2], header=F, col.names=c("id","sequence","interlen","split_count","split_pos_average","split_min_average"));

no_prediction <- split_data$split_count == 0
split_data$split_count[no_prediction] <- 0

split_pos_pvalue <- 2 * pnorm(-1.0 * abs(split_data$split_pos_average - 0.5)/(sqrt(cov_stats$split_pos_covariance+(1/(12*split_data$split_count)))))
split_min_pvalue <- pnorm((split_data$split_min_average - 0.5)/(sqrt(cov_stats$split_min_covariance+(1/(12*split_data$split_count)))))

split_pos_pvalue[no_prediction] <- 0
split_min_pvalue[no_prediction] <- 0

pvalue_data <- data.frame(id=split_data$id, split_pos_pvalue=split_pos_pvalue, split_min_pvalue=split_min_pvalue);

write.table(pvalue_data, args[3], quote=F, sep="\t", eol="\n", row.names=F, col.names=F)

