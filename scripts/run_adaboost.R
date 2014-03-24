require(kernlab)

args <- commandArgs(TRUE)

controls <- read.table(args[1], header=T);
data <- read.table(args[2], header=T);

features <- c(
"concordant_ratio",
"break_adj_entropy_min",
"cdna_breakseqs_percident",
"genome_breakseqs_percident",
"est_breakseqs_percident",
"splitr_span_pvalue",
"splitr_pos_pvalue",
"splitr_min_pvalue",
"breakpoint_homology",
"span_coverage_min",
"breakseqs_estislands_percident",
"num_splice_variants",
"splice_score",
"max_repeat_proportion",
"mean_map_count"
)

controls_features <- controls[,match(features,names(controls))]
controls_class <- controls$validated == "Y"

data_features <- data[,match(features,names(data))]

sig = sigest(controls_class~., cbind(controls_features,controls_class), frac = 0.5, na.action = na.omit, scaled = TRUE)[2]
model = ksvm(controls_class~., cbind(controls_features,controls_class), type = "C-svc", kernel = "rbfdot", kpar = list(sigma = sig), C = 1, prob.model = TRUE)
predict = predict(model, data_features, type = "probabilities")
data_prob = data.frame(data, probability=predict[,2])

write.table(data_prob, file=args[3], quote=F, sep="\t", eol="\n", row.names=F)

