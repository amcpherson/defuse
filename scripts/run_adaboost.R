require(ada)

args <- commandArgs(TRUE)

controls <- read.table(args[1], header=T);
data <- read.table(args[2], header=T);

if (length(rownames(data)) == 0) {
	write.table(data, file=args[3], quote=F, sep="\t", eol="\n", row.names=F)
	q()
}

features <- c(
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

model = ada(controls_features,controls_class)

data_features <- data[,match(features,names(data))]

if (length(rownames(data)) == 1) {
	data_features = rbind(data_features,data_features)
	predict = predict(model,data_features,type="prob")
	data_prob = data.frame(data,probability=predict[1,2])
} else {
	predict = predict(model,data_features,type="prob")
	data_prob = data.frame(data,probability=predict[,2])
}

write.table(data_prob, file=args[3], quote=F, sep="\t", eol="\n", row.names=F)
