measurement <- function(data)
{
	comb <- function(n, x) {
		tmp = n*(n-1)
  		return(tmp/2)
	}

	rand_index <- function(data)
	{

		matrix = table(data$V1, data$V2)
		sum_row = apply(matrix, 1, sum)
		sum_col = apply(matrix,2,sum)

		b = 0
		a = 0
		n = comb(nrow(data),2)
		nij = 0

		for(i in 1:nrow(matrix)){
			if(sum_row[[i]] > 1){
				a = a + comb(sum_row[[i]],2)
			}

			for(j in 1:ncol(matrix)){
				if(sum_col[[j]] > 1 && i == 1){
					b = b + comb(sum_col[[j]],2)
				}

				if(matrix[i,j] > 1){
					nij = nij + comb(matrix[i,j],2)
				}
			}
		}

		cat("Rand Index = ",a/n)

		tmp = (nij - (a*b)/n)/( ((a+b)/2) - (a*b)/n)

		cat("\nAdjusted Rand Index = ",tmp)


	}

	purity <- function(classes,clusters)
	{
		matrix = table(classes, clusters)
		tmp = sum(apply(matrix, 2, max)/apply(matrix,2,sum)) / length(table(factor(clusters)))
		cat("\nPurity = ",tmp,"\n")
	}

	NMI <- function(classes,clusters)
	{
		n = nrow(data)
		matrix = table(classes, clusters)
		sum_row = apply(matrix, 1, sum)
		sum_col = apply(matrix, 2, sum)

		MI = 0
		entropy_class = 0
		entropy_cluster = 0
		print(matrix)

		tbl = table(clusters, classes)
		conv_df = as.data.frame.matrix(tbl)
		mutual_information = 0
		for (i in 1:nrow(conv_df)) {
			for (j in 1:ncol(conv_df)) {
				if (conv_df[i, j] > 0) {
					mutual_information = mutual_information + ((conv_df[i, j]/sum(tbl)) * log2( (sum(tbl) / sum(conv_df[i, ]) )*( conv_df[i, j] / sum(conv_df[,j]) ) ) )
				}
			}
		}


		# print(-entropy_class)
		# print(-entropy_cluster)
		# NMI = MI/sqrt(entropy_cluster*entropy_class)
		# NMI = (2*MI)/((-entropy_cluster)+(-entropy_class))

		entr_cluster = sum(apply(conv_df, 1, function(x) -(sum(x)/sum(tbl)) * 
			log2(sum(x)/sum(tbl))))
		entr_class = sum(apply(conv_df, 2, function(x) -(sum(x)/sum(tbl)) * 
			log2(sum(x)/sum(tbl))))
		NMI = (mutual_information/((entr_cluster + entr_class)/2))


		cat("NMI = ",NMI)
	}


	all <- function (true_labels, clusters, method = "adjusted_rand_index",summary_stats = FALSE) 
	{
		if (is.integer(true_labels)) 
			true_labels = as.numeric(true_labels)
		if (is.integer(clusters)) 
			clusters = as.numeric(clusters)
		if (!is.vector(true_labels) || !is.numeric(true_labels)) 
			stop("true_labels should be a numeric vector")
		if (!is.vector(clusters) || !is.numeric(clusters)) 
			stop("clusters should be a numeric vector")
		if (length(true_labels) != length(clusters)) 
			stop("the length of the true_labels vector should equal the length of the clusters vector")
		if (!method %in% c("rand_index", "adjusted_rand_index", "jaccard_index", "fowlkes_mallows_index", "mirkin_metric", "purity", "entropy", "nmi", "var_info")) 
			stop("supported methods are 'rand_index', 'adjusted_rand_index', 'jaccard_index', 'fowlkes_mallows_index', 'mirkin_metric', 'purity', 'entropy', 'nmi', 'var_info'")
		tbl = table(clusters, true_labels)
		# print(tbl)
		conv_df = as.data.frame.matrix(tbl)
		conv_df = as.data.frame.matrix(tbl)
		tp_plus_fp = sum(gmp::asNumeric(gmp::chooseZ(rowSums(conv_df), 2)))
		tp_plus_fn = sum(gmp::asNumeric(gmp::chooseZ(colSums(conv_df), 2)))
		tp = sum(gmp::asNumeric(gmp::chooseZ(as.vector(as.matrix(conv_df)), 2)))
		fp = tp_plus_fp - tp
		fn = tp_plus_fn - tp
		tn = gmp::asNumeric(gmp::chooseZ(sum(as.vector(as.matrix(conv_df))), 2)) - tp - fp - fn
		if (summary_stats || method == "adjusted_rand_index") {
			prod_comb = (tp_plus_fp * tp_plus_fn)/gmp::asNumeric(gmp::chooseZ(length(true_labels), 2))
			mean_comb = (tp_plus_fp + tp_plus_fn)/2
		}
		if (summary_stats || method == "purity") {
			tmp_pur = apply(conv_df, 1, max)
			res_purity = sum(tmp_pur)/length(true_labels)
		}
		if (summary_stats || method == "entropy") {
			# tmp_entropy = sum(apply(conv_df, 2, function(x) entropy_formula(x)))
			# res_entropy = -(1/(sum(tbl) * log2(length(unique(true_labels))))) * tmp_entropy
			res_entropy = 0
		}
		if (summary_stats || method == "nmi" || method == "var_info") {
			mutual_information = 0
			for (i in 1:nrow(conv_df)) {
				for (j in 1:ncol(conv_df)) {
					if (conv_df[i, j] > 0) {
						mutual_information = mutual_information + ((conv_df[i, j]/sum(tbl)) * log2( (sum(tbl) / sum(conv_df[i, ]) )*( conv_df[i, j] / sum(conv_df[,j]) ) ) )
					}
				}
			}
			entr_cluster = sum(apply(conv_df, 1, function(x) -(sum(x)/sum(tbl)) * log2(sum(x)/sum(tbl))))
			entr_class = sum(apply(conv_df, 2, function(x) -(sum(x)/sum(tbl)) * log2(sum(x)/sum(tbl))))
			NMI = (mutual_information/((entr_cluster + entr_class)/2))
			VAR_INFO = (entr_cluster + entr_class) - 2 * mutual_information
		}
		if (summary_stats) {
			prec = tp/(tp + fp)
			rec = tp/(tp + fn)
			cat("", "\n")
			cat("----------------------------------------", "\n")
			cat("purity                         :", round(res_purity, 
				4), "\n")
			cat("entropy                        :", round(res_entropy, 
				4), "\n")
			cat("normalized mutual information  :", round(NMI, 4), 
				"\n")
			cat("variation of information       :", round(VAR_INFO, 
				4), "\n")
			cat("----------------------------------------", "\n")
			cat("specificity                    :", round(tn/(tn + 
				fp), 4), "\n")
			cat("sensitivity                    :", round(tp/(tp + 
				fn), 4), "\n")
			cat("precision                      :", round(prec, 
				4), "\n")
			cat("recall                         :", round(rec, 4), 
				"\n")
			cat("F-measure                      :", round(2 * ((prec * 
				rec)/(prec + rec)), 4), "\n")
			cat("----------------------------------------", "\n")
			cat("accuracy OR rand-index         :", round((tp + 
				tn)/(tp + fp + fn + tn), 4), "\n")
			cat("adjusted-rand-index            :", round((tp - 
				prod_comb)/(mean_comb - prod_comb), 4), "\n")
			cat("jaccard-index                  :", round(tp/(tp + 
				fp + fn), 4), "\n")
			cat("fowlkes-mallows-index          :", round(sqrt((tp/((tp + 
				fp))) * (tp/(tp + fn))), 4), "\n")
			cat("mirkin-metric                  :", round(2 * (fp + 
				fn), 4), "\n")
			cat("----------------------------------------", "\n")
		}
		
	}




	# rand_index(data)
	# purity(data$V1,data$V2)
	# NMI(data$V1,data$V2)
	all(as.numeric(data$V1),as.numeric(data$V2),"nmi",TRUE)
}