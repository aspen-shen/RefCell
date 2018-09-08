#' Find peaks of a 1D distribution
#'
#' This function finds peaks for a 1D density distribution and return positions of these peaks in a list
#'
#' This function is taken from flowDensity package at ()

.find_peaks <- function(x, y = NULL, num_peaks = NULL, adjust = 2, plot = FALSE, ...) {
  x <- as.vector(x)

  if (length(x) < 2) {
    warning("At least 2 observations must be given in 'x' to find peaks.")
    return(NA)
  }

  if (is.null(y)) {
    dens <- density(x, adjust = adjust, ...)
  } else {
    y <- as.vector(y)
    if (length(x) != length(y)) {
      stop("The lengths of 'x' and 'y' must be equal.")
    }
    dens <- list(x = x, y = y)
  }

  # Discrete analogue to a second derivative applied to the KDE. See details.
  second_deriv <- diff(sign(diff(dens$y)))
  which_maxima <- which(second_deriv == -2) + 1

  # The 'density' function can consider observations outside the observed range.
  # In rare cases, this can actually yield peaks outside this range.  We remove
  # any such peaks.
  which_maxima <- which_maxima[findInterval(dens$x[which_maxima], range(x)) == 1]

  # Next, we sort the peaks in descending order based on the density heights.
  which_maxima <- which_maxima[order(dens$y[which_maxima], decreasing = TRUE)]
  
  # Returns the local maxima. If there are none, we return 'NA' instead.
  if (length(which_maxima) > 0) {
    peaks <- dens$x[which_maxima]
    if (is.null(num_peaks) || num_peaks > length(peaks)) {
      num_peaks <- length(peaks)
    }
    peaks <- peaks[seq_len(num_peaks)]
  } else {
    peaks <- NA
  }
  
  peaks <- data.frame(x = peaks, y = dens$y[which_maxima][seq_len(num_peaks)])
  if(plot){
    plot(dens, main = paste("adjust =" ,  adjust))
    points(peaks, ,col = "red")  
  }
  
  peaks  
}

.deriv_density <- function(x, deriv = 1, bandwidth = NULL, adjust = 1,
    num_points = 10000, ...) {
  
  
  if (is.null(bandwidth)) {
    bandwidth <- hpi(x, deriv.order = deriv)
  }
  deriv_x <- drvkde(x = x, drv = deriv, bandwidth = adjust * bandwidth,
      gridsize = num_points, ...)
  list(x = deriv_x$x.grid[[1]], y = deriv_x$est)
}


#' selects typical cells
#'
#' This function selects typical cells and returns selected typical cells as a matrix
#'
#' @param input Input matrix that contains all cells
#' @param number The number of typical cells to select
#' @param index Indices of measurements based on which typical cells are selected, should be a vector of integers
#' @param method Methods to define the center of typical cells, can choose from "mean" and "peak"
#' @return A list contains two elements, the first element is a matrix contain typical cells, and the second element is a vector of index for typical cells
#' @export
typical_cell_selection <- function(input,number,index,method='peak') {
	selected <- matrix(,number,length(index))
	to_use <- input[,index]
	if (method=='mean') {
		median_data <- apply(scale(to_use),2,mean)
	} else if (method=='peak') {
		median_data <- c(1:length(index))
		for (i_c in 1:length(median_data)) {
			peaks <- .find_peaks(scale(to_use[,i_c]),adjust=1)
			if (length(peaks$x>1)) {
				median_data[i_c] <- peaks$x[which(peaks$y==max(peaks$y))]
			} else if (length(peaks$x)==1) {median_data[i_c] <- peaks$x}
			else {stop("no peaks found!!!")}
		}
	}
	pos_dis <- apply(scale(to_use),1,function(x)(sum(abs(x-median_data)^2)))
	rank_dist <- rank(pos_dis,ties.method="first")
	to_select <- which(rank_dist<=number)
	selected <- as.matrix(to_use[to_select,])
	colnames(selected) <- colnames(input)[index]
	return(list(selected,to_select))
}
#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#' appends element to a list object
#'
#' This function appends new element to an existing list object at the end (the original list can be empty)
#'
#' @param lists The original list
#' @param to_add The element to be added to the list
#' @return A new list with the element added at the end
#' @export
appends <- function(lists,to_add) {
	x <- lists
	x[[length(x)+1]] <- to_add
	return(x)
}
#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#' combine all data files within a certain directory
#'
#' This function combines all data files in a given folder
#'
#' @param directory The directory of the folder that stores all the data files to be combined
#' @param pattern An optional regular expression. Only files match this expression with be combined. Default value is NULL. 
#' @return A data frame that contains all data
#' @export
combine_data <- function(directory, pattern = NULL) {
	infiles <- list.files(directory, pattern = pattern)
	outdata <- list()
	for (i_file in 1:length(infiles)) {
		to_add <- read.table(paste0(directory,infiles[i_file]),header=T,colClasses="numeric",sep="\t")
		outdata <- appends(outdata,to_add)
	}
	outdata <- data.table::rbindlist(lapply(outdata,as.data.frame))
	return(outdata)
}
#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#' performs SVM to classify data points in two matrices
#'
#' This function applies SVM in kernlab package to classify two groups of data points stored in two matrices and returns the accuracy 
#' of SVM classification, and the direction of classification boundary
#'
#' @param healthy Matrix of healthy cells to be classified
#' @param disease Matrix of diseased cells to be classified
#' @param feature_index A vector contain index of measurements used in classification, should be a vector of integers, length must be larger than 1
#' @return accuracy Accuracy of the classification
#' @return weightnorm Normalized weights of each feature. Negative weights are higher in disease matrix
#' @return center The center (average) of all data points (healthy and disease combined). Can be used to normalize test datasets
#' @return std The standard deviation of all data points (healthy and disease combined). Can be used to normalize test datasets
#' @return SVM_bn adjusted constant for calculating distance between data point to classification boundary in test data points
#' @export
SVM <- function(healthy,disease,feature_index=c(1:ncol(healthy))) {
	library(kernlab)
	library(ks)
	AvsB <- matrix(1,nrow(healthy),1)
	healthy <- cbind(healthy,AvsB)
	AvsB <- matrix(-1,nrow(disease),1)
	disease <- cbind(disease,AvsB)
	myTable <- rbind(healthy,disease)
	AvsB <- factor(myTable[,ncol(myTable)])
	top_measures <- feature_index

	n_meas <- length(top_measures)
	x <- as.matrix(myTable[,top_measures])
	center <- apply(x,2,mean)
	std <- apply(x,2,sd)
	x <- scale(x)
	myTable2 <- data.frame(x,AvsB)
	myModel <- kernlab::ksvm(AvsB ~ ., data=myTable2,type="C-svc", kernel="vanilladot", C = 10, prob.model=TRUE)

	SVM_coeff <- alpha(myModel)[[1]]
	SVM_coeff2 <- coef(myModel)[[1]]
	SVM_index <- alphaindex(myModel)[[1]]
	weight <- rep(0,each=ncol(x))
	for (i in 1:length(SVM_coeff)) {
		weight <- weight + SVM_coeff2[i]*x[SVM_index[i],]
	}
	norm <- sqrt(sum(weight**2))
	weightnorm <- weight/norm
	label <- sign(SVM_coeff/SVM_coeff2)
	SVM_b <- 0
	for (i in 1:length(SVM_coeff)) {	
		SVM_b <- SVM_b + (sum(weight*x[SVM_index[i],])-label[i])
	}
	SVM_b <- SVM_b/length(SVM_coeff)
	SVM_bn <- SVM_b/norm

	dist <- rep(0,each=nrow(x))
	for (i in 1:nrow(x)) {
		dist[i] <- sum(weightnorm*x[i,])-SVM_bn
	}
	accuracy <- 100*sum(dist*as.numeric(as.character(AvsB))>0)/nrow(x)
	# gap size
	n_pat_class_1 <- sum(AvsB == levels(AvsB)[2])
	gap_size <- min(dist[1:n_pat_class_1])-max(dist[(n_pat_class_1+1):nrow(x)])
	results <- list(accuracy, weightnorm, center, std, SVM_bn)
	names(results) <- c("accuracy","weightnorm", "center", "std", "SVM_bn")
	return(results)
}
#--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#' performs SVM classification and testing
#'
#' This function applies SVM in kernlab package to classify two groups of data points stored in two matrices, 
#' applies its classification boundary to testing dataset and outputs corresponding accuracy.
#'
#' @param controls Training matrix of healthy cells to be classified
#' @param exper Training matrix of diseased cells to be classified
#' @param feature_index A vector contain index of measurements used in classification, should be a vector of integers, length must be larger than 1
#' @param test.control Testing matrix of healthy cells
#' @param test.exp Testing matrix of diseased cells
#' @return accuracy Classificaiton accuracy of the training dataset
#' @return test Classification accuracy of the test dataset
#' @return weight Normalized weights of the features
#' @export
SVM_cv <- function(controls,exper,feature_index,test.control,test.exp) {
	library(kernlab)
	library(ks)
	AvsB <- matrix(1,nrow(controls),1)
	controls <- cbind(controls,AvsB)
	AvsB <- matrix(-1,nrow(exper),1)
	exper <- cbind(exper,AvsB)
	myTable <- rbind(controls,exper)
	AvsB <- factor(myTable[,ncol(myTable)])
	top_measures <- feature_index

	n_meas <- length(top_measures)
	x <- as.matrix(myTable[,top_measures])
	means <- apply(x,2,mean)
	sds <- apply(x,2,sd)
	x <- scale(x)
	# x[is.na(x)] <- 0
	test <- rbind(test.control,test.exp)
	y <- as.matrix(test[,top_measures])
	y <- scale(y,center=means,scale=sds)
	# y[is.na(y)] <- 0
	myTable2 <- data.frame(x,AvsB)
	myModel <- ksvm(AvsB ~ ., data=myTable2,type="C-svc", kernel="vanilladot", C = 10, prob.model=TRUE)

	SVM_coeff <- alpha(myModel)[[1]]
	SVM_coeff2 <- coef(myModel)[[1]]
	SVM_index <- alphaindex(myModel)[[1]]
	weight <- rep(0,each=ncol(x))
	for (i in 1:length(SVM_coeff)) {
		weight <- weight + SVM_coeff2[i]*x[SVM_index[i],]
	}
	norm <- sqrt(sum(weight**2))
	weightnorm <- weight/norm
	weigh <- cbind(as.matrix(weightnorm**2),as.matrix(weightnorm))
	to_weight <- weightnorm**2
	label <- sign(SVM_coeff/SVM_coeff2)
	SVM_b <- 0
	for (i in 1:length(SVM_coeff)) {	
		SVM_b <- SVM_b + (sum(weight*x[SVM_index[i],])-label[i])
	}
	SVM_b <- SVM_b/length(SVM_coeff)
	SVM_bn <- SVM_b/norm

	dist <- rep(0,each=nrow(x))
	for (i in 1:nrow(x)) {
		dist[i] <- sum(weightnorm*x[i,])-SVM_bn
	}
	accuracy <- 100*sum(dist*as.numeric(as.character(AvsB))>0)/nrow(x)
	truth <- c(rep(1,nrow(test.control)),rep(-1,nrow(test.exp)))
	test.dist <- y%*%weigh[,2] - SVM_bn
	CV <- 100*sum(test.dist*truth>0)/nrow(y)
	# gap size
	n_pat_class_1 <- sum(AvsB == levels(AvsB)[2])
	gap_size <- min(dist[1:n_pat_class_1])-max(dist[(n_pat_class_1+1):nrow(x)])
	to_delete <- which(feature_index%in%feature_index[which(to_weight==min(to_weight))])
	results <- list(accuracy, CV, weightnorm)
	names(results) <- c("accuracy","test","weight")
	return(results)
}
#--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#' project additional data points along the direction of classification boundary
#'
#' This function projects additional data points along the direction of classification boundary after normalize it based on the center and std 
#' of classification data (used in SVM), calcualtes the length of projection and percentage of cells along the same direction of classification
#' boundary
#'
#' @param to_be_normed Matrix containing single cells information to be projected
#' @param center Center of typical cells
#' @param std Standard deviations of typical cells
#' @param weight Direction of boundary plane given by SVM(), weight should be formatted as a n by 1 matrix, n is the number of dimension
#' @return projected A list of projected distance between each test data point to the classification boundary
#' @return percent percentage of data points in to_be_normed matrix that have positive distance to the classification boundary
projection <- function(to_be_normed, center, std, weight) {
	scaled <- scale(to_be_normed,center=center,scale=std)
	projected <- to_be_normed%*%weight
	percent <- 100*sum(projected>0)/nrow(to_be_normed)
	output = list(projected,percent)
	names(output) = c("projected","percent")
	return(output)
}
#--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#' batch projection of data in siRNA files to classification boundary
#'
#' This function projects data in multiple data files to the dreiction of classification boundray identified in SVM()
#'
#' @param directory Directory of the folder that contains all data files to be projected
#' @param index A vector of index for features to be used for projection
#' @param center Center of typical cells, outputs of SVM()
#' @param std Standard deviations of typical cells, outputs of SVM()
#' @param weights Weights of each feature (i.e. direction of the classification boundary), outputs of SVM()
#' @return percentages A vector of percentage of cells with positive projection for each data file
#' @return cell.counts A vector of cell numbers in each data file
#' @return names A vector of names for each data file
#' @export
batch_projection <- function(directory,index,center,std,weight) {
	infiles <- list.files(directory)
	output <- c(1:length(infiles))
	cell_count <- c(1:length(infiles))
	name <- c(1:length(infiles))
	for (i in 1:length(infiles)) {
		name[i] <- sub(".txt","",infiles[i])
		indata <- read.table(paste0(directory,infiles[i]),header=T,colClasses="numeric",sep="\t")
		to_use <- indata[,index]
		normed <- projection(to_use, center, std, weight)
		output[i] <- normed[[2]]
		cell_count[i] <- nrow(indata)
	}
	outputs = list(output,cell_count,name)
	names(outputs) = c("percentages","cell.counts","names")
	return(outputs)
}
#--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
FPR <- function(avg,std,threshold,number_of_plates) {
	value <- pt((threshold-avg)*sqrt(number_of_plates)/std,df=number_of_plates-1)
	return(value)
}
#--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#' Identify important siRNA hits
#'
#' This function identifies important siRNA hits as the ones with percentage of healthy-like cells above a threshold
#'
#' @param disease_percent A vector of healthy-like cell percentage in disease control samples
#' @param siRNA_percent A vector of healthy-like cell percentage in siRNA samples
#' @param std_siRNA_percent A vector of standard deviations of each siRNA sample across replicate experiments
#' @param siRNA_number A vector of ratio between the number of cells in each siRNA sample and mean cell numbers in healthy control samples
#' @param siRNA_name A character vector contains the names of each siRNA
#' @param number_of_plates The number of replicate experiments
#' @return A character vector containing the names of siRNA hits
#' @export
hits_identification <- function(disease_percent,siRNA_percent,std_siRNA_percent,siRNA_number,siRNA_name,number_of_plates) {
	threshold <- mean(disease_percent) + 5*sd(disease_percent)
	pre_selected <- which(siRNA_percent>threshold)
	number_selected <- which(siRNA_number>=0.5)
	pre_hit <- intersect(pre_selected,number_selected)
	screened <- rep(0,length(pre_hit))
	index <- 1
	for (i in pre_hit) {
		fdr <- FPR(siRNA_percent[i],std_siRNA_percent[i],threshold,number_of_plates)
		if (fdr<0.2) {screened[index] <- 1}
		index <- index + 1
	}
	selected <- pre_hit[which(screened)]
	hits <- siRNA_name[selected]
	return(hits)
}
#-------------------------------------------------------------- clustering and output siRNAs inside each cluster ----------------------------------------------------------------------------
cluster <- function(avg_siRNA,outdir,siRNA_name,k,master) {
	gene_name <- read.table(master,header=T,sep="\t")
	clusters <- kmeans(avg_siRNA,centers=k,nstart=50)
	output <- paste0(outdir,"siRNA clusters Gene ID.txt")
	output2 <- paste0(outdir,"siRNA clusters symbol.txt")
	for (i_k in 1:k) {
		gene <- c()
		symbol <- c()
		siRNA_name <- rownames(RNA)[which(clusters$cluster==i_k)]
		for (i_RNA in 1:length(siRNA_name)) {
			gene <- c(gene,gene_name$GeneID[grep(paste0(siRNA_name[i_RNA],"$"),as.character(gene_name$GeneSymbol))])
			symbol <- c(symbol,as.character(gene_name$GeneSymbol)[grep(paste0(siRNA_name[i_RNA],"$"),gene_name$GeneSymbol)])
		}
		write(file=output,gene,ncolumns=length(gene),sep="\t",append=T)
		write(file=output2,symbol,ncolumns=length(symbol),sep="\t",append=T)
	}
}
#--------------------------------------------------------------- correlation analysis -------------------------------------------------------------------------------------------------------
pre_correlation <- function(directory, channels, store_dir, centers, stds, weightss) {
	infiles <- list.files(directory)
	multiple <- matrix(,length(infiles),(length(channels)+1))
	num_channel <- length(channels)
	name <- c(1:length(infiles))
	for (i in 1:length(infiles)) {
		indata <- read.table(paste0(directory,infiles[i]),header=T,colClasses="numeric",sep="\t")
		name[i] <- sub(".txt","",infiles[i])
		output <- matrix(,nrow(indata),length(channels))
		for (i_channel in 1:length(channels)) {
			index <- channels[[i_channel]]
			to_use <- indata[,index]
			center <- centers[[i_channel]]
			std <- stds[[i_channel]]
			weight <- weightss[[i_channel]]
			transformed <- projection(to_use, center, std, weight)
			output[,i_channel] <- transformed[[1]]
		}
		colnames(output) <- names(channels)
		rownames(output) <- name
		write.table(output,paste0(store_dir,infiles[i]),append=F,row.names=F,sep="\t") # store transformed data
		to_calculate <- output>0
		sums <- apply(to_calculate,1,sum)
#--------------------------------------------------------------calculate percentage of cells healthy-like in multiple channels --------------------------------------------------------------
		for (j_channel in 1:length(channels)) {
			to_add <- 100*length(sums==(num_channel-1)&(!to_calculate[,j_channel]))/nrow(indata)
			multiple[i,j_channel] <- to_add
		}
		multiple[i,ncol(multiple)] <- 100*length(sums==num_channel)/nrow(indata)
	}
	colnames(multiple) <- c(paste("healthy in all channels except",names(channels)),"healthy in all channels")
	rownames(multiple) <- name
	return(multiple)
}
#--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
correlation_analysis <- function(healthy_dir,disease_dir,siRNA_dir,channels,store_dir,center,std,weight) {
	library(GGally)
	healthy <- list()
	disease <- list()
	siRNA <- list()
	for (i in 1:length(healthy_dir)) {
		h_dir <- healthy_dir[i]
		d_dir <- disease_dir[i]
		s_dir <- siRNA_dir[i]
		centers <- center[[i]]
		stds <- std[[i]]
		weightss <- weight[[i]]
		h_add <- pre_correlation(h_dir, channels, store_dir, centers, stds, weightss)
		healthy <- appends(healthy,h_add)
		d_add <- pre_correlation(d_dir, channels, store_dir, centers, stds, weightss)
		disease <- appends(disease,d_add)
		s_add <- pre_correlation(s_dir, channels, store_dir, centers, stds, weightss)
		siRNA <- appends(siRNA,s_add)
	}
	avg_healthy <- apply(simplify2array(healthy), 1:2, mean)
	avg_disease <- apply(simplify2array(disease), 1:2, mean)
	avg_siRNA <- apply(simplify2array(siRNA), 1:2, mean)
	total_data <- rbind(avg_healthy,avg_disease)
	total_data <- rbind(total_data,avg_siRNA)
	correlation <- apply(total_data,2,function(x) cor(x,total_data[,ncol(total_data)]))
	correlation <- correlation[-ncol(correlation)]
	correlation <- as.matrix(correlation)
	rownames(correlation) <- paste(colnames(total_data)[1:length(channels)],colnames(total_data)[ncol(total_data)],sep="_")
	return(list(total_data,correlation))
}