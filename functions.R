
###############
# Functions
###############

# function printOpenxlsxStyle

printOpenxlsxStyle  <- function (dat, ratios, pvals, wb, tabName = "results", hiCutoff = 1.5, lowCutoff=0.67, pvalCutoff=0.05) 
{
  
  addWorksheet(wb, sheet=tabName)
  
  upReg <- createStyle(fgFill = "violet")
  downReg <- createStyle(fgFill = "forestgreen")
  sigStyle <- createStyle(fgFill = "khaki1")
  
  writeData(wb, tabName, dat, keepNA=FALSE)
  
  for (rat in ratios) {
    up.idx <- which(!is.na(dat[, rat]) & (dat[, rat] > hiCutoff))
    if (length(up.idx) > 0) 
      addStyle(wb, tabName, style=upReg, rows = 1 + up.idx, cols = rat)
    
    down.idx <- which(!is.na(dat[, rat]) & (dat[, rat] < 
                                              lowCutoff))
    if (length(down.idx) > 0) 
      addStyle(wb, tabName, style=downReg, rows = 1 + down.idx, cols = rat)
  }
  
  for (pval in pvals) {
    sig.idx <- which(!is.na(dat[, pval]) & (dat[, pval] < 
                                              pvalCutoff))
    if (length(sig.idx) > 0) 
      addStyle(wb, tabName, style=sigStyle, rows = 1 + sig.idx, cols = pval)
  }
}


# Added on 22/08/18 to fix the overlapped labels
plotErrorBarsLines <- function (v, barSizes, lines, labels = NULL, col = "blue", ylim = c(min(lines), 
                                                                                          max(lines)), ...) 
{
  barSizes[is.na(barSizes)] <- 0
  topBars <- v + 0.5 * barSizes
  bottomBars <- v - 0.5 * barSizes
  N <- length(v)
  if (is.null(labels)) 
    labels <- 1:N
  ylims <- c(min(bottomBars, ylim[1], min(lines)), max(topBars, 
                                                       ylim[2], max(lines)))
  par(pch = 19, xaxt = "n")
  plot(as.numeric(labels), v, ylim = ylims, col = col, type = "b", 
       lwd = 3, ...)
  par(xaxt = "s")
  
  for (i in 1:N) {
    lines(c(i, i), c(topBars[i], bottomBars[i]))
  }
  for (i in 1:ncol(lines)) {
    lines(as.numeric(labels), lines[, i], lwd = 0.5, lty = "dotted", 
          col = "gray")
  }
}

panel.cor <- function(x, y, digits = 2, prefix = "", cex.cor, ...)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- abs(cor(x, y))
  txt <- format(c(r, 0.123456789), digits = digits)[1]
  txt <- paste0(prefix, txt)
  if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
  text(0.5, 0.5, txt, cex = cex.cor * r)
}

make_CVs = function(df, replicate) {
  
  cvs = matrix(NA, nrow=nrow(df), ncol=nlevels(as.factor(replicate) ))
  colnames(cvs) = levels(as.factor(replicate))
  rownames(cvs) = rownames(df)
  
  for(ii in 1:nrow(df))
    cvs[ii,] = aggregate(t(df[ii,]), by=list(replicate), function(x) sd(x)/mean(x))[,2]
  
  return(cvs)
  
}

plotClusterProfile <- function(cluster.data, clustID, group, k=4, ylab="Abundance") {
	
	gp = group
	noClusters <- k

	r.temp <- aggregate(t(cluster.data), by=list(gp=gp), FUN=mean)
	ag.sample <- r.temp[,-1]
	rownames(ag.sample) <- r.temp[,1]
	ag.genes <- aggregate(t(ag.sample), by=list(Cluster=clustID), FUN=mean)
	ag.sd <- aggregate(t(ag.sample), by=list(Cluster=clustID), FUN=sd)
	ag.matrix <- as.matrix(ag.genes[,-1])
	ag.counts <- summary(as.factor(clustID))
	ag.bars <- as.matrix(ag.sd[,-1])
	
	png("ClusterPatterns.png", 2000, 2000, res=300)
	par(bg=gray(.95), fg=gray(0.3), oma= c(5, 2, 2, 1) + 0.1, col.main="black", col.sub="black", col.lab="black", col.axis="black")
	layout(matrix(1:4, ncol=2, byrow=TRUE))
	NSig <- noClusters
	for(i in 1:NSig) {
		cols <- rainbow(4) 
		# cols <- rep("gray", 6)
		gname <-  paste("Cluster", i, "(", ag.counts[i], "proteins )")
		lines <- ag.sample[, clustID==i, drop=FALSE]
		plotErrorBarsLines(ag.matrix[i,], 2*ag.bars[i,], lines, 
			labels=1:ncol(ag.matrix), 
			col=cols[i],  main=gname, # bgcol="gray", split=split,
			ylab=ylab, xlab="",
			ylim=c(min(ag.matrix), max(ag.matrix)))
		axis(1,at=1:ncol(ag.matrix), las=2, labels=colnames(ag.matrix), col="black")
		abline(h=0, lty="dotted")
	}
	
	dev.off()

}




plotDensityNA <- function (data, group = rownames(data), xlab = "Abundance",
		main="Sample Densities", legend=TRUE) 
{
    group <- as.factor(group)
    colours <- rainbow(length(levels(group)))
    col <- colours[group]
    # par(mfrow = c(1, 1))
   

        x <- as.matrix(data)
	  y <- x[!is.na(x)]

        try(d1 <- density(y))
        if (inherits(d1, "try-error")) 
            Error("Failed to generate the density plots")
        ymx <- max(d1$y)
        plot(d1, type = "n", xlab = xlab, ylab = "PDF", main = main)
         #   ylim = c(0, 2 * ymx), yaxp = c(0, 2 * ymx, 5))
        for (i in 1:ncol(x)) {
            try(d1 <- density(x[!is.na(x[,i]), i]))
            if (inherits(d1, "try-error")) 
                Error(paste("Failed to generate the density plot for sample", 
                  i))
            # lines(d1, lty = i, col = col[i])
		lines(d1, col = col[i])

        }
	
	  if (legend) {
	    par(cex=0.5)
        legend("topright", legend = levels(group), col = colours, 
            lty = 1)
         par(cex=1)
		 }

    
}




# Jemma - IRS; needs raw data, not log transformed
IrsNormalisation = function (Mat, Run, method="total") {

	Run = as.factor(Run)
	
	list.exp = list()
	
	for(ii in 1:nlevels(Run)) list.exp[[ii]]=Mat[,Run==levels(Run)[ii] ]
	
	# Calculate the protein sum for each batch
	if(tolower(method) =="total" ) {
		list.rowsum = lapply(list.exp, function(x) apply(x,1, function(y) (sum(na.omit(y)) ) ) )
	} else {
	
		list.rowsum = lapply(list.exp, function(x) apply(x,1, function(y) (median(na.omit(y)) ) ) )
	}
	
	irs = as.data.frame(list.rowsum[[1]])
	if(length(list.rowsum) > 1)
	for(ii in 2:length(list.rowsum)) irs = cbind(irs, as.data.frame(list.rowsum[[ii]]) )
	
	# convert 0 to NA	
	irs[irs==0] = NA
	
	colnames(irs) <- paste0("sum", 1:ncol(irs))
	
	rowsum.average <- apply(irs, 1, function(x) exp(mean(log(na.omit(x)))))
		

	# compute the scaling factor vectors
	
	irs.fac = sweep(irs, 1, rowsum.average, "/")
	
	list.irs.scaled = lapply(1:length(list.exp), function(x) sweep(list.exp[[x]], 1, irs.fac[,x], "/") ) 
		
	# make new data frame with normalized data
	data_irs <- list.irs.scaled[[1]]
	if(length(list.irs.scaled) > 1)
	for(ii in 2:length(list.irs.scaled)) data_irs = cbind(data_irs, list.irs.scaled[[ii]])

	data_irs
}

		
processDataMatrix = function(mg.all, Design, DuplicateFilter=TRUE) {

# remove "scaled columns"
if (length(grep("sn.scaled", colnames(mg.all))) > 0) {
mg.all = mg.all[,-grep("sn.scaled", colnames(mg.all))]
}
colnames(mg.all) = gsub(".rq", "", colnames(mg.all))
colnames(mg.all) = gsub(".sn.sum", "", colnames(mg.all))

labelcols = grep("(12[6789])|(13[01])", colnames(mg.all))
infocols = setdiff(1:ncol(mg.all), labelcols)

cat(labelcols, "\n")

NCOL=max(infocols)
if (nrow(Design) != length(labelcols)) stop("The data file and design file don't match - different numbers of samples.")

file.labels = gsub("(.*)((12[6789])|(13[01]))(.*)", "\\2",  colnames(mg.all)[labelcols])
design.labels = gsub("((12[6789])|(13[01]))(.*)", "\\1", Design$Label)

if (sum(file.labels != design.labels) > 0) stop("Design and file labels do not match")

colnames(mg.all)[infocols] = gsub("Peptides", "Number.of.peptides", colnames(mg.all)[infocols])

# assumption: Protein.Id, Gene.Symbol, Description, Group.ID, others
mg.all = mg.all[,c(1:3, 5:NCOL, 4, (1+NCOL):ncol(mg.all))] 

# Drop duplicates from Group
if (DuplicateFilter) mg.all = mg.all[!duplicated(mg.all$Group.ID),]

return(mg.all)
}

cordist = function (x) 
{
    as.dist((1 - cor(t(x)))/2)
}

PCA = function (data, labelValue, scaleR = FALSE, scaleC = TRUE, k = min(dim(data)) - 
    1) 
{
    if (k > min(dim(data) - 1)) 
        warning("The number of components was too large compared to the data and was adjusted accordingly")
    k <- min(k, min(dim(data)) - 1)
    if (scaleR) {
        row.nrm <- apply(data, 1, sd)
        row.nrm <- pmax(row.nrm, 1e-04)
        data <- sweep(data, 1, row.nrm, FUN = "/")
    }
    result <- try(prcomp(data, retx = TRUE, scale = scaleC), 
        silent = TRUE)
    if (inherits(result, "try-error")) 
        stop("Failed to Calculate Principal Components")
    componentVariances <- result$sdev^2
    componentLoadings <- result$rotation[, 1:k]
    componentScores <- result$x[, 1:k]
    totalVariance <- sum(componentVariances)
    componentVariances <- componentVariances[1:k]
    z <- componentScores
    plot(cloud(z[, 1] ~ z[, 3] + z[, 2], groups = as.factor(labelValue), 
        auto.key = list(points = TRUE, pch = 19, space = "right"), 
        xlab = "PC 3", ylab = "PC 2", zlab = "PC 1", 
        distance = 0.1, main = "Projection in the space of the first 3 princial components"))
    value <- list(componentVariances = componentVariances, componentScores = componentScores, 
        componentLoadings = componentLoadings, summary = summary(result))
    value
}


TMTPrePro = function(mg.all, Design, Comparisons, PvalCutoff, FCCutoff, DuplicateFilter, KeepREF, SampleLoadNorm)	{

dat.para = data.frame(Parameter=c('P value', 'Fold change', 'DuplicateFilter', 'KeepREF', 'SampleLoadNorm', 'Date'), 
			Cutoff=c(PvalCutoff, FCCutoff, DuplicateFilter, KeepREF, SampleLoadNorm, date() ))	



mg.proc = processDataMatrix(mg.all, Design=Design, DuplicateFilter=TRUE)
mg.all = mg.proc


labelcols = grep("(12[6789])|(13[01])", colnames(mg.all))
infocols = setdiff(1:ncol(mg.all), labelcols)
NCOL = max(infocols)


selection.idx = !(Design$Group %in% "REF")

data_raw = mg.all[, -c(1:NCOL)][, selection.idx]	

Design = Design[selection.idx,]

data_raw[data_raw==0] = NA
rownames(data_raw) = mg.all[,"Protein.Id"]

Group = as.factor(Design$Run)


boxplot(log(data_raw))
	
	Cols = rainbow(nlevels(Group))
	
	#################################
	# Sample loading normalistation
	#################################

	# Sample loading median or total normalisation
	if(tolower(SampleLoadNorm) == "total") {
		tot = apply(data_raw, 2, function(x) sum(na.omit(x)))
	} else {
		tot = apply(data_raw, 2, function(x) median(na.omit(x)))
	}

	
	data_sl = sweep(data_raw, 2, tot/max(tot), "/")

	format(round(colSums(na.omit(data_sl)), digits = 0), big.mark = ",")
	
	
	# see what the SL normalized data look like
	png(paste("Boxplot raw and", SampleLoadNorm, "norm.png"), 2000, 2000, res=300)
	layout(matrix(1:4, ncol=2))	
	boxplot(log2(data_raw), col = Cols[Group], 
		notch = TRUE, main = "Plain", las=2, cex.names=.8)
	plotDensityNA(log2(data_raw), group = Group, main="Raw", legend=FALSE)
	
	boxplot(log2(data_sl), col = Cols[Group], 
			notch = TRUE, main = paste(SampleLoadNorm, "normalized data"), las=2, cex.names=.8)
	plotDensityNA(log2(data_sl), group = Cols[Group], main = paste(SampleLoadNorm, "normalization"), legend=FALSE)
	
	dev.off()



data_irs = IrsNormalisation(data_sl, Design$Run, method="median")


	# remove reference ?
	data_irs_full = data_irs
	
	
	# if(KeepREF == FALSE) data_irs = data_irs[,!colnames(data_irs) %in% refInRuns]
	
	# data_irs =a data_sl
	
	Group = as.factor(Design$Group)
				
	
	Replicate = Design$ID	
	
	Cols = rainbow(nlevels(Group))
	
	####################################################
	# Some overall data look metrics
	# -- correlation
	# -- boxplots and density plots
	# -- within group correlation for each level of the group
	# -- PCA
	# -- heatmap
	# -- Anova, followed by clustering of DE proteins
	# -- Diff exp proteins _to the reference_, e.g. Mod/Control, ..., etc, via 1-sample t-test,
	# -- Also combine z-scores
	# -- Venn diagrams of overlap, and also barplots
	####################################################


	# Correlation heatmap
	png("Correlation heatmap IRS.png", 2000, 2000,res=300)
	par(xpd=TRUE)
	heatmap3(cor(log(na.omit(data_irs+.5)), use="pairwise.complete.obs"), distfun=cordist,
		 col=colorRampPalette(c("green","black", "red"))(1024),
		 main="IRS correlation",
		ColSideColors=rainbow(nlevels(Group))[Group], margins=c(20,20))
	legend("topright", fill=Cols[1:nlevels(Group)], legend=levels(Group),  cex=.6)  # xpd=TRUE,
	dev.off()


	png("BoxplotDensity.png", width=3500, height=1700,res=300)
	layout(matrix(1:2, nrow=1))
	par(mar=c(13,4,4,2)+.1)
	plotDensityNA(log(na.omit(data_irs+.5)), group=Group, legend=FALSE, main=paste(SampleLoadNorm, "and IRS normalised"))
	legend('topright', fill=Cols[1:nlevels(Group)], legend=levels(Group), cex=0.5)
	
	# boxplots and density plots
	boxplot(log(data_irs[, order(Group)]+.5), las=2, col=Cols[Group[order(Group)]], 
		main=paste(SampleLoadNorm, "and IRS normalised"),
		cex.axis=0.6)	

	dev.off()

if(FALSE) {

	cat('Begin correlations for all groups\n')
	


	# Correlations for all levels of the group
	for (lev in levels(Group)) {
	
		dd = na.omit(data_irs[,Group == lev, drop=FALSE])
		if(ncol(dd) > 1) {
		png(paste("Cor", lev, ".png", sep=""), 2000, 2000,res=300)
		pairs(log(dd+.5), lower.panel = panel.smooth, upper.panel = panel.cor, main=lev)
		dev.off()
		}
	}	

	cat('End correlations for all groups\n')
	
	}
	
	
	############################################
	# unsupervised analysis: clustering and PCA
	############################################

		cat('Begin unsupervised\n')
		
	png("HeatmapAll.png", 2000, 2000,res=300)
	heatmap3(as.matrix(log(na.omit(data_irs+.5))), distfun=cordist,col=colorRampPalette(c("green","black","red"))(100),
		ColSideColors=rainbow(nlevels(Group))[Group], margins=c(20,20))
	legend("topright", fill=Cols[1:nlevels(Group)], legend=levels(Group),  xpd=TRUE, cex=.6 )
	dev.off()


	#pca.res <- PCA(log(t(na.omit(data_irs+.5))), Group, k=5, scaleC = FALSE)
	pca.res <- PCA(log(t(na.omit(data_irs+.5))), Group, k=5)
	z <- pca.res$componentScores

	ld = pca.res$componentLoadings

	# proportion of variance of the top 3 components
	props = round(100*pca.res$summary$importance[2,1:3], 1)
	
	png("PCA3dPlot.png", 2000, 2000, res=300)
	s3d <- scatterplot3d(z[,1:3], color = Cols[Group], col.axis=gray(0.85),col.grid="lightblue",
		box = T, angle = 26, pch=20 )
	s3d$points3d(z[,1:3], pch=21)
	legend("topleft", fill=Cols[1:nlevels(Group)], legend=levels(Group), cex=.6)
	text(s3d$xyz.convert(3+z[,1], 3+z[,2], z[,3]), labels = colnames(data_irs), cex=0.4)
	dev.off()

	ord.list = list()

	png("PCA2DAll.png", 2000, 2000, res=300)
	layout(matrix(1:4, ncol=2))
	plot(z[,1], z[,2], col=Cols[Group], pch=20, xlab=paste0("PC1(", props[1],"%)"), 
		ylab=paste0("PC2(", props[2], "%)"))
	points(z[,1], z[,2], pch=21, cex=1.1, lwd=1.3)
	text(z[,1], z[,2], colnames(data_irs), pos=3, cex=.5)
	
	plot(z[,1], z[,3], col=Cols[Group], pch=20, xlab=paste0("PC1(", props[1],"%)"), 
		ylab=paste0("PC3(", props[3],"%)"))
	points(z[,1], z[,3], pch=21, cex=1.1, lwd=1.3)
	text(z[,1], z[,3], colnames(data_irs), pos=3, cex=.5)
	
	plot(z[,2], z[,3], col=Cols[Group], pch=20, xlab=paste0("PC2(", props[2], "%)"), 
		ylab=paste0("PC3(", props[3],"%)"))
	points(z[,2], z[,3], pch=21, cex=1.1, lwd=1.3)
	text(z[,2], z[,3], colnames(data_irs), pos=3, cex=.5)
	
	plot(z[,2], z[,3], col=Cols[Group], pch=20, xlab="", ylab="", axes=FALSE, type='n')
		
	legend("center", fill=Cols[1:nlevels(Group)], legend=levels(Group), cex=0.5)	

	dev.off()
	
	png("PCATopLoadings.png", width=2000, height=700, res=300)
	par(oma=c(2,1,1,1))
	layout(matrix(1:3, nrow=1))

	for (ii in 1:3) {
	 ord = order(abs(ld[,ii]), decreasing=TRUE)[1:5]
	 barplot(sort(ld[ord, ii]), las=2, main=paste("Top loadings PC", ii))
	 ord.list[[ii]]=ord
	}
	dev.off()

	png("PCATopLoadingsProteinPatterns.png", width=2500, height=2500, res=300)
	par(mar=c(5,2,3,1))
	layout(matrix(1:15, nrow=3, byrow=T))
	for (ii in 1:3) {
		ord = ord.list[[ii]]
		for (xx in 1:5) {
			boxplot(as.vector(as.matrix(data_irs[match(rownames(ld)[ord[xx]], rownames(data_irs)),])) ~ Group, 
			boxwex=0.5, main=rownames(ld)[ord[xx]], col="gray", las=2)
		}
	}
	dev.off()


	cat('End unsupervised \n')
	############################################
	# End unsupervised analysis
	############################################


	################
	# ANOVA
	################
	
	Group = as.factor(Design$Group)
	# Group = as.factor(Design$Sex)

		cat('Begin anova\n')
	Anova = rep(NA, nrow(data_irs))

	# compute Group means (in log space, geometric)
	data.ag = aggregate(t(data_irs), by=list(Group=Group), 
			FUN=function(v){exp(mean(log(na.omit(v))))} )
	Means = t( data.ag[,-1])
	colnames(Means) = paste("Means",data.ag[,1])
	MaxFC = apply(Means, 1, FUN=function(v){max(v)/min(v)})


	for (i in 1:nrow(data_irs)) {
		
		v= t(data_irs[i,])
		nna.idx = !is.na(v)

		an.res =  try(anova(lm( log(v[nna.idx]+.5) ~ Group[nna.idx, drop=TRUE] ))[1,"Pr(>F)"])

		if (!inherits(an.res, "try-error")) Anova[i] = an.res;

	}



	Anova.adj  = p.adjust(Anova, method = "fdr")
	Anova.idx =  !is.na(MaxFC) & ( MaxFC > FCCutoff ) & !is.na(Anova) & (Anova < PvalCutoff)


	if(nrow(data_irs[Anova.idx,]) > 3) {
  	png("Heatmap - Anova DE.png", 2000, 2000, res=300)
  	hm1 <- heatmap3(as.matrix(na.omit(log(data_irs[Anova.idx,]+.5))),  margins=c(20,20), cexRow=1,
  		 col=colorRampPalette(c("green", "black", "red"))(120), 
  		ColSideColors=Cols[Group]  )
  	legend("topright", fill=Cols[1:nlevels(Group)], legend=levels(Group),
  	 xpd=TRUE,cex=.6 )
  	dev.off()
	}
		cat('End anova\n')

		cat('Begin clustering\n')
	#cluster.data = na.omit(log(data_irs[Anova.idx,]))
	cluster.data = na.omit(log(data_irs_full[Anova.idx,]))
	NotOmitted = setdiff(1:sum(Anova.idx), attr(cluster.data, "na.action"))

	gp = Group


	Cluster = rep(NA, sum(Anova.idx))
	
	res1 <- try(HClust(cluster.data, metric="pearsonCorrelation", method="complete", cutNumber=4))
	
	if(!inherits(res1, "try-error")) {
		clustID <- res1$clustID

		Cluster[NotOmitted] = clustID

		noClusters <- 4

		# scale cluster data for visulisation only	
		# mat.ref = cluster.data[,colnames(cluster.data)%in%refInRuns, drop=FALSE]
		mat.ref = cluster.data
		
		list.runs = list()
		for(ii in 1:ncol(mat.ref)) list.runs[[ii]] = cluster.data[,grep(ii, Design$Run)]
		
		list.scaledruns = lapply(1:length(list.runs), function(x) list.runs[[x]]/mat.ref[,x])
		
		scaled.cluster.data = list.scaledruns[[1]]
		if(length(list.scaledruns) > 1) 
		for(ii in 2:length(list.scaledruns)) scaled.cluster.data = cbind(scaled.cluster.data, list.scaledruns[[ii]])
		
		#if(KeepREF == FALSE) scaled.cluster.data = 
		#	scaled.cluster.data[,!colnames(scaled.cluster.data)%in%refInRuns]
		
		plotClusterProfile(log(scaled.cluster.data), clustID, Group, k=noClusters, ylab="Average log ratio")

	}

	Clusters = rep(NA, nrow(data_irs))
	Clusters[Anova.idx] = Cluster

	# reorder abundance data
	data_irs_orderbygroup = data_irs[,order(Group)]
	colnames(data_irs_orderbygroup) = gsub("\\.", "\\.Abundance.", colnames(data_irs_orderbygroup))
	
	Uniprot = gsub("(.*)\\|(.*)\\|(.*)", "\\2", rownames(data_irs_orderbygroup))
	Uniprot = gsub("\\-.", "", Uniprot)
	full.res = data.frame(Accession=rownames(data_irs_orderbygroup),
	Uniprot ,
  	data_irs_orderbygroup, 
		Means, MaxFC, Anova, Anova.adj, Clusters)
		
	
	full.res = data.frame(full.res, mg.all)
	
	# Output all samples and corresponding groups in the new order
	write.csv(data.frame(Sample=colnames(data_irs_orderbygroup), Group=Group[order(Group)]), "samplegroup.csv")
	
	
	
	cat('end clustering\n')
	################
	# End ANOVA
	################

	
	# Output results
	
	wb <- createWorkbook()
	full.res[full.res=='NaN'] = NA

	# add parameter cutoffs
	addWorksheet(wb, sheet='Parameter')
	writeData(wb, 'Parameter', dat.para)
		
	
	
	ps <- try(printOpenxlsxStyle(full.res, ratios=grep("MaxFC", names(full.res)), 
					pvals=c(grep("Anova", names(full.res))), 
					wb = wb, 
					tabName = "AllData", hiCutoff=FCCutoff, lowCutoff=1/FCCutoff, pvalCutoff=PvalCutoff) )

	if(inherits(ps, 'try-error') ) warning('Error with print overall xlsx file')
					
	ps <- try(printOpenxlsxStyle(data.frame(rownames(pca.res$componentScores),pca.res$componentScores), ratios=NULL, pvals=NULL, wb = wb, tabName = "PCAScores"))

	if(inherits(ps, 'try-error') ) warning('Error with print overall PCA component scores tab')

	ps <- try(printOpenxlsxStyle(data.frame(rownames(ld), ld), ratios=NULL, pvals=NULL, wb = wb, tabName = "PCALoadings"))

	if(inherits(ps, 'try-error') ) warning('Error with print overall PCA loading tab')


	saveWorkbook(wb, file="ResultsOverall.xlsx", overwrite=TRUE)

	
	
	##############
	# Pairwise 
	##############

	# if(sum(grepl("comp", tolower(names(designSheets)) )) == 1 ) {
	
	dat.comparisons = Comparisons

	tarcompres.list = list()

	for(idx.comp in 1:nrow(dat.comparisons) ) {
		comp = dat.comparisons[idx.comp,-1]
		
		idx.group1 = which(Group %in% comp[1,1])

		idx.group2 = which(Group %in% comp[1,2])

		idx.mean1 = which(comp[1,1] == gsub('Means ', '', colnames(Means), fixed=TRUE) ) 
		  
		idx.mean2 = which(comp[1,2] == gsub('Means ', '', colnames(Means), fixed=TRUE) )

		FC = Means[,idx.mean1]/Means[,idx.mean2]
		
		TwoSplTTest = rep(NA, nrow(data_irs))
		
		for(ii in 1:nrow(data_irs)) {	
			v = t(data_irs[ii,])
			
			# TWo sample t test
			temp = try(t.test(log(na.omit(v[idx.group1])), 
				log(na.omit(v[idx.group2])), var.equal=TRUE) )
				
			if(!inherits(temp, "try-error")) TwoSplTTest[ii] = temp$p.value		
		
		}
		
		TwoSplTTestAdj = p.adjust(TwoSplTTest, 'fdr')
		
		idx.sig = (!is.na(TwoSplTTest)) & (!is.na(FC)) & (TwoSplTTest < PvalCutoff) & (abs(log(FC))) > log(FCCutoff)
	
		png(paste('Volcano plot for targeted', dat.comparisons[idx.comp,1], '.png', sep=''), 2000, 2000, res=300)
		plot(log(FC), -log(TwoSplTTest), xlab='log FC', ylab='-log p value', 
					main=paste('Protein volcano plot', dat.comparisons[idx.comp,1]))
		abline(h=-log(PvalCutoff), col="red")
		abline(v=log(FCCutoff), col="blue")
		abline(v=-log(FCCutoff), col="blue")
		
		if(sum(idx.sig) > 0) 			
			points(log(FC[idx.sig]), -log(TwoSplTTest[idx.sig]), col='red', pch=20)

		dev.off()	

		
		dat.annot = mg.all[match(rownames(data_irs), mg.all$Protein.Id), c(grep("Description", colnames(mg.all)), grep("Number.of.peptides", colnames(mg.all)))]
		
		dat.comp = data.frame(Accession=rownames(data_irs), dat.annot, data_irs[,c(idx.group1, idx.group2)], Means[,c(idx.mean1, idx.mean2)],
					FC, TwoSplTTest, TwoSplTTestAdj, Significant = idx.sig)

		tarcompres.list[[idx.comp]] = list(Result=dat.comp, ComparisonName=dat.comparisons[idx.comp,1] )

	}
	
	wb <- createWorkbook()
	
	# add parameter cutoffs
	addWorksheet(wb, sheet='Parameter')
	writeData(wb, 'Parameter', dat.para)
	
	for(idx.comp in 1:length(tarcompres.list)) {
	
		dd = tarcompres.list[[idx.comp]][[1]]
		dd[dd=='NaN'] = NA
		
		colnames(dd)[1] = 'Accession'
		
		# sort by significant
		dd = dd[order(dd$Significant, decreasing=TRUE), ]
		
		printOpenxlsxStyle(dd, ratios=grep('FC', names(tarcompres.list[[idx.comp]][[1]])),
				pvals=grep('TwoSplTTest', names(tarcompres.list[[idx.comp]][[1]])), wb=wb,
				tabName = tarcompres.list[[idx.comp]]$ComparisonName, hiCutoff=FCCutoff, lowCutoff=1/FCCutoff, pvalCutoff=PvalCutoff)
	}
	
		
	# images
	addWorksheet(wb, sheet='images')	
	startCol = 1
	for(idx.comp in 1:length(tarcompres.list)) {
		if(file.exists(paste('Volcano plot for targeted', dat.comparisons[idx.comp,1], '.png', sep=''))) {			
			insertImage(wb, sheet='images', 
				file=paste('Volcano plot for targeted', dat.comparisons[idx.comp,1], '.png', sep=''), 
				width=12, height=15, startRow=2, startCol=startCol, units='cm')
			
			startCol = startCol + 10
		}
	}
	
	# comparisons
	addWorksheet(wb, sheet='comparisons')
	writeData(wb, sheet='comparisons', dat.comparisons)
	
	saveWorkbook(wb, file="ResultsTargetedPaiwise.xlsx", overwrite=TRUE)
	
	

	# generate MD5 checksums
		
	flist = list.files()

	md5list = lapply(flist, FUN=md5sum)

	# format •	MD5 checksum
	#        •	two spaces ('  ')
	#        •	the filename

	write.table(data.frame(unlist(md5list), names(unlist(md5list))), sep="  ", quote=F, row.names=F, file = "checksums.txt")
	
		
}	




