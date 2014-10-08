#!/usr/bin/env Rscript
library("lattice");
library("grid");
#options(error=recover);

# Number of decimal places to round non-carrier phenotype values to when drawing on the plot.
# This helps reduce the number of unnecessary points plotted.
ROUND_NONCARRIER = 3;

#function qqplot (taken from the sph wiki)
qqunif.plot = function(pvalues,
	should.thin=T, thin.obs.places=2, thin.exp.places=2, 
	xlab=expression(paste("Expected (",-log[10], " p-value)")),
	ylab=expression(paste("Observed (",-log[10], " p-value)")), 
	draw.conf=TRUE, conf.points=1000, conf.col="lightgray", conf.alpha=.05,
	already.transformed=FALSE, pch=20, aspect="iso", 
	par.settings=list(superpose.symbol=list(pch=pch)), ...) {

	#error checking
	if (length(pvalues)==0) stop("pvalue vector is empty, can't draw plot")
	if(!(class(pvalues)=="numeric" || 
		(class(pvalues)=="list" && all(sapply(pvalues, class)=="numeric"))))
		stop("pvalue vector is not numeric, can't draw plot")
	if (any(is.na(unlist(pvalues)))) stop("pvalue vector contains NA values, can't draw plot")
	if (already.transformed==FALSE) {
		if (any(unlist(pvalues)==0)) stop("pvalue vector contains zeros, can't draw plot")
	} else {
		if (any(unlist(pvalues)<0)) stop("-log10 pvalue vector contains negative values, can't draw plot")
	}

	grp<-NULL
	n<-1
	exp.x<-c()
	if(is.list(pvalues)) {
		nn<-sapply(pvalues, length)
		rs<-cumsum(nn)
		re<-rs-nn+1
		n<-min(nn)
		if (!is.null(names(pvalues))) {
			grp=factor(rep(names(pvalues), nn), levels=names(pvalues))
			names(pvalues)<-NULL
		} else {
			grp=factor(rep(1:length(pvalues), nn))
		}
		pvo<-pvalues
		pvalues<-numeric(sum(nn))
		exp.x<-numeric(sum(nn))
		for(i in 1:length(pvo)) {
			if (!already.transformed) {
				pvalues[rs[i]:re[i]] <- -log10(pvo[[i]])
				exp.x[rs[i]:re[i]] <- -log10((rank(pvo[[i]], ties.method="first")-.5)/nn[i])
			} else {
				pvalues[rs[i]:re[i]] <- pvo[[i]]
				exp.x[rs[i]:re[i]] <- -log10((nn[i]+1-rank(pvo[[i]], ties.method="first")-.5)/(nn[i]+1))
			}
		}
	} else {
		n <- length(pvalues)+1
		if (!already.transformed) {
			exp.x <- -log10((rank(pvalues, ties.method="first")-.5)/n)
			pvalues <- -log10(pvalues)
		} else {
			exp.x <- -log10((n-rank(pvalues, ties.method="first")-.5)/n)
		}
	}


	#this is a helper function to draw the confidence interval
	panel.qqconf<-function(n, conf.points=1000, conf.col="gray", conf.alpha=.05, ...) {
		require(grid)
		conf.points = min(conf.points, n-1);
		mpts<-matrix(nrow=conf.points*2, ncol=2)
        	for(i in seq(from=1, to=conf.points)) {
            		mpts[i,1]<- -log10((i-.5)/n)
            		mpts[i,2]<- -log10(qbeta(1-conf.alpha/2, i, n-i))
            		mpts[conf.points*2+1-i,1]<- -log10((i-.5)/n)
            		mpts[conf.points*2+1-i,2]<- -log10(qbeta(conf.alpha/2, i, n-i))
        	}
        	grid.polygon(x=mpts[,1],y=mpts[,2], gp=gpar(fill=conf.col, lty=0), default.units="native")
    	}

	#reduce number of points to plot
	if (should.thin==T) {
		if (!is.null(grp)) {
			thin <- unique(data.frame(pvalues = round(pvalues, thin.obs.places),
				exp.x = round(exp.x, thin.exp.places),
				grp=grp))
			grp = thin$grp
		} else {
			thin <- unique(data.frame(pvalues = round(pvalues, thin.obs.places),
				exp.x = round(exp.x, thin.exp.places)))
		}
		pvalues <- thin$pvalues
		exp.x <- thin$exp.x
	}
	gc()
 
	prepanel.origin = function(x,y,...) {
		A = list()
		A$xlim = range(x, y)*1.02
		A$xlim[1]=0
		A$ylim = A$xlim
		return(A)
	}
 
	#draw the plot
	xyplot(pvalues~exp.x, groups=grp, xlab=xlab, ylab=ylab, aspect=aspect,
		prepanel=prepanel.origin, scales=list(axs="i"), pch=pch,
		panel = function(x, y, ...) {
			if (draw.conf) {
				panel.qqconf(n, conf.points=conf.points, 
					conf.col=conf.col, conf.alpha=conf.alpha)
			};
			panel.xyplot(x,y, ...);
			panel.abline(0,1);
		}, par.settings=par.settings, ...
	)
}

#function for manhattan plot: from the sph wiki
manhattan.plot = function(chr, pos, pvalue,
	sig.level=NA, annotate=NULL, ann.default=list(),
	should.thin=T, thin.pos.places=2, thin.logp.places=2, 
	xlab="Chromosome", ylab=expression(-log[10](p-value)),
	col=c("gray","darkgray"), panel.extra=NULL, pch=20, cex=0.8,...) {
 
	if (length(chr)==0) stop("chromosome vector is empty")
	if (length(pos)==0) stop("position vector is empty")
	if (length(pvalue)==0) stop("pvalue vector is empty")
 
	#make sure we have an ordered factor
	if(!is.ordered(chr)) {
		chr <- ordered(chr)
	} else {
		chr <- chr[,drop=T]
	}
 
	#make sure positions are in kbp
	if (any(pos>1e6)) pos<-pos/1e6;
 
	#calculate absolute genomic position
	#from relative chromosomal positions
	posmin <- tapply(pos,chr, min);
	posmax <- tapply(pos,chr, max);
	posshift <- head(c(0,cumsum(posmax)),-1);
	names(posshift) <- levels(chr)
	genpos <- pos + posshift[chr];
	getGenPos<-function(cchr, cpos) {
		p<-posshift[as.character(cchr)]+cpos
		return(p)
	}
 
	#parse annotations
	grp <- NULL
	ann.settings <- list()
	label.default<-list(x="peak",y="peak",adj=NULL, pos=3, offset=0.5, 
		col=NULL, fontface=NULL, fontsize=NULL, show=F)
	parse.label<-function(rawval, groupname) {
		r<-list(text=groupname)
		if(is.logical(rawval)) {
			if(!rawval) {r$show <- F}
		} else if (is.character(rawval) || is.expression(rawval)) {
			if(nchar(rawval)>=1) {
				r$text <- rawval
			}
		} else if (is.list(rawval)) {
			r <- modifyList(r, rawval)
		}
		return(r)
	}
 
	if(!is.null(annotate)) {
		if (is.list(annotate)) {
			grp <- annotate[[1]]
		} else {
			grp <- annotate
		} 
		if (!is.factor(grp)) {
			grp <- factor(grp)
		}
	} else {
		grp <- factor(rep(1, times=length(pvalue)))
	}
 
	ann.settings<-vector("list", length(levels(grp)))
	ann.settings[[1]]<-list(pch=pch, col=col, cex=cex, fill=col, label=label.default)
 
	if (length(ann.settings)>1) { 
		lcols<-trellis.par.get("superpose.symbol")$col 
		lfills<-trellis.par.get("superpose.symbol")$fill
		for(i in 2:length(levels(grp))) {
			ann.settings[[i]]<-list(pch=pch, 
				col=lcols[(i-2) %% length(lcols) +1 ], 
				fill=lfills[(i-2) %% length(lfills) +1 ], 
				cex=cex, label=label.default);
			ann.settings[[i]]$label$show <- T
		}
		names(ann.settings)<-levels(grp)
	}
	for(i in 1:length(ann.settings)) {
		if (i>1) {ann.settings[[i]] <- modifyList(ann.settings[[i]], ann.default)}
		ann.settings[[i]]$label <- modifyList(ann.settings[[i]]$label, 
			parse.label(ann.settings[[i]]$label, levels(grp)[i]))
	}
	if(is.list(annotate) && length(annotate)>1) {
		user.cols <- 2:length(annotate)
		ann.cols <- c()
		if(!is.null(names(annotate[-1])) && all(names(annotate[-1])!="")) {
			ann.cols<-match(names(annotate)[-1], names(ann.settings))
		} else {
			ann.cols<-user.cols-1
		}
		for(i in seq_along(user.cols)) {
			if(!is.null(annotate[[user.cols[i]]]$label)) {
				annotate[[user.cols[i]]]$label<-parse.label(annotate[[user.cols[i]]]$label, 
					levels(grp)[ann.cols[i]])
			}
			ann.settings[[ann.cols[i]]]<-modifyList(ann.settings[[ann.cols[i]]], 
				annotate[[user.cols[i]]])
		}
	}
 	rm(annotate)
 
	#reduce number of points plotted
	if(should.thin) {
		thinned <- unique(data.frame(
			logp=round(-log10(pvalue),thin.logp.places), 
			pos=round(genpos,thin.pos.places), 
			chr=chr,
			grp=grp)
		)
		logp <- thinned$logp
		genpos <- thinned$pos
		chr <- thinned$chr
		grp <- thinned$grp
		rm(thinned)
	} else {
		logp <- -log10(pvalue)
	}
	rm(pos, pvalue)
	gc()
 
	#custom axis to print chromosome names
	axis.chr <- function(side,...) {
		if(side=="bottom") {
			panel.axis(side=side, outside=T,
				at=((posmax+posmin)/2+posshift),
				labels=levels(chr), 
				ticks=F, rot=0,
				check.overlap=F
			)
		} else if (side=="top" || side=="right") {
			panel.axis(side=side, draw.labels=F, ticks=F);
		}
		else {
			axis.default(side=side,...);
		}
	 }
 
	#make sure the y-lim covers the range (plus a bit more to look nice)
	prepanel.chr<-function(x,y,...) { 
		A<-list();
		maxy<-ceiling(max(y, ifelse(!is.na(sig.level), -log10(sig.level), 0)))+.5;
		A$ylim=c(0,maxy);
		A;
	}
 
	xyplot(logp~genpos, chr=chr, groups=grp,
		axis=axis.chr, ann.settings=ann.settings, 
		prepanel=prepanel.chr, scales=list(axs="i"),
		panel=function(x, y, ..., getgenpos) {
			if(!is.na(sig.level)) {
				#add significance line (if requested)
				panel.abline(h=-log10(sig.level), lty=2);
			}
			panel.superpose(x, y, ..., getgenpos=getgenpos);
			if(!is.null(panel.extra)) {
				panel.extra(x,y, getgenpos, ...)
			}
		},
		panel.groups = function(x,y,..., subscripts, group.number) {
			A<-list(...)
			#allow for different annotation settings
			gs <- ann.settings[[group.number]]
			A$col.symbol <- gs$col[(as.numeric(chr[subscripts])-1) %% length(gs$col) + 1]    
			A$cex <- gs$cex[(as.numeric(chr[subscripts])-1) %% length(gs$cex) + 1]
			A$pch <- gs$pch[(as.numeric(chr[subscripts])-1) %% length(gs$pch) + 1]
			A$fill <- gs$fill[(as.numeric(chr[subscripts])-1) %% length(gs$fill) + 1]
			A$x <- x
			A$y <- y
			do.call("panel.xyplot", A)
			#draw labels (if requested)
			if(gs$label$show) {
				gt<-gs$label
				names(gt)[which(names(gt)=="text")]<-"labels"
				gt$show<-NULL
				if(is.character(gt$x) | is.character(gt$y)) {
					peak = which.max(y)
					center = mean(range(x))
					if (is.character(gt$x)) {
						if(gt$x=="peak") {gt$x<-x[peak]}
						if(gt$x=="center") {gt$x<-center}
					}
					if (is.character(gt$y)) {
						if(gt$y=="peak") {gt$y<-y[peak]}
					}
				}
				if(is.list(gt$x)) {
					gt$x<-A$getgenpos(gt$x[[1]],gt$x[[2]])
				}
				do.call("panel.text", gt)
			}
		},
		xlab=xlab, ylab=ylab, 
		panel.extra=panel.extra, getgenpos=getGenPos, ...
	);
	#abline(h=log(0.0000025))
}

#function addspaces: adds spaces to variables in order to output columns of fixed width in the pdf
addspaces <- function(string, maxstring){
	if(nchar(string) < nchar(maxstring)){
		for(i in 1:(nchar(maxstring) - nchar(string))){
			string = paste(string, " ", sep="")
		}
	}
	return(string)
}

##############################################################################
# MAIN
##############################################################################

#get arguments from input
arguments = commandArgs(trailingOnly=TRUE)
#variants = read.table(arguments[1], header=T,stringsAsFactors=FALSE,sep="\t")
#variants$PVALUE = as.numeric(variants$PVALUE)
phenotype = arguments[2]

prefix = arguments[3]				#this is the argument for the prefix of the filename to read
tests = arguments[4]
pvaluethreshold = as.numeric(arguments[5])
tests = strsplit(tests, ',')
covariates = arguments[6]

#read files
phenofile = read.table(paste(prefix, '.pheno.ped',sep=''),stringsAsFactors=F, header=T,comment.char="")
summaryfile = read.table(paste(prefix,'.summary.txt',sep=''),sep="\t",stringsAsFactors=FALSE)

allvariants = read.table(paste(prefix, '.allvariants.txt',sep=''),sep="\t",stringsAsFactors=FALSE,header=T)
allvariants$PVALUE = as.numeric(allvariants$PVALUE)

#phenotypes
phenos = read.table(paste(prefix,".pheno.ped",sep=""), comment.char="", header=T,stringsAsFactors=F)

phenotypes =  as.numeric(phenos[!is.na(phenos[phenotype]),phenotype])
phenotypemean = mean(phenotypes)

#get phenotype residuals if covariates are given, else set phenotype residuals to phenotypes
if(covariates != 'NA'){
	model = paste(phenotype , covariates, sep=' ~ ')
	phenotyperesiduals = lm(model, phenofile)$residuals
	j = 1
	for(i in 1:nrow(phenos)){
		if(is.na(phenos[i, phenotype])){TEMPORARY = 1;}
		else{phenos[i,phenotype] = phenotyperesiduals[j]; j = j +1}
	}
}else{
	phenotyperesiduals = phenotypes
}

merged = merge(allvariants, phenos, by.x = 'IND', by.y = 'IND_ID')

#number of columns to add to the pdf on the right.
cols_to_annotate = allvariants
cols_to_annotate$GROUP = NULL
cols_to_annotate$MARKER_ID = NULL
cols_to_annotate$MARKERNAME = NULL
cols_to_annotate$GENOTYPE = NULL
cols_to_annotate$IND = NULL
cols_to_annotate$MAF = NULL
cols_to_annotate$MAC = NULL
cols_to_annotate$BETA = NULL
cols_to_annotate$PVALUE = NULL

#after deleting these columns, whatever is left in columnstoannotate has to be added to the info on the right. We have to compute how many such columns there are, and how much space it is likely to take up
extralengthtoadd = 0
extracolumns = NULL
if(ncol(cols_to_annotate) > 0){
	extracolumns = colnames(cols_to_annotate);
	lengthextracolumns = NA
	for(column in extracolumns){
		len = cols_to_annotate[which.max(nchar(cols_to_annotate[column][[1]])),column];
		lengthextracolumns = append(lengthextracolumns, len)
		extralengthtoadd = extralengthtoadd + max(nchar(cols_to_annotate[column][[1]]))
	}
	lengthextracolumns = lengthextracolumns[-1]
}

#for each additional groupwise-test performed, add extra columns
if(length(tests[[1]]) > 1){
	for(i in 2:length(tests[[1]])){
		extralengthtoadd = extralengthtoadd + 7
	}
}

#for each test, get the epacts burden test results and put them all in a data frame called epacts1
for(i in 1:length(tests[[1]])){
	epacts1 = read.table(paste(prefix, '.', tests[[1]][i] ,".epacts",sep=""),header=T, comment.char="",stringsAsFactors=F,sep="\t")
	if(!('BEGIN' %in% colnames(epacts1))){
		colnum = which(colnames(epacts1) == 'BEG')
		colnames(epacts1)[colnum] = 'BEGIN'
	}
	if(i == 1){
		epacts = epacts1;
		genepvalues = epacts1[,c('MARKER_ID', 'PVALUE')]
		genepvalues['TEST'] = tests[[1]][i]
		epacts[tests[[1]][i]] = epacts1['PVALUE']
	}else{
		genepvalues1 = epacts1[,c('MARKER_ID', 'PVALUE')]
		genepvalues1['TEST'] = tests[[1]][i]
		genepvalues = rbind(genepvalues, genepvalues1)
		epacts[tests[[1]][i]] = epacts1['PVALUE']
		
	}
}

merged = merge(merged, epacts, by.x = 'GROUP', by.y = 'MARKER_ID')
merged$'MAC' = round(merged$'MAC',2)
#sort by BETA, then by pvalue for the gene, then by pvalue of the individual markers
if('BETA' %in% colnames(merged)){
	merged = merged[order(as.numeric(merged$'PVALUE.y'), merged$'PVALUE.x', -abs(merged$'BETA')),];
}

for(i in 1:length(tests[[1]])){
	merged[tests[[1]][i]] = apply(merged[tests[[1]][i]], 1, function(x) sprintf("%0.3g",as.numeric(x)))
}

merged$'BETA' = format(round(merged$'BETA',2),nsmall=2)
merged$'PVALUE.x' = sprintf("%0.3g",merged$'PVALUE.x')
merged$'PVALUE.y' = sprintf("%0.3g",merged$'PVALUE.y')
#merged$'BETA' = format(round(merged$'BETA',2),nsmall=2)
merged$'MAF' = sprintf("%0.3g",merged$'MAF')

allgenes = unique(merged[,c('GROUP','PVALUE.y')])
allgenes = allgenes[which(as.numeric(allgenes[,2]) <= pvaluethreshold),]
allgenes = allgenes[order(allgenes[,2]),1]
allgenes = unique(allgenes)

genes = unique(merged[,c('GROUP', 'PVALUE.y')])
genes = genes[order(genes[,2]),]
genes = genes[which(as.numeric(genes[,2]) <= pvaluethreshold),]
genes = unique(genes[,1])
totalgenes = length(genes)

genes_passing_filters = read.table(paste(prefix,'.genes_passing_filters.txt',sep=""),header=T,stringsAsFactors=F)
merged['GENENAME'] = unlist(strsplit(merged$GROUP, "_"))[seq(2,2*nrow(merged),2)]
filteredresults = merge(merged, genes_passing_filters, by.x='GENENAME', by.y=1,sort=F)
genes = filteredresults[,c('GROUP', 'PVALUE.y')]
genes = genes[order(genes[,2]),]
genes = genes[which(as.numeric(genes[,2]) <= pvaluethreshold),]
genes = unique(genes[,1])
totalgenes = length(genes)

# Create a version of the PDF's table for a file
test_cols = unlist(tests)
# ddply version (don't want to have to include this as dependency)
#merged_out = ddply(merged,"GROUP",function(x) x[!duplicated(x$MARKER_ID),c("GROUP",test_cols,"MARKER_ID","BETA","PVALUE.x","MAF","MAC",extracolumns)])
#merged_out = rename(merged_out,c("PVALUE.x" = "VAR.P"))
merged_out = do.call(rbind,by(merged,merged$GROUP,function(x) x[!duplicated(x$MARKER_ID),c("GROUP",test_cols,"MARKER_ID","BETA","PVALUE.x","MAF","MAC",extracolumns)]))
names(merged_out)[names(merged_out) == "PVALUE.x"] = "VAR.P"
for (itest in test_cols) {
  names(merged_out)[names(merged_out) == itest] = paste(itest,".P",sep="")
}
write.table(merged_out,file=paste(prefix,".topgenes",'.tab',sep=""),sep="\t",row.names=F,quote=F)

categorical = 0
if (length(unique(phenotypes)) == 2){
	categorical = 1;
}

#write out information for carriers to a file

if(categorical == 0){				#quantitative phenotype
		print("Quantitative phenotype")
		numgenestoinclude = 0
		totalvariants = 0
		length_longest_variant = 0
		longest_variant = ""
		longest_beta = ""
		longest_mac = ""

		print(genes)
		#this for loop gets the number of genes to plot as a histogram, and how many pages to split it into. maximum of 40 rows per page.
		if(length(genes) > 0){
			for(genenum in 1:length(genes)){
				gene = genes[genenum]
				rowstoplot = merged[merged$'GROUP' == gene,]
				num_variants = length(unique(rowstoplot$'MARKER_ID'))
				
				totalvariants = totalvariants + num_variants
				markersingene = unique(rowstoplot$'MARKER_ID')
				if((numgenestoinclude == 0) || (totalvariants <= 40)){
					numgenestoinclude = numgenestoinclude + 1;
					l = nchar(markersingene[which.max(nchar(markersingene))])
					if (l > length_longest_variant){
						length_longest_variant = l; 
						longest_variant = markersingene[which.max(nchar(unique(rowstoplot$'MARKER_ID')))];
					}
					if('BETA' %in% colnames(rowstoplot)){
						l = rowstoplot$'BETA'[which.max(nchar(as.numeric(rowstoplot$'BETA')))]
						if(nchar(l) > nchar(longest_beta)){longest_beta = l;}
					}
					l = rowstoplot$'MAC'[which.max(nchar(rowstoplot$'MAC'))]
					if(nchar(l) > nchar(longest_mac)){longest_mac = l;}
				}
			}
		}
		print(numgenestoinclude)
		if(numgenestoinclude > 1){
			totalvariants = min(totalvariants, 40)
		}

		numgenesplotted = 0
		numiter = 1
		print("Creating PDF...")
		pdf(paste(prefix,".topgenes",'.pdf',sep=""), width = 11, height = 8.5)
		widthinches = unit(11, "inches")
		maxwidth = convertUnit(widthinches, "char")
		totalsummary = paste(summaryfile[,1],summaryfile[,2],sep=" : ")
		maxchars1 = nchar(summaryfile[which.max(nchar(summaryfile[,1])),1])
		maxchars2 = nchar(summaryfile[which.max(nchar(summaryfile[,2])),2])
		maxchars = nchar(totalsummary[which.max(nchar(totalsummary))])
		cextouse = 1
		#if(maxchars > as.numeric(maxwidth)){
		#	cextouse = 0.5;
		#}
		
		outer_layout = grid.layout(
		nrow=nrow(summaryfile)+1,
		ncol = 1,
		widths =  unit(1, 'npc'),
		heights = unit.c(rep(unit(1, 'char'), nrow(summaryfile)), unit(1,'null'))
		);

		pushViewport(viewport(
		layout = outer_layout,
		x = unit(0.05,'npc'),
		y = unit(0,'npc'),
		height = unit(0.975,'npc'),
		width = unit(0.95,'npc'),
		just = c("left","bottom"),
		name = "outer"
		));
		
		for(i in 1:nrow(summaryfile)){
			pushViewport(viewport(
				layout.pos.row = i,
				layout = grid.layout(
						nrow = 1,
						ncol= 2,
						widths = unit.c(unit(0.1, 'npc'), unit(0.9, 'npc'))
					)
			));

			pushViewport(viewport(
				layout.pos.col=1
			));
			
			grid.text(
				paste(summaryfile[i,1],sep=""),
				x= 0.00,
				y= 0.5,
				gp = gpar(cex=0.9, fontsize=8),
				just = 'left'
			);
			#grid.rect();
			upViewport(1);

			pushViewport(viewport(
				layout.pos.col=2
			));

			grid.text(
				paste(summaryfile[i,2],sep=""),
				x= 0.02,
				y= 0.5,
				gp = gpar(cex=0.9),
				just='left'
			);
		#	grid.rect();
			upViewport(1);
			upViewport(1);
		}
		
		for(i in 1:length(tests[[1]])){
			test = tests[[1]][i]
			#plot qq plot
			#par(mfrow = c(1,2))
			if(i == 1){
				masterlayout = grid.layout(
					nrow = 1, 
					ncol = 3,
					widths = unit.c(unit(0.47, "npc"), unit(0.07, "npc"), unit(0.47, "npc"))
				)
				vp1 = viewport(layout.pos.col = 1, name="vp1")
				vp2 = viewport(layout.pos.col = 3, name="vp2")
				vp3 = viewport(layout.pos.col = 2, name="vp3")
			}
			epactsfortest = read.table(paste(prefix,'.', tests[[1]][i] ,".epacts",sep=""),header=T, comment.char="",stringsAsFactors=F,sep="\t")
			if (!('BEGIN' %in% colnames(epactsfortest))){
				colnum = which(colnames(epactsfortest) == 'BEG')
				colnames(epactsfortest)[colnum] = 'BEGIN'
			}
			#get genes
			##If there is no filter, then genes_passing_filters = epactsfortest['MARKER_ID']
			genes_passing_filters = read.table(paste(prefix,'.genes_passing_filters.txt',sep=""),header=T,stringsAsFactors=F)
			epactsfortest['GENENAME'] = unlist(strsplit(epactsfortest$MARKER_ID, "_"))[seq(2,2*nrow(epactsfortest),2)]
			filteredresults = merge(epactsfortest, genes_passing_filters, by.x='GENENAME', by.y=1,sort=F)
			
			#plotting all variants
			variantstoplot = epactsfortest['PVALUE']
			toplot = variantstoplot['PVALUE']
			toplot = toplot[!is.na(toplot)]
			qq1plot = toplot
			qq1.plot = qqunif.plot(toplot, main=paste("Test: ", test, " Unfiltered QQ Plot\n", phenotype, "\nNo. Samples=", length(phenotypes) ,", No. Genes=",length(toplot),sep=""),cex.main=0.5)
			#update(qq.plot,par.settings = list(fontsize = list(text = 10, points = 8)))

			#plotting filtered variants
			toplot = filteredresults[,c('MARKER_ID','PVALUE')]
			toplot = unique(toplot)
			toplot = toplot[,2]
			toplot = toplot[!is.na(toplot)]
			qq2.plot = qqunif.plot(toplot, main=paste("Test: ", test, " Filtered QQ Plot\n", phenotype, "\nNo. Samples=", length(phenotypes) ,", No. Genes=",length(toplot),sep=""),cex.main=0.5)

			grid.newpage()
			pushViewport(vpTree(viewport(layout = masterlayout,name="master"), vpList(vp1, vp2)))
			seekViewport("master")
			print(qq1.plot, draw.in = "vp1")
			print(qq2.plot, draw.in = "vp2")
			upViewport(1);
			#update(qq.plot,par.settings = list(fontsize = list(text = 10, points = 8)))

			variantstoplot = epactsfortest[,c('X.CHROM', 'BEGIN', 'PVALUE')]
			if (length(which(variantstoplot[,1] == 'X')) > 0){variantstoplot[which(variantstoplot[,1] == 'X'),1] = '23';}
			if (length(which(variantstoplot[,1] == 'Y')) > 0){variantstoplot[which(variantstoplot[,1] == 'Y'),1] = '24';}

			variantstoplot = variantstoplot[which(!is.na(variantstoplot[,3])),]
			toplot = cbind(as.numeric(variantstoplot[,1]), variantstoplot[,2], variantstoplot['PVALUE'])
			toplot = toplot[order(toplot[,1], toplot[,2]),]
			temppp = toplot
			mhtplot = manhattan.plot(toplot[[1]], toplot[[2]], toplot[[3]], sig.level=2e-6,main=paste("Test: ", test, " Unfiltered. ",phenotype, ", No. Samples=",length(phenotypes),", No. Genes=",nrow(toplot),sep=""))
			print(mhtplot)
			
			##draw filtered manhattan plot
			toplot = filteredresults[,c('X.CHROM', 'BEGIN','PVALUE','MARKER_ID')]
			toplot = toplot[which(!is.na(toplot[,3])),]
			toplot = unique(toplot)
			if (length(which(toplot[,1] == 'X')) > 0){toplot[which(toplot[,1] == 'X'),1] = '23';}
			if (length(which(toplot[,1] == 'Y')) > 0){toplot[which(toplot[,1] == 'Y'),1] = '24';}
			toplot[,1] = as.numeric(toplot[,1])
			toplot[,2] = as.numeric(toplot[,2])
			toplot = toplot[order(toplot[,1], toplot[,2]),]
			mhtplot = manhattan.plot(toplot[,1], toplot[,2], toplot[,3], sig.level=2e-6,main=paste("Test: ", test, " Filtered ", phenotype, ", No. Samples=",length(phenotypes),", No. Genes=",nrow(toplot),sep=""))

			print(mhtplot)
		}

		while(numgenesplotted < totalgenes){
			numgenesplotted = numgenesplotted + numgenestoinclude
			genes = genes[1:numgenestoinclude]

			longestgenename = gsub(".*_","",genes[which.max(nchar(unique(gsub(".*_","",genes,perl=T))))],perl=T)
			grid.newpage();

			xscale = lattice:::extend.limits(range(phenotyperesiduals,finite=TRUE));
			nint = if (is.factor(phenotyperesiduals)) nlevels(phenotyperesiduals) else round(log2(length(phenotyperesiduals)) + 3);
			if (is.factor(phenotypes)) {
			breaks = seq_len(1 + nlevels(phenotyperesiduals)) - 0.5;
			} else {
			breaks = do.breaks(xscale,nint);
			}

			# Calculate y scale assuming we're showing a histogram in percent mode. 
			max_y = 100 * max(table(cut(phenotyperesiduals,breaks))) / length(phenotyperesiduals);
			yscale = lattice:::extend.limits(c(0,max_y));

			# Setup our layouts for the plotting area. 
			axis_height = 1; # in units of 'lines', this is how tall the variant lines are

			widthtouse = strwidth(longest_variant,units='in',cex=0.6) + 
				strwidth(longestgenename,units='in',cex=0.6) + 
				strwidth(paste(rep('a', (extralengthtoadd+20)),collapse=""),units='in',cex=0.6) + 
				strwidth(paste(rep('a', max(length(longest_beta)+1, nchar('BETA '))),collapse=""),units='in',cex=0.6) + 
				strwidth(paste(rep('a', max(nchar('4.0e-07'), nchar('VAR.P '))),collapse=""),units='in',cex=0.6) + 
				strwidth(paste(rep('a', max(length('4.0e-07')+1, nchar('MAF '))),collapse=""),units='in',cex=0.6) +
				strwidth(paste(rep('a', max(nchar(longest_mac)+1, length('MAC '))),collapse=""),units='in',cex=0.6) 
			for(testnum in 1:length(tests[[1]])){
				widthtouse = widthtouse + strwidth('4.0e-07',units='in',cex=0.6);
			}

			# Create our outer layout. This is just 3 columns, where the left column
			# leaves room for the y-axis label, the middle column is for the histogram
			# and variant ticks/axes, and the last column will have the variant names. 
			outer_layout = grid.layout(
				nrow = 1,
				ncol = 4,
				widths = unit.c(
					unit(3,'lines'),
					unit(1,'null'),
					unit(widthtouse, 'in'),
					unit(2, 'lines')
				)
			);
			# Create a viewport using the layout above. 
			pushViewport(viewport(
				layout = outer_layout,
				x = unit(0,'npc'),
				y = unit(0,'npc'),
				height = unit(0.975,'npc'),
				width = unit(1,'npc'),
				just = c("left","bottom"),
				name = "outer"
			));

			# Now the inner layout (within the 2nd column of the outer layout). 
			# This has 1 row for the histogram, a spacer row, and then however many 
			inner_layout = grid.layout(
			# rows are needed for each variant. 
				nrow = 1 + 1 +1 + totalvariants  + length(allgenes)+ 1,
				ncol = 1,
				height = unit.c(
					unit(13,'lines'),
					unit(3,'lines'),
					unit(1,'lines'),
					rep(unit(axis_height,'lines'),totalvariants),
					rep(unit(axis_height,'lines'),length(allgenes)),
					unit(1,'null')
				)
			);

			# Create a viewport using the above layout, and positioned within the 
			# 2nd column of the outer viewport.
			pushViewport(viewport(
				layout.pos.col = 2,
				layout = inner_layout,
				name = "inner",
				gp = gpar(cex=0.6)
			));

			pushViewport(viewport(
				layout.pos.row = 2,
				xscale = xscale,
				name = "header"
			));
			
			name = "VARIANT"
			name = addspaces(name, longest_variant)
			betaheader = ''
			if(length(longest_beta >= 1)){betaheader = addspaces("BETA", longest_beta);}
			
			geneheader= addspaces("GENE",longestgenename)
			macheader = addspaces("MAC",longest_mac)
			
			testheader = ''
#			longest_test.p = NA
			for(testnum in 1:length(tests[[1]])){
				curtest = paste(tests[[1]][testnum],'.P',sep="")
				curtest = addspaces(curtest, '4.0e-07 ')
				testheader = paste(testheader, curtest, sep=" ")
			}
			
			geneheader = paste(geneheader, testheader, name, betaheader, addspaces("VAR.P", '4.0e-07'), "MAF     ", macheader, sep=" ")
			for(column in extracolumns){
				maxcol = allvariants[which.max(nchar(allvariants[,column])),column]
				print(maxcol)
				geneheader = paste(geneheader, addspaces(column,maxcol), sep=" ")
			}
			
			grid.text(
			  geneheader,
			  x = unit(1,'npc') + unit(1,'char'),
			  y =  unit(0.5,'npc'),
			  just = 0, 
			  check.overlap=FALSE,
			  gp = gpar(fontfamily="mono",cex=1)
			);

				tempx = convertUnit(unit(1,'npc')+convertUnit(unit(widthtouse,'in'), 'npc'), 'native')
				tempx = convertUnit(tempx, 'npc')
				grid.lines(
					x = unit(c(0,tempx),'npc'),
					y = unit(c(1,1), 'npc')
				)
				grid.lines(
					x = unit(c(0,tempx),'npc'),
					y = unit(c(0,0), 'npc')
				)

			grid.rect()
			upViewport(1);
			#now you are back in the inner layout viewport
			j = 0

			for(genenum in 1:length(genes)){
				thisgene = genes[genenum]
				print(thisgene)
				rowstoplot = merged[merged$'GROUP' == thisgene,]
				num_variants = length(unique(rowstoplot$'MARKER_ID'))
				markersingene = unique(rowstoplot$'MARKER_ID')
				thisgene = gsub(".*_","",thisgene, perl=T)
				genepvalue = ''

				for(testnum in 1:length(tests[[1]])){
					thistestgenepvalue = rowstoplot[1,tests[[1]][testnum]]
					thistestgenepvalue = addspaces(thistestgenepvalue, "4.0e-07 ")
					genepvalue = paste(genepvalue, thistestgenepvalue,sep=" ")
				}

				rowswithminorvariants = c()
				for (i in 1:num_variants){
					j = j + 1
					this_variant = markersingene[i];
					print(this_variant)
					markerrows = rowstoplot[rowstoplot$'MARKER_ID' == this_variant,]
					
					this_pvalue = markerrows$'PVALUE.x'[1]
					this_maf = as.numeric(markerrows$'MAF'[1])
					if(this_maf < 0.0005){this_maf = format(this_maf,scientific=TRUE);}
					else{this_maf = round(this_maf, 3);}
					this_mac = sprintf("%0.4g",as.numeric(markerrows$'MAC'[1]))
					this_beta = ""

					if ('BETA' %in% colnames(markerrows)) {
						this_beta = markerrows$'BETA'[1];
					}

					#What is the rarest genotype for this variant? If it is a het, discard it because we extract hets separately
					rare_geno = names(tail(sort(table(markerrows$'GENOTYPE'),decreasing=T),1))
					if (rare_geno == '0/1') {
						rare_geno = '-999';
					}

					if (nchar(this_pvalue) < 8) {
						chars_to_add = 8 - nchar(this_pvalue)
						for(chartoadd in 1:chars_to_add){
							this_pvalue = paste(this_pvalue, " ", sep="")
						}
					}

					if (nchar(this_maf) < 5) {
						chars_to_add = 5 - nchar(this_maf)
						for(chartoadd in 1:chars_to_add){
							this_maf = paste(this_maf, " ", sep="")
						}
					}
					
					#only print gene and gene-pvalue if it's the first variant
					if (i > 1) {
						thisgene = paste(rep(" ",nchar(thisgene)),collapse="")
						genepvalue = paste(rep(" ",nchar(genepvalue)),collapse="")
					}
					
					print(j)

					# Create the viewport. 
					pushViewport(viewport(
					layout.pos.row = 2 + j,
					xscale = xscale,
					name = paste(thisgene, this_variant,this_pvalue,this_maf, this_mac, sep=" ")
					));
					#draw marker axes boundaries
					
					# Get the phenotype values for just those people with the rare genotype and hets. 
					rare_pheno = as.numeric(markerrows[markerrows['GENOTYPE'] == rare_geno,phenotype])
					hets_pheno = as.numeric(markerrows[markerrows['GENOTYPE'] == '0/1',phenotype])
					
					notrarerows = which(markerrows['GENOTYPE'] != rare_geno) 
					nothetrows = which(markerrows['GENOTYPE'] != '0/1') 
					notmissing = which(markerrows['GENOTYPE'] != 'NA')
					
					if (length(rowswithminorvariants) == 0){
						rowswithminorvariants = markerrows[intersect(intersect(notrarerows, nothetrows),notmissing),]
					}else{
						rowswithminorvariants = rbind(rowswithminorvariants, markerrows[intersect(intersect(notrarerows, nothetrows),notmissing),])
					}
					
					rare_pheno = rare_pheno[!is.na(rare_pheno)]
					hets_pheno = hets_pheno[!is.na(hets_pheno)]			
					if (length(rare_pheno) > 0){
						meanpheno_variant = (sum(rare_pheno) + sum(hets_pheno))/(length(rare_pheno) + length(hets_pheno))
					}else{
						meanpheno_variant = (sum(hets_pheno))/(length(hets_pheno))
					}
					# Ignore the variant if it only had 1 genotype.
					if (!length(unique(markerrows$'GENOTYPE')) == 1) {
						if(length(hets_pheno) > 0){
							grid.points(
							x = unit(hets_pheno,'native'),
							y = unit(rep(0.5,length(hets_pheno)),'npc'),
							pch = 1,
							gp = gpar(
								cex = 0.3,
								col="blue"
							)
							);
						}
						if(length(rare_pheno)> 0){
							grid.points(
							x = unit(rare_pheno,'native'),
							y = unit(rep(0.5,length(rare_pheno)),'npc'),
							pch = 6,
							gp = gpar(
							cex = 0.5,
							col="red"
							)
							);
						}				
					};
					
					# Write the name of the SNP on the right hand side.
					thisgene = addspaces(thisgene, longestgenename)
					this_variant = addspaces(this_variant, longest_variant)
#					if(substr(this_beta,1,1) != '-'){ this_beta = paste('+', this_beta,sep=""); }
					this_beta = addspaces(this_beta, longest_beta)
					this_mac = addspaces(this_mac, longest_mac)
					this_maf = addspaces(this_maf, "4.0e-07")
					this_pvalue = addspaces(this_pvalue, "4.0e-07")
					infoonright = paste(thisgene, genepvalue, this_variant,this_beta, this_pvalue,this_maf,this_mac, sep=" ")
					colnum = 0
					for(column in extracolumns){
						colnum = colnum + 1
						this_col = markerrows[column]
						this_col = this_col[[1]][1]
						this_col = addspaces(this_col, lengthextracolumns[colnum])
						infoonright = paste(infoonright, this_col,sep=" ")
					}
					grid.text(
					  infoonright,
					  x = unit(1,'npc') + unit(1,'char'),
					  #y = unit(1,'npc') - unit(1,'lines'),
					  y =  unit(0.5,'npc'),
					  just = 0, 
					  check.overlap=FALSE,
					  gp = gpar(fontfamily="mono",cex=1)
					);
					grid.lines(
						x=unit(phenotypemean, 'native'),
						y = unit(c(0,1),'npc'),
						gp = gpar(col="grey")
					);
					grid.lines(
						x=unit(meanpheno_variant, 'native'),
						y = unit(c(0,1),'npc'),
						gp = gpar(col="dark green")
					);			
					grid.lines(
						x=unit(0, 'npc'),
						y = unit(c(0,1),'npc'),
						gp = gpar(col="black")
					);
					grid.lines(
						x=unit(1, 'npc'),
						y = unit(c(0,1),'npc'),
						gp = gpar(col="black")
					);

					tempx = convertUnit(unit(1,'npc')+convertUnit(unit(widthtouse,'in'), 'npc'), 'native')
					tempx = convertUnit(tempx, 'npc')
					
					if(i <= num_variants){
						grid.lines(
							x = unit(c(1,tempx),'npc'),
							y = unit(c(0,0), 'npc'),
							gp= gpar(col="grey", lty="dotdash")
						)
						grid.lines(
							x = unit(c(0,1), 'npc'),
							y = unit(c(0,0), 'npc'),
							gp = gpar(col="grey", lty = "dotdash")
						)
					}
					# Go back up so we can push another viewport within the layout (the viewport higher up in the tree.) 	
					if(i <= num_variants){
						upViewport(1);
					}
				}

				j = j+1
				#additional row for non-carriers of the gene
				pushViewport(viewport(
					layout.pos.row = 2 + j,
					xscale = xscale,
					name = paste(thisgene, "noncarriers", sep=" ")
					));
				noncarrierspheno = as.numeric(rowswithminorvariants[,phenotype])
				grid.lines(
					x=unit(phenotypemean, 'native'),
					y = unit(c(0,1),'npc'),
					gp = gpar(col="grey")
				);
				grid.lines(
					x=unit(0, 'npc'),
					y = unit(c(0,1),'npc'),
					gp = gpar(col="black")
				);
				grid.lines(
					x=unit(1, 'npc'),
					y = unit(c(0,1),'npc'),
					gp = gpar(col="black")
				);

				# Thin out the non-carrier phenotype values a little to avoid overplotting (and making the PDF take
				# a long time to render.)
				noncarrierspheno = unique(round(noncarrierspheno,ROUND_NONCARRIER))
				grid.points(
                    x = unit(noncarrierspheno,'native'),
                    y = unit(rep(0.5,length(noncarrierspheno)),'npc'),
                    pch = 1,
                    gp = gpar(
                        cex = 0.2,
                        col="grey",
                        pch=3
                    )
				);
				
				#extend the gridline on the right
				tempx = convertUnit(unit(1,'npc')+convertUnit(unit(widthtouse,'in'), 'npc'), 'native')
				tempx = convertUnit(tempx, 'npc')
				grid.lines(
					
					x = unit(c(0,tempx),'npc'),
					y = unit(c(0,0), 'npc')
				)
				
				
				upViewport(1);
			}
			
			print("here")
			print(current.viewport())
			# Now draw the phenotype distribution. 
			pushViewport(viewport(
			layout.pos.row = 1,
			xscale = xscale,
			yscale = yscale,
			name = "hist"
			))
			print(current.viewport())

			# Now draw the histogram. 
			panel.histogram(phenotyperesiduals,breaks=breaks,type="percent");
			grid.rect();
			grid.xaxis();
			grid.yaxis();
			grid.text("Percent of Total",x = unit(-4,'lines'),rot = 90);
			if(length(covariates) > 0){adjustedpheno = paste("Adjusted ", phenotype,sep="");}
			else{adjustedpheno = phenotype}
			grid.text(paste(adjustedpheno, ": ", length(phenotypes), "samples", sep=" "), x=unit(16, 'char'), y=unit(0.97, 'npc'))
			panel.rug(x=max(phenotyperesiduals),y=0)
			panel.rug(x=min(phenotyperesiduals),y=0)
			upViewport(1);
			
			###THIS PLOT IS DONE. NOW START PREPARING INPUT FOR NEXT PLOT
			if(numgenesplotted < length(allgenes)){
				genes = allgenes[(numgenesplotted+1):length(allgenes)]
				numgenestoinclude = 0
				totalvariants = 0
				length_longest_variant = 0
				longest_variant = ""
				longest_beta = ""
				for(genenum in 1:length(genes)){
					gene = genes[genenum]
					rowstoplot = merged[merged$'GROUP' == gene,]
					num_variants = length(unique(rowstoplot$'MARKER_ID'))
					totalvariants = totalvariants + num_variants
					markersingene = unique(rowstoplot$'MARKER_ID')
					if((numgenestoinclude == 0) || (totalvariants <= 40)){
						numgenestoinclude = numgenestoinclude + 1;
						l = nchar(markersingene[which.max(nchar(markersingene))])
						if (l > length_longest_variant){
							length_longest_variant = l;
							longest_variant = markersingene[which.max(nchar(markersingene))];
						}
						l = rowstoplot$'BETA'[which.max(nchar(as.numeric(rowstoplot$'BETA'),))]
						if(nchar(l) > nchar(longest_beta)){longest_beta = l;}
						l = rowstoplot$'MAC'[which.max(nchar(rowstoplot$'MAC'))]
						if(nchar(l) > nchar(longest_mac)){longest_mac = l;}
					}
				}
				if(numgenestoinclude >= 1){totalvariants = min(totalvariants, 40);}
				numiter  = numiter + 1
				
			}
			else{
				donothing = 1
				}
		}
		dev.off()
}

#if the phenotype is categorical, then there is no need to create a PDF. 
if(categorical == 1){
	pheno1 = unique(phenotypes)[1]
	pheno2 = unique(phenotypes)[2]
	tabletowrite = matrix(ncol =10  , nrow = length(genes)*length(unique(merged$'MARKER')))
	
	j = 0
	for(genenum in 1:length(genes)){
		thisgene = genes[genenum]
		rowstouse = merged[merged$'GROUP' == thisgene,]
		num_variants = length(unique(rowstouse$'MARKER_ID'))
		markersingene = unique(rowstouse$'MARKER_ID')
		thisgene = gsub(".*_","",thisgene, perl=T)
		genepvalue = ''
		for(testnum in 1:length(tests)){
			thistestgenepvalue = sprintf("%0.3g", rowstoplot$'PVALUE.y'[1])
			genepvalue = paste(genepvalue, thistestgenepvalue,sep="\t")
		}
		for (i in 1:num_variants) {
			j = j + 1
			this_variant = markersingene[i];
			markerrows = rowstoplot[rowstouse$'MARKER_ID' == this_variant,]
			this_pvalue = sprintf("%0.3g", markerrows$'PVALUE.x'[1])
			this_maf = round(markerrows$'MAF'[1],digits=3)
			this_mac = markerrows$'MAC'[1]
			this_beta = ""
			if('BETA' %in% colnames(markerrows)){this_beta = round(markerrows$'BETA'[1], digits=2);}
			#What is the rarest genotype for this variant? If it is a het, discard it because we extract hets separately
			rare_geno = names(tail(sort(table(markerrows$'GENOTYPE'),decreasing=T),1))
			if(rare_geno == '0/1'){
				rare_geno = '-999';
			}
		
			case_rare = markerrows[markerrows['GENOTYPE'] == rare_geno, phenotype]
			case_rare = length(case_rare[which(case_rare == pheno1)])
			
			case2_rare = markerrows[markerrows['GENOTYPE'] == rare_geno, phenotype]
			case2_rare = length(case2_rare[which(case2_rare == pheno2)])
			
			case_het = markerrows[markerrows['GENOTYPE'] == '0/1', phenotype]
			case_het = length(case_het[which(case_het == pheno1)])
			
			case2_het = markerrows[markerrows['GENOTYPE'] == '0/1', phenotype]
			case2_het = length(case2_het[which(case2_het == pheno2)])
			
			tabletowrite[j,1] = thisgene
			tabletowrite[j,2] = genepvalue
			tabletowrite[j,3] = this_variant
			tabletowrite[j,4] = this_pvalue
			tabletowrite[j,5] = this_maf
			tabletowrite[j,6] = this_beta
			tabletowrite[j,7] = case_rare
			tabletowrite[j,8] = case_het
			tabletowrite[j,9] = case2_rare
			tabletowrite[j,10] = case2_het
		}
		write.table(tabletowrite, file=paste(prefix, ".topgenes.txt",sep=''), quote=FALSE, col.names=F, row.names=F)
	}
}
