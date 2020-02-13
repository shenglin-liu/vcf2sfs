# Functions for creating site frequency spectra (SFS) from a VCF file.
# Functions for manipulating and plotting SFS.
# By Shenglin Liu, Feb 12, 2020.

## Requirements:
# Only for diploids.
# Convention for file name extension (not mandatory):
#	fs (dadi), txt (dadi SNP data), fsc (fastSimCoal).
# Format of dadi:
#	strictly two lines;
#	dimension vector for the first line;
#	SFS for the second line;
#	no information for fold, mask or population names.
# Format of dadi SNP data:
#	refIn	refOut	Allele1	Aan	Aro	Allele2	Aan	Aro
#	AAA	AAA	A	49	56	T	7	0

## Frequent variables in the script:
# gt: a genotype object created from a VCF file and a popmap file.
# sfs: a table or array object representing SFS.

## Function list:
#vcf2gt(f.vcf,f.popmap)
#choose.pops(gt,pops)
#samSize(gt)
#minSamSize(gt)
#viewMissing(gt)
#filter.gt(gt,filter.indi=NA,filter.snp=NA)
#gt2sfs.raw(gt,pops)
#gt2sfs.impute(gt,pops,sampleSizes)
#project.sfs(sfs,sampleSizes)
#write.sfs.dadi(sfs,f.output)
#read.sfs.dadi(f.dadi,popnames=NULL)
#fold.sfs(sfs)
#sample.sfs(sfs)
#gt2snp(gt,pops)
#snp2sfs.raw(snp,pops,sampleSizes=NA)
#plot.sfs(sfs,cr=c(1,20000),colScheme=function(x)rainbow(x,s=0.8),...)
#plot.2D.pairwise(pops,sep="&",ext.o="fs",ext.e=NA,fold=T,mask=NA,
#	cr=c(1,20000),colScheme=function(x)rainbow(x,s=0.8))
#write.1D.fsc(sfs,f.output)
#write.2D.fsc(sfs,f.output)
#read.1D.fsc(f.fsc,popname=NULL)
#read.2D.fsc(f.fsc,popnames=NULL)

## Read VCF file and popmap file, and create a gt (genotype) object.
# It will generate a named list of two elements, i.e., a popmap vector and a genotype matrix.
# The FORMAT column of the VCF file should start with GT tag;
#	only biallelic, encoded with 0 and 1; "." as missing.
# The popmap file lists the population IDs of individuals;
#	two columns; tab delimited; individual ID as 1st column, pop ID as 2nd column;
#	1st column can be arbitrary; 2nd column should be an integer or string vector;
#	minimum sample size per population is 2.
vcf2gt<-function(f.vcf,f.popmap)
{
	oldw<-getOption("warn")
	options(warn=-1)
	
	# Read VCF file and popmap file.
	vcf.gt<-as.matrix(read.table(f.vcf,sep="\t",stringsAsFactors=F)[,-c(1:9)])
	popmap<-read.table(f.popmap,sep="\t",stringsAsFactors=F)[,2]
	
	nrow.vcf<-nrow(vcf.gt)
	ncol.vcf<-ncol(vcf.gt)
	
	# Parse genotypes.
	chrom1<-substring(vcf.gt,1,1)
	chrom2<-substring(vcf.gt,3,3)
	chrom<-matrix(as.integer(chrom1)+as.integer(chrom2),nrow.vcf,ncol.vcf)
	
	options(warn=oldw)
	
	list(popmap=popmap,genotype=chrom)
}

## Choose populations (subsetting the gt object by popultions).
# pops: a character or integer vector; IDs of the chosen populations.
choose.pops<-function(gt,pops)
{
	popmap<-gt$popmap
	chrom<-gt$genotype
	
	# Choose populations.
	index<-sapply(pops,function(x){which(popmap==x)})
	index<-sort(unlist(index))
	popmap<-popmap[index]
	chrom<-chrom[,index]
	
	list(popmap=popmap,genotype=chrom)
}

## Calculate sample sizes of the populations.
samSize<-function(gt)
{
	popmap<-gt$popmap
	chrom<-gt$genotype
	
	pops<-unique(popmap)
	sampleSizes<-sapply(pops,function(x){sum(popmap==x)})
	names(sampleSizes)<-pops
	sampleSizes
}

## Calculate minimum sample sizes of the populations accounting for the missing values.
minSamSize<-function(gt)
{
	popmap<-gt$popmap
	chrom<-gt$genotype
	
	pops<-unique(popmap)
	sampleSizes<-sapply(pops,function(x){min(rowSums(!is.na(chrom[,popmap==x])))})
	names(sampleSizes)<-pops
	sampleSizes
}

## View the distribution of the missing values.
# Help decide whether some individuals or SNPs should be filtered out.
viewMissing<-function(gt)
{
	popmap<-gt$popmap
	chrom<-gt$genotype
	
	missing<-is.na(chrom)
	layout(matrix(c(1,2,1,3),2,2))
	image(1:nrow(missing),1:ncol(missing),missing,
		xlab=paste("SNPs: ",dim(missing)[1],sep=""),
		ylab=paste("Individuals: ",dim(missing)[2],sep=""),
		col=c("olivedrab3","firebrick3"),main="Missing values (red)")
	plot(sort(colSums(missing)),main="By individual",
		xlab="Individuals",ylab="Number of missing values")
	plot(table(rowSums(missing)),main="By SNP",
		xlab="Number of missing values",ylab="Number of SNPs")
}

## Filter the gt object for the missing values.
# filter.indi: integer; when the number of missing values is lager than this, 
#	the individual will be filtered out; when NA, no filtering.
# filter.snp: integer; when the number of missing values is lager than this, 
#	the SNP will be filtered out; when NA, no filtering.
filter.gt<-function(gt,filter.indi=NA,filter.snp=NA)
{
	popmap<-gt$popmap
	chrom<-gt$genotype
	
	# Filter individuals according to missing values.
	if(!is.na(filter.indi))
	{
		missingVal<-is.na(chrom)
		delete.indi<-which(colSums(missingVal)>filter.indi)
		if(length(delete.indi)>0)
		{
			chrom<-chrom[,-delete.indi]
			popmap<-popmap[-delete.indi]
		}
		rm(missingVal)
	}
	
	# Filter SNPs according to missing values.
	if(!is.na(filter.snp))
	{
		missingVal<-is.na(chrom)
		delete.snp<-which(rowSums(missingVal)>filter.snp)
		if(length(delete.snp)>0)
		{
			chrom<-chrom[-delete.snp,]
		}
		rm(missingVal)
	}
	
	list(popmap=popmap,genotype=chrom)
}

## Generate a SFS (table object) from the gt object.
# It will output a SFS based on raw count without accounting for the missing values.
# pops: a character or integer vector; IDs of populations to be included in the SFS.
gt2sfs.raw<-function(gt,pops)
{
	popmap<-gt$popmap
	chrom<-gt$genotype
	
	nrow.vcf<-nrow(chrom)
	ncol.vcf<-ncol(chrom)
	
	n.pop<-length(pops)
	
	# Number of chromosomes.
	ns.chr<-sapply(pops,function(x){sum(popmap==x)})*2
	
	# SFS based on raw count.
	cnt<-matrix(0,nrow.vcf,n.pop)
	ext<-list()
	for(i in 1:n.pop)
	{
		index<-which(popmap==pops[i])
		cnt[,i]<-rowSums(chrom[,index],na.rm=T)
		ext<-c(ext,list(0:ns.chr[i]))
	}
	ext<-as.matrix(expand.grid(ext))
	cnt<-data.frame(rbind(cnt,ext))
	sfs.raw<-table(cnt)-1
	names(dimnames(sfs.raw))<-pops
	
	sfs.raw
}

## Generate a SFS (table object) from the gt object.
# It will output an imputed SFS accounting for the missing values through binomial distribution.
# pops: a character or integer vector; IDs of populations to be included in the SFS.
# sampleSizes: an integer vector; downsized number of individuals in each population.
gt2sfs.impute<-function(gt,pops,sampleSizes)
{
	popmap<-gt$popmap
	chrom<-gt$genotype/2
	
	nrow.vcf<-nrow(chrom)
	ncol.vcf<-ncol(chrom)
	
	n.pop<-length(pops)
	
	# Number of chromosomes.
	ns.chr<-sampleSizes*2
	
	# SFS based on count imputed through binomial distribution.
	cnt<-matrix(0,nrow.vcf,n.pop)
	ext<-list()
	for(i in 1:n.pop)
	{
		index<-which(popmap==pops[i])
		p<-rowMeans(chrom[,index],na.rm=T)
		cnt[,i]<-rbinom(nrow.vcf,ns.chr[i],p)
		ext<-c(ext,list(0:ns.chr[i]))
	}
	ext<-as.matrix(expand.grid(ext))
	cnt<-data.frame(rbind(cnt,ext))
	sfs.imp<-table(cnt)-1
	names(dimnames(sfs.imp))<-pops
	
	sfs.imp
}

## Project a SFS to smaller sample sizes according to hypergeometric distribution.
# sampleSizes: an integer vector; downsized number of individuals in each population.
project.sfs<-function(sfs,sampleSizes)
{
	project<-function(a,m,n)
	{
		b<-numeric(m+1)
		for(i in 0:m)for(j in i:(n-m+i))
		{
			b[i+1]<-b[i+1]+choose(m,i)*choose(n-m,j-i)/choose(n,j)*a[j+1]
		}
		b
	}
	dims<-dim(sfs)
	l.dim<-length(dims)
	sfs.pro<-sfs
	if(l.dim==1)
	{
		m<-sampleSizes[1]*2
		n<-dims[1]-1
		sfs.pro<-array(project(sfs.pro,m,n))
	}else
	{
		for(k in 1:l.dim)
		{
			m<-sampleSizes[k]*2
			n<-dims[k]-1
			sfs.pro<-apply(sfs.pro,c(1:l.dim)[-k],function(x)project(x,m,n))
			perm<-integer(l.dim)
			perm[k]<-1
			perm[-k]<-2:l.dim
			sfs.pro<-aperm(sfs.pro,perm)
		}
	}
	dimnames(sfs.pro)<-lapply(sampleSizes*2,function(x)0:x)
	names(dimnames(sfs.pro))<-names(dimnames(sfs))
	sfs.pro<-as.table(sfs.pro)
	sfs.pro
}

## Write a SFS to a file in dadi format.
write.sfs.dadi<-function(sfs,f.output)
{
	sink(f.output)
	cat(dim(sfs))
	cat("\n")
	cat(aperm(sfs))
	sink()
}

## Read a SFS from a dadi format file.
read.sfs.dadi<-function(f.dadi,popnames=NULL)
{
	dims<-scan(f.dadi,nlines=1,comment.char="#",quiet=T)
	sfs<-scan(f.dadi,skip=1,nlines=1,comment.char="#",quiet=T)
	sfs<-aperm(array(sfs,rev(dims)))
	sfs<-as.table(sfs)
	dimnames(sfs)<-lapply(dims-1,function(x)0:x)
	names(dimnames(sfs))<-popnames
	sfs
}

## Fold a SFS.
fold.sfs<-function(sfs)
{
	sfs[]<-sfs+rev(sfs)
	dims<-dim(sfs)
	cnt.pool<-rowSums(expand.grid(lapply(dims-1,function(x)0:x)))
	index<-cnt.pool>(sum(dims-1)/2)
	sfs[index]<-0
	index<-cnt.pool==(sum(dims-1)/2)
	sfs[index]<-sfs[index]/2
	sfs
}

## Resampling from a SFS for the purpose of bootstrapping.
sample.sfs<-function(sfs)
{
	n.site<-sum(sfs)
	sfs[1]<-0
	sfs[]<-rpois(length(sfs),sfs)
	sfs[1]<-n.site-sum(sfs)
	sfs
}

## Transform the gt object to a headered data.frame of dadi SNP data format.
# The function will output a headered data.frame as such:
#	refIn	refOut	Allele1	Aan	Aro	Allele2	Aan	Aro
#	AAA	AAA	A	49	56	T	7	0
# pops: a character or integer vector; IDs of populations to be included in the SFS.
gt2snp<-function(gt,pops)
{
	popmap<-gt$popmap
	chrom<-gt$genotype
	
	nrow.vcf<-nrow(chrom)
	ncol.vcf<-ncol(chrom)
	
	n.pop<-length(pops)
	
	ref<-integer(0)
	alt<-integer(0)
	for(i in 1:n.pop)
	{
		index<-which(popmap==pops[i])
		temp<-rowSums(chrom[,index],na.rm=T)
		ref<-cbind(ref,rowSums(!is.na(chrom[,index]))*2-temp)
		alt<-cbind(alt,temp)
	}
	snp<-data.frame("AAA","AAA","A",ref,"T",alt)
	names(snp)<-c("refIn","refOut","Allele1",pops,"Allele2",pops)
	
	snp
}

## Transform dadi SNP data format to SFS (table object).
# It will output a SFS based on raw count.
# snp: a headered data.frame as such:
#	refIn	refOut	Allele1	Aan	Aro	Allele2	Aan	Aro
#	AAA	AAA	A	49	56	T	7	0
# pops: a character vector; IDs of populations to be included in the SFS.
# sampleSizes: an integer verctor of the same length as pops;
#	actual sample sizes (number of diploid individuals) of populations;
#	if NA, calculate from the SNP data.
snp2sfs.raw<-function(snp,pops,sampleSizes=NA)
{
	i.ref<-sapply(pops,function(x)which(names(snp)==x)[1])
	i.alt<-sapply(pops,function(x)which(names(snp)==x)[2])
	names(i.ref)<-names(i.alt)<-pops
	
	# Number of chromosomes.
	if(is.na(sampleSizes[1]))
	{
		ns.chr<-sapply(pops,function(x)max(snp[,i.ref[x]]+snp[,i.alt[x]]))
	}else
	{
		ns.chr<-sampleSizes*2
	}
	
	# SFS based on raw count.
	ext<-list()
	for(n.chr in ns.chr)ext<-c(ext,list(0:n.chr))
	ext<-as.matrix(expand.grid(ext))
	cnt<-data.frame(rbind(as.matrix(snp[i.alt]),ext))
	sfs.raw<-table(cnt)-1
	names(dimnames(sfs.raw))<-pops
	
	sfs.raw
}

## Plot a SFS (barplot for 1D, image for 2D).
# cr: color range.
plot.sfs<-function(sfs,cr=c(1,20000),colScheme=function(x)rainbow(x,s=0.8),...)
{
	pops<-names(dimnames(sfs))
	dims<-dim(sfs)
	n.pop<-length(dims)
	sfs[1]<-sfs[length(sfs)]<-0
	
	# Color coding for number of SNPs.
	cr<-log10(cr)			## Color range
	nb<-101				## Number of breaks
	breaks<-10^seq(cr[1],cr[2],length.out=nb)
	
	if(n.pop==1)
	{
		barplot(sfs,names.arg=c(0,rep(NA,dims-2),dims-1),main=pops,...)
	}
	if(n.pop==2)
	{
		image(sfs,axes=F,col=colScheme(nb-1),breaks=breaks,
			xlab=pops[1],ylab=pops[2],...)
		axis(1,at=c(0,1),labels=c(0,dims[1]-1))
		axis(2,at=c(0,1),labels=c(0,dims[2]-1))
		box()
	}
	if(n.pop>2)stop("Cannot plot SFS with dimension higher than 2!")
}

## Plot pairwise 2D-SFSs among populations.
# Reads in from SFS files in dadi format.
# Each file should have names of population pairs separated by certain symbol.
# pops: names of the populations as in the file names; order matters.
# sep: the separation symbol in the file names.
# ext.o: extension name of the files for the observed SFSs.
# ext.e: extension name of the files for the expected SFSs; if NA, not plotted.
# mask: a two-column matrix; coordinates of the bins to be masked; if NA, only mask fixed.
# cr: color range.
# colScheme: color scheme function.
plot.2D.pairwise<-function(pops,sep="&",ext.o="fs",ext.e=NA,fold=T,mask=NA,
	cr=c(1,20000),colScheme=function(x)rainbow(x,s=0.8))
{
	n<-length(pops)
	
	# Color coding for number of SNPs.
	cr<-log10(cr)			## Color range
	nb<-101				## Number of breaks
	breaks<-10^seq(cr[1],cr[2],length.out=nb)
	
	# Layout of the plotting region.
	a<-matrix(0,n,n);k<-1;for(i in 1:(n-1))for(j in (i+1):n){a[i,j]<-k;k<-k+1}
	temp<-a*2-1;temp[temp[,]==-1]<-0;a<-temp+t(a)*2;a[a[,]==0]<-n*(n-1)+1:n
	a<-cbind(a,n*n+1)
	layout(a,widths=c(rep(8/9/n,n),1/9))
	
	# Plot SFSs.
	par(mar=c(0,0.5,0.5,0),oma=c(1,1,1,1))	##
	for(i in 1:(n-1))
	{
		for(j in (i+1):n)
		{
			# Plot observed.
			f.dadi<-sprintf("%s%s%s.%s",pops[i],sep,pops[j],ext.o)
			sfs<-read.sfs.dadi(f.dadi)
			if(fold)sfs<-fold.sfs(sfs)
			sfs[1]<-sfs[length(sfs)]<-0
			if(!is.na(mask[1]))sfs[mask]<-0
			image(t(sfs),col=colScheme(nb-1),axes=F,breaks=breaks)
			temp<-sum(sfs)
			
			# Plot expected.
			if(is.na(ext.e))
			{
				plot(0,type="n",axes=F,xlab=NA,ylab=NA)
			}else
			{
				f.dadi<-sprintf("%s%s%s.%s",pops[i],sep,pops[j],ext.e)
				sfs<-read.sfs.dadi(f.dadi)
				sfs[sfs>2]<-0
				if(fold)sfs<-fold.sfs(sfs)
				sfs[1]<-sfs[length(sfs)]<-0
				if(!is.na(mask[1]))sfs[mask]<-0
				sfs<-round(sfs/sum(sfs,na.rm=T)*temp)
				image(sfs[nrow(sfs):1,ncol(sfs):1],
					col=colScheme(nb-1),axes=F,breaks=breaks)
			}
		}
	}
	
	# Plot diagonal.
	for(i in 1:n)
	{
		plot(0,0,type="n",axes=F,xlab=NA,ylab=NA)
		text(0,0,pops[i],cex=1.5,font=2)
	}
	
	# Color bar.
	h<-dev.size()[2]
	par(mar=c(12,2,12,3)/9*h,cex.axis=1.5)	##
	image(1,breaks[-nb],matrix(1:(nb-1),1),
		col=colScheme(nb-1),log="y",axes=F,xlab=NA,ylab=NA)
	axis(4)
}

## Write a 1D-SFS to a file in fastSimCoal format.
write.1D.fsc<-function(sfs,f.output)
{
	fsc<-matrix(sfs,1)
	pop<-names(dimnames(sfs))
	# add column names
	colnames(fsc)<-paste("d",pop,"_",1:ncol(fsc)-1,sep="")
	# output
	cat("1 observations\n",file=f.output)
	oldw<-getOption("warn")
	options(warn=-1)
	write.table(fsc,file=f.output,append=T,sep="\t",col.names=T,row.names=F,quote=F)
	options(warn=oldw)
}

## Write a 2D-SFS to a file in fastSimCoal format.
write.2D.fsc<-function(sfs,f.output)
{
	fsc<-as.matrix(sfs)
	pops<-names(dimnames(sfs))
	# add row names and column names
	rownames(fsc)<-paste("d",pops[1],"_",1:nrow(fsc)-1,sep="")
	colnames(fsc)<-paste("d",pops[2],"_",1:ncol(fsc)-1,sep="")
	# output
	cat("1 observations\n\t",file=f.output)
	oldw<-getOption("warn")
	options(warn=-1)
	write.table(fsc,file=f.output,append=T,sep="\t",col.names=T,row.names=T,quote=F)
	options(warn=oldw)
}

## Read a 1D-SFS from a fastSimCoal format file.
read.1D.fsc<-function(f.fsc,popname=NULL)
{
	fsc<-read.table(f.fsc,header=T,sep="\t",skip=1)
	sfs<-array(unlist(fsc),ncol(fsc))
	sfs<-as.table(sfs)
	dimnames(sfs)<-lapply(ncol(fsc)-1,function(x)0:x)
	names(dimnames(sfs))<-popname
	sfs
}

## Read a 2D-SFS from a fastSimCoal format file.
read.2D.fsc<-function(f.fsc,popnames=NULL)
{
	fsc<-read.table(f.fsc,header=T,sep="\t",skip=1,row.names=1)
	sfs<-array(unlist(fsc),dim(fsc))
	sfs<-as.table(sfs)
	dimnames(sfs)<-lapply(dim(fsc)-1,function(x)0:x)
	names(dimnames(sfs))<-popnames
	sfs
}

