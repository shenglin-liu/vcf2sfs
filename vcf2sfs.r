# Functions for creating site frequency spectra from a VCF file.
# Modified from "vcf2fs.r".
# By Shenglin Liu, Apr 13, 2016.


## Check the distribution of the missing values.
# Help decide whether some individuals or SNPs should be filtered out.
checkMissing<-function(f.vcf,f.popmap,pops)
{
	vcf<-as.matrix(read.table(f.vcf,sep="\t",stringsAsFactors=F)[,-c(3:9)])
	snpid<-paste(vcf[,1],as.integer(vcf[,2]),sep="_")
	vcf<-vcf[,-c(1:2)]
	
	popmap<-read.table(f.popmap,sep="\t",stringsAsFactors=F)[,2]
	index<-sapply(pops,function(x){which(popmap==x)})
	index<-sort(unlist(index))
	popmap<-popmap[index]
	vcf<-vcf[,index]
	
	nrow.vcf<-nrow(vcf)
	ncol.vcf<-ncol(vcf)
	chrom1<-matrix(substring(vcf,1,1),nrow.vcf,ncol.vcf)
	chrom2<-matrix(substring(vcf,3,3),nrow.vcf,ncol.vcf)
	
	missing<-chrom1=="."
	par(mfrow=c(1,3))
	image(missing,
		xlab=paste("SNPs: ",dim(missing)[1],sep=""),
		ylab=paste("Individuals: ",dim(missing)[2],sep=""))
	plot(sort(colSums(missing)),main="By individual",
		xlab="Individuals",ylab="Number of missing values")
	plot(table(rowSums(missing)),main="By SNP",
		xlab="Number of missing values",ylab="Number of SNPs")
}

## Transform VCF data to SFS (dadi format)
# The SNPs must be biallelic.
# f.popmap: popmap file; integers for the second column, pop IDs.
# pops: integers; IDs (same as in popmap) of pops to be included in the SFS.
# filter.indi: integer; when the number of missing values is lager than this, 
#	the individual will be filtered out; when NA, no filtering.
# filter.snp: integer; when the number of missing values is lager than this, 
#	the SNP will be filtered out; when NA, no filtering.
# n.digit: number of digits to keep when rounding up the allele frequencies; 
#	VERY IMPORTANT!!!
vcf2dadi<-function(f.vcf, f.popmap, f.output, pops,
		ploidy=2, n.digit=4,
		filter.indi=NA, filter.snp=NA)
{
	vcf<-as.matrix(read.table(f.vcf,sep="\t",stringsAsFactors=F)[,-c(3:9)])
	snpid<-paste(vcf[,1],as.integer(vcf[,2]),sep="_")
	vcf<-vcf[,-c(1:2)]
	
	# Choose populations
	popmap<-read.table(f.popmap,sep="\t",stringsAsFactors=F)[,2]
	index<-sapply(pops,function(x){which(popmap==x)})
	index<-sort(unlist(index))
	popmap<-popmap[index]
	vcf<-vcf[,index]
	
	# Parse genotypes
	nrow.vcf<-nrow(vcf)
	ncol.vcf<-ncol(vcf)
	chrom1<-matrix(substring(vcf,1,1),nrow.vcf,ncol.vcf)
	chrom2<-matrix(substring(vcf,3,3),nrow.vcf,ncol.vcf)
	
	# Filter individuals according to missing values
	if(!is.na(filter.indi))
	{
		missingVal<-chrom1=="."
		delete.indi<-which(colSums(missingVal)>filter.indi)
		if(length(delete.indi)>0)
		{
			chrom1<-chrom1[,-delete.indi]
			chrom2<-chrom2[,-delete.indi]
			popmap<-popmap[-delete.indi]
		}
	}
	
	# Filter SNPs according to missing values
	if(!is.na(filter.snp))
	{
		missingVal<-chrom1=="."
		delete.snp<-which(rowSums(missingVal)>filter.snp)
		if(length(delete.snp)>0)
		{
			chrom1<-chrom1[-delete.snp,]
			chrom2<-chrom2[-delete.snp,]
			snpid<-snpid[-delete.snp]
		}
	}
	
	nrow.vcf<-nrow(chrom1)
	ncol.vcf<-ncol(chrom1)
	
	n.pop<-length(pops)
	
	# Calculate alternative allele frequencies
	freq<-matrix(0,nrow.vcf,n.pop)
	colnames(freq)<-paste("Pop",pops,sep="_")
	rownames(freq)<-snpid
	for(i in 1:n.pop)
	{
		index<-which(popmap==pops[i])
		sub.chrom1<-chrom1[,index]
		sub.chrom2<-chrom2[,index]
		ref<-rowSums(sub.chrom1=="0")+rowSums(sub.chrom2=="0")
		alt<-rowSums(sub.chrom1=="1")+rowSums(sub.chrom2=="1")
		freq[,i]<-alt/(ref+alt)
	}
	
	# creating a list containing the allele frequencies of the target populations
	data<-vector("list",n.pop)
	for(i in 1:n.pop)
	{
		data[[i]]<-round(as.numeric(freq[,i]),n.digit)	## n.digit!!! Assigning allele frequencies for the ith population
		names(data[[i]])<-snpid	# the names are very important!!
	}
	names(data)<-pops	# these names are also very important; have to be numerics
	
	sampleSizes<-sapply(pops,function(x){sum(popmap==x)})
	
	write.jafs(data,f.output,sampleSizes,ploidy,n.digit)
}

## Generate a joint allele frequency spectrum (dadi format) from a list.
# data: list object; the name of each element is a numeric value as population ID;
#	each element is a numeric vector containing the alternative allele frequencies (p) of a population;
#	the name of each element in the vector is SNP ID.
# file: character value; output file name.
# sampleSizes: numeric vector; the number of individuals in each population.
# ploidy: numeric value; the ploidy level of the individuals.
# n.digit: number of effective digits to keep after decimal; very important!!!
write.jafs<-function(data,file,sampleSizes,ploidy,n.digit)
{
	pops=as.numeric(names(data))
	n.pop=length(pops)
	n.bin=sampleSizes*ploidy+1
	
	zz<-file(file,"w")
	writeLines(as.character(n.bin),zz,sep=" ")
	writeLines("",zz)
	
	recur<-function(data,pops,n.bin,n.digit)
	{
		n.pop=length(pops)
		n.bin1=n.bin[1]
		bound=round(c(-1,seq(0,1,1/(n.bin1-1))),n.digit)	## n.digit!!! Setting boundaries for the bins
		temp=data[[1]]
		if(length(pops)>1)
		{
			for(i.bin in 1:n.bin1)
			{
				index=names(temp)[(temp<=bound[i.bin+1])==(temp>bound[i.bin])]
				newdata=vector("list",n.pop-1)
				for(i in 2:n.pop) newdata[[i-1]]=na.omit(data[[i]][index])
				recur(newdata,pops[-1],n.bin[-1],n.digit)
			}
		} else
		{
			for(i.bin in 1:n.bin1)
			{
				writeLines(as.character(sum((temp<=bound[i.bin+1])==(temp>bound[i.bin]))),zz,sep=" ")
			}
		}
	}
	
	recur(data,pops,n.bin,n.digit)
	close(zz)
}

## Change from dadi format to fastSimCoal format, and fold
# Only for fastSimCoal (one dimension)
# A plot will be generated as well.
# pop: ID for the population; doesn't have to be the same as in popmap
dadi2fsc.1pop<-function(f.fs,f.output,pop,fold=T)
{
	dimention<-as.integer(read.table(f.fs,nrows=1))
	fs<-scan(f.fs,skip=1,quiet=T)
	
	# format into matrix
	fs<-matrix(fs,1,dimention[1])
	
	# fold (with overall minor allele frequency)
	if(fold)
	{
		fs.turn<-fs[,dimention[1]:1]
		fs<-fs+fs.turn
		index<-(dimention[1]+1)/2
		fs[,index]<-fs[,index]/2
		index<-1:index
		fs[,-index]<-0
	}
	
	# add row names and column names
	colnames(fs)<-paste("d",pop[1],"_",1:dimention[1]-1,sep="")
	
	# output
	cat("1 observations\n",file=f.output)
	oldw<-getOption("warn")
	options(warn=-1)
	write.table(fs,file=f.output,append=T,sep="\t",col.names=T,row.names=F,quote=F)
	options(warn=oldw)
	barplot(as.integer(fs))
}

## Change from dadi format to fastSimCoal format, and fold.
# Only for fastSimCoal (two dimensions).
# A plot will be generated as well.
# pops: ID's for the populations; don't have to be the same as in popmap.
dadi2fsc.2pop<-function(f.fs,f.output,pops,fold=T)
{
	dimention<-as.integer(read.table(f.fs,nrows=1))
	fs<-scan(f.fs,skip=1,quiet=T)
	
	# format into matrix (row for pop1, col for pop2)
	fs<-t(matrix(fs,dimention[2],dimention[1]))
	
	# fold (with overall minor allele frequency)
	if(fold)
	{
		half<-sum(dimention)/2
		index<-t(matrix(1:dimention[2],dimention[2],dimention[1]))
		index<-index<c(half:(half-dimention[1]+1))
		index.turn<-index[dimention[1]:1,dimention[2]:1]
		fs.turn<-fs[dimention[1]:1,dimention[2]:1]
		fs[index]<-fs[index]+fs.turn[index]
		fs[index.turn]<-0
	}
	fs[1,1]<-0
	
	# add row names and column names
	rownames(fs)<-paste("d",pops[1],"_",1:dimention[1]-1,sep="")
	colnames(fs)<-paste("d",pops[2],"_",1:dimention[2]-1,sep="")
	
	# output
	cat("1 observations\n\t",file=f.output)
	oldw<-getOption("warn")
	options(warn=-1)
	write.table(fs,file=f.output,append=T,sep="\t",col.names=T,row.names=T,quote=F)
	options(warn=oldw)
	image(log(t(fs[dimention[1]:1,])),col=rainbow(100),axes=F)
	axis(2,at=c(1,0),labels=c(0,dimention[1]-1))
	axis(3,at=c(0,1),labels=c(0,dimention[2]-1))
	mtext(paste("d",pops[2],sep=""),3,cex=1.5,line=1)
	mtext(paste("d",pops[1],sep=""),2,cex=1.5,line=1)
	box()
}

## Resampling from a SFS with a fastSimCoal format for the purpose of bootstrapping
# f.fsc: input file (fastSimCoal format)
# n.rep: number of repeats for the bootstrap
# nameroot: nameroot for the output files
fsc.sample<-function(f.fsc,n.rep,nameroot)
{
	fsc<-read.table(f.fsc,header=T,sep="\t",skip=1,row.names=1)
	n.site<-sum(fsc)
	n.row<-nrow(fsc)
	n.col<-ncol(fsc)
	n.entry<-n.row*n.col
	
	lamdas<-unlist(fsc)[-1]
	
	for(i in 1:n.rep)
	{
		temp<-rpois(n.entry-1,lamdas)
		temp<-c(n.site-sum(temp),temp)
		temp<-matrix(temp,n.row,n.col)
		rownames(temp)<-rownames(fsc)
		colnames(temp)<-colnames(fsc)
		f.output<-paste("%s-%0",nchar(n.rep),"d.fsc",sep="")
		f.output<-sprintf(f.output,nameroot,i)
		cat("1 observations\n\t",file=f.output)
		oldw<-getOption("warn")
		options(warn=-1)
		write.table(temp,file=f.output,append=T,
			sep="\t",col.names=T,row.names=T,quote=F)
		options(warn=oldw)
	}
}

