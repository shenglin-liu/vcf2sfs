# "vcf2sfs.r" contains a set of functions that can create SFS files from a VCF file and 
#	implementing further treatments on the SFS files.
# This file provides instructions as to how to use these functions.

# Load the functions.
source("vcf2sfs.r")

# Create a dadi format SFS from a VCF file.
vcf2dadi(f.vcf, f.popmap, f.output, pops,
		ploidy=2, n.digit=4,
		filter.indi=NA, filter.snp=NA)
# f.vcf: name of the VCF file.
# f.popmap: name of the popmap file (the same format as for Stacks).
# f.output: name of the output (SFS file).
# pops: the populations to be include in the SFS; population IDs are the same as in popmap.
# ploidy: ploidy of the species; normally default is fine.
# n.digit: number of digits after decimal for allele frequencies; normally default is fine.
# filter.indi: filter individuals according to missing values; 
#	when the number of missing values is larger than this, 
#	the individual will be filtered out; 
#	if NA, no filtering will be carried out.
# filter.snp: filter SNPs according to missing values; 
#	when the number of missing values is larger than this, 
#	the SNP will be filtered out; 
#	if NA, no filtering will be carried out.

# When you don't know how to set filter.indi and filter.snp in the previous function,
#	use this function to find out; it creats three plots.
checkMissing(f.vcf,f.popmap,pops)

# Change the SFS (one population) from dadi format to fastsimcoal format, and fold; 
#	meanwhile, a plot will be generated for review; 
#	the entry for 0 is not reliable.
dadi2fsc.1pop(f.fs,f.output,pop,fold=T)

# Change the SFS (two populations) from dadi format to fastsimcoal format, and fold; 
#	meanwhile, a plot will be generated for review; 
#	the entry for (0,0) is not reliable.
dadi2fsc.2pop(f.fs,f.output,pops,fold=T)

# Create bootstrapped SFS files.
fsc.sample(f.fsc,n.rep,nameroot)
# f.fsc: input file (fastSimCoal format).
# n.rep: number of repeats for the bootstrap.
# nameroot: nameroot for the output files.

