# vcf2sfs
R functions for calculating site frequency spectrum (SFS) or joint SFS of more than one population from a VCF file.

# Input and output formats
The main transformation function is vcf2dadi(). It requires two input files. One is the VCF file, and the other is a so-called popmap file. The popmap file is a tab-delimited two-column text file. Its first column specifies the IDs of the individuals included in the VCF file (the individuals should be arranged in the same order as in the VCF file; the IDs do not need to be the same as those in the header line of the VCF file). The second column specifies the populaton ID of each individual (coded in integers). The output of vcf2dadi() is in dadi format. But there are extra functions to transform dadi format into fastSimCoal format.

# Instructions
For more detailed instructions, please see "readme.r"

# Example
I provide an example VCF file and an example popmap file. The dataset includes 197 individuals belonging to 11 populations.

To create the SFS of Population 1, you simply need to run:
  source("vcf2sfs.r")
  vcf2dadi("example.vcf","example_popmap.txt","pop1.fs",1)

To create the joint SFS between Populations 1 and 2:
  vcf2dadi("example.vcf","example_popmap.txt","pop1_pop2.fs",c(1,2))

To create the joint SFSs for all possible pairs or populations:
  popPairs<-combn(1:11,2)
  for(i in 1:ncol(popPairs))
  vcf2dadi("example.vcf","example_popmap.txt",sprintf("pop%d_pop%d.fs",popPairs[1,i],popPairs[2,i]),popPairs[,i])

# Citation
Liu S, Ferchaud AL, Groenkjaer P, Nygaard R, Hansen MM (2018) Genomic parallelism and lack thereof in contrasting systems of three-spine sticklebacks. Molecular Ecology, Vol. 27 (in press).
