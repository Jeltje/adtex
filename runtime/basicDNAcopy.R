library(DNAcopy)

args <- commandArgs(TRUE)

# arguments are input, SD, and sampleID

# Alway use the same random seed for reproducible results
set.seed(0xcafe)	# cafe is a hex number

cn <- read.table(args[1],header=TRUE)
cn$logscore <- log2(cn$rationormalized_after_smoothing) 
sd <- as.double(args[2])
sample <- args[3]
CNA.object <-CNA(genomdat = cn$logscore, chrom = cn$chr, maploc = cn$exon_end, 
	data.type = 'logratio', sampleid = sample)

smoothed.CNA.object <- smooth.CNA(CNA.object)

segment.smoothed.CNA.object <- segment(smoothed.CNA.object, undo.splits="sdundo", undo.SD=sd, verbose=0)
p.segment.smoothed.CNA.object <- segments.p(segment.smoothed.CNA.object)

#outfile <- paste(args[1], "SD", sd, "dnacopy.out", sep=".")

write.table(p.segment.smoothed.CNA.object[,1:6], quote=F, row.names=F, sep="\t")

detach(package:DNAcopy)
