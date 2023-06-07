# The script used for splitting chromosome into bins with size of 200bp
# Example: Rscript 1_setBin.R ChrPos.txt
# Format of ChrPos.txt just as follows:
# chr1 274330532
# chr2 151935994
# chr3 132848913

args <- commandArgs(trailingOnly = TRUE)
# read the position of each chromosome
ChrPos <- read.table(args[1], stringsAsFactors=F)
ChrBinBed <- NULL
for(i in 1:nrow(ChrPos)){
	num = as.numeric(ChrPos[i,2])
	a1 <- seq(401, num-400, 200) # interval of 200 bp without begin and end
	ChrBinBed <- rbind(ChrBinBed, cbind(ChrPos[i,1], a1[-length(a1)], a1[-1]-1))	
}
write.table(ChrBinBed, file = paste0(args[1],".bin.bed")
	,sep="\t",  quote=F, row.names=F,col.names=F)



