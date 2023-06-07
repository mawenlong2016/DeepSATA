args <- commandArgs(trailingOnly = TRUE)
motif <- read.table('motif.name', header=F)[,1]
options(stringsAsFactors=F)
options(scipen=999)
suppressPackageStartupMessages(library(doParallel))
suppressPackageStartupMessages(library(foreach))
cl <-  makeCluster(as.numeric(args[1]))
registerDoParallel(cl)
#do parallel computation
tempList <- foreach(i = 1:length(motif)) %dopar%{
	system(paste0('grep ',motif[i],' ',args[2],' | wc -l | awk ','\'{print $1"\t""',motif[i],'"}\'',' > ',args[3],'_motif/',motif[i],'_motif.count'))
}
stopCluster(cl)
