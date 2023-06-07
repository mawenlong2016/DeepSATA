# The script used for combining feature matrix and label into a single H5 file
# Example usage: Rscript 7_labelH5.R valid 7_valid.fasta.ref.h5 5_validLabel.bed
# Format of 7_valid.fasta.ref.h5 just as follows:
# 1	0	0	0
# 0	1	0	0
# 0	0	1	0
# 0	0	0	1
# 0	0	1	0
# 0	1	0	0
# 1	0	0	0
# 0	1	0	0
# 0	0	1	0
# 0	0	0	1

# Format of 5_validLabel.bed just as follows:
# chr7_30508801_30509800  D8_xiaochang_duluoke    1
# chr7_30508801_30509800  D9_xiaochang_duluoke    0
# chr7_30509001_30510000  D8_danao_duluoke        0

# process input arguments
args <- commandArgs(trailingOnly = TRUE)

label = read.table(args[3], header = F, stringsAsFactors=F)
colnames(label) <- c("bin_names","label_names","label")
require(reshape2)
b1 <- dcast(label[1:(nrow(label)*1),], bin_names~label_names, value.var="label")
rownames(b1) <- b1[,1]
b1 <- b1[,-1]
b1 <- as.matrix(b1)
b1[which(is.na(b1) == T)] = 0

system(command=paste0("grep '>' 4_", args[1], "Bin.bed.wt1000.fasta > 7_",args[1],"bin_name.txt"))
bin_name_tmp <- read.table(paste0("7_",args[1],"bin_name.txt"), header = F, stringsAsFactors=F)[,1]
bin_name <- do.call(rbind, strsplit(bin_name_tmp, "_"))[,3:5]
label_name <- read.table("label_name.txt", header = F, stringsAsFactors=F)[,1]

flag1 <- setdiff(paste0(bin_name[,1],"_",bin_name[,2],"_",bin_name[,3]), rownames(b1))
if(length(flag1) > 0){
  c1<- rbind(b1, matrix(0, length(flag1),ncol(b1)))
  rownames(c1)[(nrow(b1)+1):nrow(c1)] <- flag1
  b1 <- c1
}
flag2 <- setdiff(label_name, colnames(b1))
if(length(flag2) > 0){
	c1<- cbind(b1, matrix(0, nrow(b1),length(flag2)))
	colnames(c1)[(ncol(b1)+1):ncol(c1)] <- flag2
	b1 <- c1
}

b1 <- b1[paste0(bin_name[,1],"_",bin_name[,2],"_",bin_name[,3]),label_name]
b1 <- rbind(b1,b1)
require(rhdf5)
# if(file.exists(paste0("7_",args[1],".combine.ref.h5"))){
# 	system(command=paste0("rm 7_",args[1],".combine.ref.h5"))
# }

if(file.exists(paste0("7_",args[1],".label.ref.h5"))){
	system(command=paste0("rm 7_",args[1],".label.ref.h5"))
}

h5write(b1, file=paste0("7_",args[1],".label.ref.h5"), name="testdata")

#h5write(b1, file=paste0("7_",args[1],".combine.ref.h5"), name="testdata")
#c1=h5read(file=paste0("7_", args[1], '.fasta.ref.h5'), name="testxdata")
#h5write(c1, file=paste0("7_",args[1],".combine.ref.h5"), name="testxdata")

write.table(rownames(b1), file = paste0("7_",args[1],".bin.name")
	, quote=F, row.names=F,col.names=F)
write.table(colnames(b1), file = paste0('7_',args[1],".label.name")
	, quote=F, row.names=F,col.names=F)


