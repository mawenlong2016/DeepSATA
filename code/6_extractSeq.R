# The script used for extracting sequence from vcf or bed files
# Example: Rscript extractSeq.R example.bed pig
# Format of 4_trainBin.bed just as follows:
# chr1    3601    4600
# Example: Rscript extractSeq.R example.vcf pi
# Format of infile.vcf just as follows:
# chr1  109817590 [known_CEBP_binding_increase] G T

# The following should be optimized by automatically extract genome information
##########################################################
# install BSgenome to download Sscrofa reference genome
# BiocManager::install("BSgenome")
# library(BSgenome)
# (ag <- available.genomes())
# unique(gsub("BSgenome\\.([^\\.]+).*", "\\1", ag))
# BiocManager::install("BSgenome.Sscrofa.UCSC.susScr11")
##########################################################
require('rtracklayer')
#process input arguments
args <- commandArgs(trailingOnly = TRUE)

if(args[2] == 'pig'){
  require('BSgenome.Sscrofa.UCSC.susScr11')  
}else if(args[2] == 'mouse'){
  require('BSgenome.Mmusculus.UCSC.mm10')
  }else if(args[2] == 'chicken'){
    require('BSgenome.Ggallus.UCSC.galGal6')
    }else if(args[2] == 'cattle'){
      require('BSgenome.Btaurus.UCSC.bosTau9')
      }else{
        print("Currently, we only support pig, mouse, chicken, and cattle")
        print("Other species please wait for our updates!")
        break
      }

process.simple.offsets.python<-function(DF,prefix="./",window=1000,filterUnique=T){
  colnames(DF)[1:6]<-c("A","B","C","D","E","F")
  DF<-DF[!is.na(DF$D),]
  if(filterUnique)
    DF=DF[!duplicated(paste(as.character(DF$C),DF$D,DF$A,DF$B)),]
  DF$D=as.numeric(as.character(as.matrix(DF$D)))  
  halfw=window/2  

  if(args[2] == 'pig'){
    #filter by chromosome length
    chrends<-seqlengths(BSgenome.Sscrofa.UCSC.susScr11)[match(DF$C,names(seqlengths(BSgenome.Sscrofa.UCSC.susScr11)))]
    DF=DF[ (as.numeric(as.character(DF$D))-(halfw-1))>0 &  (as.numeric(as.character(DF$D))+halfw)<chrends & (!is.na(chrends)),]
    DF.rd1000<-GRanges(ranges=IRanges(start=as.numeric(as.character(DF$D))-(halfw-1),as.numeric(as.character(DF$D))+halfw),seqnames=as.character(DF$C),ori=DF$A,mut=DF$B,name=DF$F)
    #retrieve reference genome sequences
    DF.wt1000<-as.character(as.matrix(as.data.frame(getSeq(x=Sscrofa,names=as(DF.rd1000,"GRanges")))))    
  }
  if(args[2] == 'mouse'){
    #filter by chromosome length
    chrends<-seqlengths(BSgenome.Mmusculus.UCSC.mm10)[match(DF$C,names(seqlengths(BSgenome.Mmusculus.UCSC.mm10)))]
    DF=DF[ (as.numeric(as.character(DF$D))-(halfw-1))>0 &  (as.numeric(as.character(DF$D))+halfw)<chrends & (!is.na(chrends)),]
    DF.rd1000<-GRanges(ranges=IRanges(start=as.numeric(as.character(DF$D))-(halfw-1),as.numeric(as.character(DF$D))+halfw),seqnames=as.character(DF$C),ori=DF$A,mut=DF$B,name=DF$F)
    #retrieve reference genome sequences
    DF.wt1000<-as.character(as.matrix(as.data.frame(getSeq(x=Mmusculus,names=as(DF.rd1000,"GRanges")))))    
  }
  if(args[2] == 'chicken'){
    #filter by chromosome length
    chrends<-seqlengths(BSgenome.Ggallus.UCSC.galGal6)[match(DF$C,names(seqlengths(BSgenome.Ggallus.UCSC.galGal6)))]
    DF=DF[ (as.numeric(as.character(DF$D))-(halfw-1))>0 &  (as.numeric(as.character(DF$D))+halfw)<chrends & (!is.na(chrends)),]
    DF.rd1000<-GRanges(ranges=IRanges(start=as.numeric(as.character(DF$D))-(halfw-1),as.numeric(as.character(DF$D))+halfw),seqnames=as.character(DF$C),ori=DF$A,mut=DF$B,name=DF$F)
    #retrieve reference genome sequences
    DF.wt1000<-as.character(as.matrix(as.data.frame(getSeq(x=Ggallus,names=as(DF.rd1000,"GRanges")))))    
  }
  if(args[2] == 'cattle'){
    #filter by chromosome length
    chrends<-seqlengths(BSgenome.Btaurus.UCSC.bosTau9)[match(DF$C,names(seqlengths(BSgenome.Btaurus.UCSC.bosTau9)))]
    DF=DF[ (as.numeric(as.character(DF$D))-(halfw-1))>0 &  (as.numeric(as.character(DF$D))+halfw)<chrends & (!is.na(chrends)),]
    DF.rd1000<-GRanges(ranges=IRanges(start=as.numeric(as.character(DF$D))-(halfw-1),as.numeric(as.character(DF$D))+halfw),seqnames=as.character(DF$C),ori=DF$A,mut=DF$B,name=DF$F)
    #retrieve reference genome sequences
    DF.wt1000<-as.character(as.matrix(as.data.frame(getSeq(x=Btaurus,names=as(DF.rd1000,"GRanges")))))    
  }    
  #write to fasta format
  temp<-DNAStringSet(DF.wt1000)
  names(temp)<-paste(DF.rd1000$ori,DF.rd1000$mut,seqnames(DF.rd1000),start(DF.rd1000),end(DF.rd1000),DF.rd1000$name,sep = "_")
  writeXStringSet(temp,filepath=paste(prefix,".wt",window,".fasta",sep=""))
  print(DF.rd1000)
}
process.simple.offsets.python.bed<-function(data,prefix="./",window=1000){
  halfw=window/2
  
  pos = round((as.numeric(as.character(data[,2]))+as.numeric(as.character(data[,3])))/2)
  chrs = data[,1]
  chrs = paste("chr",gsub("chr","",chrs),sep="")
  DF.rd1000<-GRanges(ranges=IRanges(start=pos-(halfw-1),end=pos+halfw),seqnames=chrs,ori="",mut="")

  if(args[2] == 'pig'){
    #retrieve reference genome sequences
    DF.wt1000<-as.character(as.matrix(as.data.frame(getSeq(x=Sscrofa,names=as(DF.rd1000,"GRanges")))))
  }
  if(args[2] == 'mouse'){
    #retrieve reference genome sequences
    DF.wt1000<-as.character(as.matrix(as.data.frame(getSeq(x=Mmusculus,names=as(DF.rd1000,"GRanges")))))
  }
  if(args[2] == 'chicken'){
    #retrieve reference genome sequences
    DF.wt1000<-as.character(as.matrix(as.data.frame(getSeq(x=Ggallus,names=as(DF.rd1000,"GRanges")))))
  }
  if(args[2] == 'cattle'){
    #retrieve reference genome sequences
    DF.wt1000<-as.character(as.matrix(as.data.frame(getSeq(x=Btaurus,names=as(DF.rd1000,"GRanges")))))
  }
  
  #write to fasta format
  temp<-DNAStringSet(DF.wt1000)
  names(temp)<-paste(DF.rd1000$ori,DF.rd1000$mut,seqnames(DF.rd1000),start(DF.rd1000),end(DF.rd1000),sep = "_")
  writeXStringSet(temp,filepath=paste(prefix,".wt",window,".fasta",sep=""))
  print(DF.rd1000)
}
#read bed format
if(grepl("vcf$",args[1])){
  data=read.csv(args[1],sep='\t',header=F, comment.char = "#",colClasses = c("character"))
  data[,2]=as.numeric(data[,2])
  data=data[,c(4,5,1,2,2,3)]
  colnames(data)[1:6]<-c("A","B","C","D","E","F")
  data=data[order(data$C,data$D),]
  #for vcf format, we retrieve 1100bp sequences (instead of 1000bp) to allow for small deletions (<100bp)
  process.simple.offsets.python(data,window = 1000,prefix=args[1], filterUnique = F)
}else if(grepl("bed$",args[1])){
  data=read.csv(args[1],sep='\t',header=F, comment.char = "#",colClasses = c("character"))
  data[,2]=as.numeric(data[,2])
  data[,3]=as.numeric(data[,3])
  colnames(data)[1:3]<-c("A","B","C")
  data=data[order(data$A,data$B),]
  process.simple.offsets.python.bed(data,window = 1000,prefix=args[1])
}

