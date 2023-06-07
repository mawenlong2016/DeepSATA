
# Rscript code/10_evaluation.R
predict <- read.csv(paste0("./4_testBin.bed.out"), header = T, stringsAsFactors=F)
test_label = read.table(paste0("./5_testLabel.bed"), header = F, stringsAsFactors=F)
colnames(test_label) <- c("bin_names","label_names","label")
require(reshape2)
b1 <- dcast(test_label[1:(nrow(test_label)*1),], bin_names~label_names, value.var="label")
rownames(b1) <- b1[,1]
b1 <- b1[,-1]
b1 <- as.matrix(b1)
b1[which(is.na(b1) == T)] = 0
rownames(predict) <- paste0(predict[,2],"_",predict[,3],"_",predict[,4])
flag1 <- rownames(b1)
if(nrow(b1) < nrow(predict)){
	b1 <- rbind(b1, matrix(0, nrow = nrow(predict) - nrow(b1), ncol = ncol(b1)))
	rownames(b1)[(length(flag1)+1) : nrow(b1)] <- setdiff(rownames(predict), flag1)
}
b1 <- b1[rownames(predict),]

predict <- predict[,c(colnames(predict)[1:4],colnames(b1))]
library(ROCR)
sum2 = apply(b1, 2, sum)
num2= which(apply(b1, 2, sd) > 0)
AUC1=sapply(num2, function(ii){
	pred <- prediction(as.numeric(predict[,ii+4]), b1[,ii])	
	auc <- performance(pred,'auc')
	return(unlist(slot(auc,"y.values")))
})

fea_auroc <- cbind(min=min(AUC1),
					max=max(AUC1),
					mean=mean(AUC1),
					median=median(AUC1),
					low0.5=length(which(AUC1 <= 0.5)),
					high0.9=length(which(AUC1 >= 0.9)),
					effect=length(AUC1),
					number=nrow(b1),
					num0=length(which(sum2 == 0)),
					num1=length(which(sum2 > 0)),
					t(AUC1)
					)


sum1 = apply(b1, 1, sum)
num1= which(apply(b1, 1, sd) > 0)
AUC3=sapply(num1, function(ii){
	pred <- prediction(as.numeric(predict[ii,5:ncol(predict)]), b1[ii,])	
	auc <- performance(pred,'auc')
	return(unlist(slot(auc,"y.values")))
})

bin_auroc <- cbind(min=min(AUC3),
					max=max(AUC3),
					mean=mean(AUC3),
					median=median(AUC3),
					low0.5=length(which(AUC3 <= 0.5)),
					high0.9=length(which(AUC3 >= 0.9)),
					effect=length(AUC3),
					number=nrow(b1),
					num0=length(which(sum1 == 0)),
					num1=length(which(sum1 > 0))
					)


library(PRROC)
AUC2=sapply(num2, function(ii){
	pr <- pr.curve(scores.class0 = as.numeric(predict[which(b1[,ii] == 1),ii+4])
		           , scores.class1 = as.numeric(predict[which(b1[,ii] == 0),ii+4]))
	return(pr$auc.integral)
})

fea_auprc <- cbind(min=min(AUC2),
					max=max(AUC2),
					mean=mean(AUC2),
					median=median(AUC2),
					low0.5=length(which(AUC2 <= 0.5)),
					high0.9=length(which(AUC2 >= 0.9)),
					effect=length(AUC2),
					number=nrow(b1),
					num0=length(which(sum2 == 0)),
					num1=length(which(sum2 > 0)),
					t(AUC2)
					)

AUC4=sapply(num1, function(ii){
	pr <- pr.curve(scores.class0 = as.numeric(predict[ii,which(b1[ii,] == 1)+4])
		           , scores.class1 = as.numeric(predict[ii,which(b1[ii,] == 0)+4]))
	return(pr$auc.integral)
})
bin_auprc <- cbind(min=min(AUC4),
					max=max(AUC4),
					mean=mean(AUC4),
					median=median(AUC4),
					low0.5=length(which(AUC4 <= 0.5)),
					high0.9=length(which(AUC4 >= 0.9)),
					effect=length(AUC4),
					number=nrow(b1),
					num0=length(which(sum1 == 0)),
					num1=length(which(sum1 > 0))
					)
bin_result <- rbind(auroc=bin_auroc,auprc=bin_auprc)
fea_result <- rbind(auroc=fea_auroc,auprc=fea_auprc)
write.table(bin_result, file = paste0("./10_bin_result.txt"),
			append=F, sep="\t",quote=F,row.names=T,col.names=T)
write.table(fea_result, file = paste0("./10_fea_result.txt"),
			append=F, sep="\t",quote=F,row.names=T,col.names=T)
save(list=ls(), file = paste0("./10_evaluation_final.RData"))




