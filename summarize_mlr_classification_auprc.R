library("ROCR")
library("PRROC")
###############################
rsqrd <- function(observe, model){
  ss_res <- sum((observe-model)^2)
  y_bar <- mean(observe)
  ss_tot <- sum((observe-y_bar)^2)
  return (1-ss_res/ss_tot)
}

classification_accuracy <- function(label, prediction){
  prediction[which(prediction >= 0.5)] = 1
  prediction[which(prediction < 0.5)] = 0
  accuracy <- length(which(label==prediction)) / length(label)
  return(accuracy)
}
################################ 
args <- commandArgs(trailingOnly = TRUE)

#prefix <- '../tianyin-method/data/GATA/GATA3_TCATGC20NCG_GATA_10'
#outFile <- '../output/GATA/summary.txt'

prefix <- args[1]
outPath <- args[2]
featureList <- args[3]

outFile <- paste0(outPath, 'summary_auprc.txt') # fixed filename "summary.txt"

# settings
## what combinations of features are to be made ...
cmb <- read.table(featureList, header=F, colClasses=c('character', 'character'))
combinations <- cmb[, 2]
#combinations <- c('10000000000', '10011110000', '11000000000', '11011110000', '10100000000', '10111110000', '11100000000')
## ...
subs <- c('.fold01','.fold02','.fold03','.fold04','.fold05','.fold06','.fold07','.fold08','.fold09','.fold10')

r2 <- c()

for (i in 1 : nrow(cmb)){
    cb <- cmb[i,2]
    title <- cmb[i,1]
    test <- c()
    predict <- c()
    for (sb in subs){
    	t <- as.matrix(read.table(paste0(prefix, '.', cb, sb, '.test.y'), header=F, colClasses='numeric'))
    	p <- as.matrix(read.table(paste0(prefix, '.', cb, sb, '.test.p'), header=F, colClasses='numeric'))
	test <- c(test, t)
	predict <- c(predict, p)
    }
    
    pred <- prediction(predict, test)
    perf <- performance(pred, "prec", "rec")
    #aucpr <- performance(pred, "aucpr")
	scores <- data.frame(predict, test)
	aucpr <- pr.curve(scores.class0=scores[scores$test=="1",]$predict,
                scores.class1=scores[scores$test=="0",]$predict,
               curve=T)
    jpeg(paste0(outPath, basename(prefix), '_', cb, '.jpg'), width=6, height=6, units='in', res=600)
    plot(perf, main=paste0(title, " AUPRC: ", round(as.numeric(aucpr$auc.integral), 4)))
    dev.off()
	
	jpeg(paste0(outPath, basename(prefix), '_', cb, '_new.jpg'), width=6, height=6, units='in', res=600)
    plot(aucpr, main=paste0(title, " AUPRC: ", round(as.numeric(aucpr$auc.integral), 4)))
    dev.off()
	

    r2 <- c(r2, round(as.numeric(aucpr$auc.integral), 4))    
}

num_seq <- length(predict)
r2 <- c(num_seq, r2)
identifier <- basename(prefix)
r2 <- c(identifier, r2)

write(r2, outFile, append=T, ncolumns=length(r2))
