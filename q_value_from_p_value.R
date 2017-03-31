args <- commandArgs(trailingOnly = TRUE)
library("MASS")
library("qvalue")

folder <- args[1]
tf <- args[2]

subs <- unlist(strsplit(folder, '/'))
folder <- paste0(paste0(subs, collapse='/'),'/')


bind_file <- paste0(folder, tf, '/BS_10_HM.txt')   #/home/cmb-04/rr/bxin/hi-c-human/$tf/BS_10_HM.txt
nobind_file <- paste0(folder, tf, '/non-BS_10_HM.txt')   #/home/cmb-04/rr/bxin/hi-c-human/$tf/non-BS_10_HM.txt

bind <- as.matrix(read.table(bind_file, header=FALSE))
unbind <- as.matrix(read.table(nobind_file, header=FALSE))
p_values <- c()
for (i in 1:10){
	p_values <- c(p_values,wilcox.test(as.numeric(bind[,i]), as.numeric(unbind[,i]), alternative = "less")$p.value)
	
}
qobj_greater <- p.adjust(p_values, method = "bonferroni")


p_values <- c()
for (i in 1:10){
  p_values <- c(p_values,wilcox.test(as.numeric(bind[,i]), as.numeric(unbind[,i]), alternative = "greater")$p.value)
}
qobj_less <- p.adjust(p_values, method = "bonferroni")

write.table(-log(t(qobj_greater)) + log(t(qobj_less)), paste0(folder, tf, '/qvalue_delta.txt') ,quote=FALSE, sep= "\t", append = TRUE, row.names=FALSE, col.names=FALSE)

