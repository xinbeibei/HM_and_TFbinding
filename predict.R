rsqrd <- function(observe, model){
  ss_res <- sum((observe-model)^2)
  y_bar <- mean(observe)
  ss_tot <- sum((observe-y_bar)^2)
  return (1-ss_res/ss_tot)
}
################################

args <- commandArgs(trailingOnly=T)

encoded_file <- args[1]
test_index_file <- args[2]
model_file <- args[3]
outDirection <- args[4]

encoded_data <- as.matrix(read.table(encoded_file, header=F))
test_indice <- as.matrix(read.table(test_index_file, header=F, colClasses='numeric'))
test <- encoded_data[test_indice,]


n <- ncol(test)  
y <- test[, 1] # real y's
model <- as.matrix(read.table(model_file, header=F)) # load omega 
p <- test[, 2:n] %*% model # calculate predicted y's

write.table(y, paste0(outDirection, '.y'), col.names=F, row.names=F, quote=F, sep=" ")
write.table(p, paste0(outDirection, '.p'), col.names=F, row.names=F, quote=F, sep=" ")  

