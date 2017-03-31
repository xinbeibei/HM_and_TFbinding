rsqrd <- function(observe, model){
  ss_res <- sum((observe-model)^2)
  y_bar <- mean(observe)
  ss_tot <- sum((observe-y_bar)^2)
  return (1-ss_res/ss_tot)
}
################################

args <- commandArgs(trailingOnly=T)

encoded_file <- args[1]
train_index_file <- args[2]
outDirection <- args[3]
k <- as.numeric(args[4])   # k-fold cv

encoded_data <- as.matrix(read.table(encoded_file, header=F))
train_indice <- as.matrix(read.table(train_index_file, header=F, colClasses='numeric'))
train <- encoded_data[train_indice,]

m <- nrow(train)
n <- ncol(train)
s <- floor(m/k)

best_lmbd <- -1
best_acc <- -1
lambdas <- 2 ^ c(-15:15)
info <- matrix(data=NA, nrow=length(lambdas), ncol=4)
colnames(info) = c('lambda', 'acc_train_mean', 'acc_dev_mean', 'acc_dev_all')

lmbd_ind <- 1 
for (lmbd in lambdas){
y_dev_all_pred <- c()
acc_train <- rep(-1, k)
acc_dev <- rep(-1, k)    

for (i in 1 : k){ ## calculate training/developing accuracy for each fold; collect k-fold prediction ...
  if(i != k){
    sub <- c((s * (i-1) + 1) : (s * i))
  }else{
    sub <- c((s * (i-1) + 1) : m)
  }
  comp_sub <- setdiff(c(1:m), sub)
  x_dev_foldi <- train[sub, 2:n]
  y_dev_foldi <- train[sub, 1]      
  x_train_foldi <- train[comp_sub, 2:n]
  y_train_foldi <- train[comp_sub, 1]      
  omg <- solve(t(x_train_foldi) %*% x_train_foldi + lmbd*diag(n-1), t(x_train_foldi) %*% y_train_foldi)
  y_train_foldi_pred <- x_train_foldi %*% omg      
  y_dev_foldi_pred <- x_dev_foldi %*% omg
  acc_train[i] <- rsqrd(y_train_foldi, y_train_foldi_pred)
  acc_dev[i] <- rsqrd(y_dev_foldi, y_dev_foldi_pred)      
  y_dev_all_pred <- c(y_dev_all_pred, y_dev_foldi_pred)
}## ...

## calculate mean training/dev accuracy for the current lambda
info[lmbd_ind, 1] <- lmbd
info[lmbd_ind, 2] <- mean(acc_train)
info[lmbd_ind, 3] <- mean(acc_dev)    
## ...

## calculate overall performance for the current lambda ...
curr_acc <- rsqrd(train[, 1], y_dev_all_pred)
info[lmbd_ind, 4] <- curr_acc
lmbd_ind <- lmbd_ind + 1    
if(curr_acc > best_acc){
  best_acc <- curr_acc
  best_lmbd <- lmbd
}
## ...
}

omg <- solve(t(train[,2:n]) %*% train[,2:n] + best_lmbd*diag(n-1), t(train[,2:n]) %*% train[,1])

write.table(info, paste0(outDirection, '.log'), col.names=T, row.names=F, quote=F, sep="\t")
write.table(omg, paste0(outDirection, '.omg'), col.names=F, row.names=F, quote=F, sep=" ")
write.table(best_lmbd, paste0(outDirection, '.best_lambda'), col.names=F, row.names=F, quote=F, sep=" ")    
 
