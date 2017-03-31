# some functions ||
encode_1mer <- function(seq, map){
  features <- c()
  n <- nchar(seq)
  seq <- toupper(seq)
  for (i in 1 : n)
    features <- c(features, map[substr(seq, i, i),])

  return(features)
}


complDNA <- function(seq){
    tbl <- list()
    tbl[["a"]] <- "t"
    tbl[["c"]] <- "g"
    tbl[["g"]] <- "c"
    tbl[["t"]] <- "a"
    tbl[["A"]] <- "T"
    tbl[["C"]] <- "G"
    tbl[["G"]] <- "C"
    tbl[["T"]] <- "A"
    complseq <- ""
    for( i in nchar(seq) : 1 ){
        complseq <- paste(complseq, tbl[[substr(seq, i, i)]], sep="")
    }
    return(complseq)
}
# ||

# read arguments
args <- commandArgs(trailingOnly = TRUE)
seqFile <- args[1]
outFolder <- args[2]
palinCore <- F
coreStart <- 3  #if keepcore=T, of course, we need features for all seqs
coreLen <- as.numeric(args[3])  #
keepCore <- T
featureList <- args[4]

# some settings about the input file
cc <- c('character', 'numeric')
seqCol <- 1
sigCol <- 2
logarithm <- F
shuffleShape <- F
normalize <- T
threholdLow <- 8
excludeLow <- F
minSampleSize <- 1

# some glabal settings that don't need customization
global_mu_mgw <- 5.072129
global_mu_roll <- -0.6979492
global_mu_prot <- -6.791621
global_mu_helt <- 34.32604

global_sd_mgw <- 0.5213046
global_sd_roll <- 3.296222
global_sd_prot <- 3.237483
global_sd_helt <- 1.541358

# read data in ...
data <- read.table(seqFile, header=F, colClasses=cc)
## get sequences and signals ++
seqs <- toupper(data[, seqCol])
signals <- data[, sigCol]
## ++
# ...

# test if the sample size is big enough to continue |||
if(excludeLow){
  usedDataIdx <- which(signals > threholdLow) # only use data that has signal > threholdLow
  if(length(usedDataIdx) < minSampleSize){
    stop('Sample size too small')
  }
  ## update seqs and signals to store only data above the threshold ---
  seqs <- seqs[usedDataIdx]
  signals <- signals[usedDataIdx]
  ## ---
}else{
  if(length(signals) < minSampleSize){
    stop('Sample size too small')
  }
}
# |||

# stop when there is signal < 1, if logarithm is set to be true ...
if(logarithm){
  if(sum(signals < 1) > 0){
    stop("You chose to logarize, but there were signal < 1 in the data.")
  }
}
# ...

# get identifier --
identifier <- ""
subs <- unlist(strsplit(seqFile, '/'))
identifierWithExt <- subs[length(subs)]
identifierSubs <- unlist(strsplit(identifierWithExt, '[.]'))
n <- length(identifierSubs)
if(n > 1){
  identifier <- paste0(identifierSubs[1:(n-1)], collapse='.')
}else{
  identifier <- identifierWithExt # there was no extension anyway
}
# --

# re-assemble outFolder (to avoid problem with tailing '/') --
subs <- unlist(strsplit(outFolder, '/'))
outFolder <- paste0(subs, collapse='/')
# --

# create filename for the temp fasta file +++
tmpFaFile <- paste0(outFolder, '/', identifier, '.fa')
# +++

# encode features
## construct mapping ...
map_1mer <- rbind(diag(4), diag(4))
row.names(map_1mer) <- c("A", "C", "G", "T", "a", "c", "g", "t")

## encode 1mer/2mer/3mer features...
if(!palinCore){
  ### pre-alloc ...
  n <- nchar(seqs[1])
  m <- length(seqs)
  feature_1mer <- matrix(data=0, nrow=m, ncol=4*n)
  ### ...
  for (i in 1 : m){
    feature_1mer[i,] <- encode_1mer(seqs[i], map_1mer)
    cat(paste0('>', i), file=tmpFaFile, append=T, fill=T)
    cat(seqs[i], file=tmpFaFile, append=T, fill=T)
    #### ---
  }
  if(!keepCore){
    feature_1mer <- feature_1mer[, -(((coreStart-1)*4+1):((coreStart+coreLen-1)*4))]
}  
}else{
  ### pre-alloc ...
  n <- nchar(seqs[1])
  m <- length(seqs)
  feature_1mer <- matrix(data=0, nrow=m, ncol=4*n)
  ### ...
  for (i in 1 : m){
    feature_1mer[i,] <- encode_1mer(seqs[i], map_1mer) + encode_1mer(complDNA(seqs[i]), map_1mer)
    #### --- prepare fasta file
    cat(paste0('>', i), file=tmpFaFile, append=T, fill=T)
    cat(seqs[i], file=tmpFaFile, append=T, fill=T)
    #### ---
  }
  if(keepCore){
    feature_1mer <- feature_1mer[, 1:(n/2*4)]
	}else{
    feature_1mer <- feature_1mer[, 1:((n/2-coreLen/2)*4)]
}
}
## ...

# encode shape features ++
## read in shape prediciton ...
mgw <- as.matrix(read.table(paste0(tmpFaFile, ".MGW"), header=F, comment.char=">", na.strings="NA"))
prot <- as.matrix(read.table(paste0(tmpFaFile, ".ProT"), header=F, comment.char=">", na.strings="NA"))
roll <- as.matrix(read.table(paste0(tmpFaFile, ".Roll"), header=F, comment.char=">", na.strings="NA"))
helt <- as.matrix(read.table(paste0(tmpFaFile, ".HelT"), header=F, comment.char=">", na.strings="NA"))
## ...

if(!palinCore){
  if(!keepCore){ # if remove core
    if(coreLen >= 5){
      mgw <- mgw[, -c((coreStart+2):(coreStart+coreLen-1-2))]
      prot <- prot[, -c((coreStart+2):(coreStart+coreLen-1-2))]
      if(coreLen >= 6){
        roll <- roll[, -c((coreStart+2):(coreStart+coreLen-2-2))]
        helt <- helt[, -c((coreStart+2):(coreStart+coreLen-2-2))]
      }      
    }
  }  
  ## remove NAs ...
  endCol <- ncol(mgw)
  feature_mgw <- mgw[, -c(1,2,endCol-1,endCol)]
  feature_prot <- prot[, -c(1,2,endCol-1,endCol)]
  
  endCol <- ncol(roll)
  feature_roll <- roll[, -c(1,endCol)]
  feature_helt <- helt[, -c(1,endCol)]
  ## ...
}else{
  if(!keepCore){ # if remove core
    if(coreLen >= 5){
      mgw <- mgw[, -c((coreStart+2):(coreStart+coreLen-1-2))]
      prot <- prot[, -c((coreStart+2):(coreStart+coreLen-1-2))]
      if(coreLen >= 6){
        roll <- roll[, -c((coreStart+2):(coreStart+coreLen-2-2))]
        helt <- helt[, -c((coreStart+2):(coreStart+coreLen-2-2))]
      }      
    }
  }
  ## symmetrize and remove NAs and take first half...
  endCol <- ncol(mgw)
  mgw <- (mgw + mgw[, endCol:1]) /2
  mgw <- mgw[, -c(1,2,endCol-1,endCol)]
  feature_mgw <- mgw[, 1:ceiling(endCol/2-2)]
  prot <- (prot + prot[, endCol:1]) /2
  prot <- prot[, -c(1,2,endCol-1,endCol)]
  feature_prot <- prot[, 1:ceiling(endCol/2-2)]

  endCol <- ncol(roll)
  roll <- (roll + roll[, endCol:1]) /2
  roll <- roll[, -c(1,endCol)]
  feature_roll <- roll[, 1:ceiling(endCol/2-1)]
  helt <- (helt + helt[, endCol:1]) /2
  helt <- helt[, -c(1,endCol)]
  feature_helt <- helt[, 1:ceiling(endCol/2-1)]  
  ## ...
}
# ++

# clean up ...
junkFiles <- paste0(outFolder, '/', identifier, '.fa*')
cmd <- paste0('/bin/rm -f ', junkFiles)
system(cmd)
# ...

# shuffle or not ...
if(shuffleShape){
  n <- ncol(feature_mgw)
  for(i in 1 : n){
    feature_mgw[, i] <- sample(feature_mgw[, i])
    feature_prot[, i] <- sample(feature_prot[, i])    
  }
  n <- ncol(feature_roll)
  for(i in 1 : n){
    feature_roll[, i] <- sample(feature_roll[, i])
    feature_helt[, i] <- sample(feature_helt[, i])    
  }  
}
# ...

# normalize or not ...
if(normalize){

  #avr <- mean(as.vector(feature_mgw, mode='numeric'))
  #std <- sd(as.vector(feature_mgw, mode='numeric'))  
  #feature_mgw <- (feature_mgw - global_mu_mgw) / global_sd_mgw
  feature_mgw <- (feature_mgw - 2.85) / (6.2-2.85)
  #avr <- mean(as.vector(feature_roll, mode='numeric'))
  #std <- sd(as.vector(feature_roll, mode='numeric'))  
  #feature_roll <- (feature_roll - global_mu_roll) / global_sd_roll
  feature_roll <- (feature_roll - (-8.57)) / (8.64-(-8.57))
  #avr <- mean(as.vector(feature_prot, mode='numeric'))
  #std <- sd(as.vector(feature_prot, mode='numeric'))  
  #feature_prot <- (feature_prot - global_mu_prot) / global_sd_prot
  feature_prot <- (feature_prot - (-16.51)) / (-0.03-(-16.51))
  #avr <- mean(as.vector(feature_helt, mode='numeric'))
  #std <- sd(as.vector(feature_helt, mode='numeric'))  
 # feature_helt <- (feature_helt - global_mu_helt) / global_sd_helt
  feature_helt <- (feature_helt - 30.94)/(38.05 - 30.94)
  
}
# ...


# write combinations of features to files ++
library(bitops)

## mapping of features to binary barcodes ...
features <- list(feature_1mer, feature_2mer, feature_3mer, feature_mgw, feature_roll, feature_prot, feature_helt, feature_2nd_mgw, feature_2nd_roll, feature_2nd_prot, feature_2nd_helt)
pool <- c('10000000000', '01000000000', '00100000000', '00010000000', '00001000000', '00000100000', '00000010000', '00000001000', '00000000100', '00000000010', '00000000001')
## ...

## what combinations of features are to be made ...
cmb <- read.table(featureList, header=F, colClasses=c('character', 'character'))
combinations <- cmb[, 2]
#combinations <- c('10000000000', '10011110000', '11000000000', '11011110000', '10100000000', '10111110000', '11100000000')
## ...

if(logarithm){
  y <- log2(signals)
}else{
  y <- signals
}

for(i in combinations){
  pool_bin <- strtoi(pool, base=2)
  cb <- strtoi(i, base=2)
  features_selected <- do.call(cbind, features[which(bitAnd(pool_bin, cb) == pool_bin)])
  features_selected <- cbind(y,rep(c(1), nrow(features_selected)), features_selected) # add y column and constant column
  outFile <- paste0(outFolder, '/', identifier, '.', i)
  write.table(features_selected, outFile, sep=" ", quote=F, row.names=F, col.names=F)  
}  

listFile <- paste0(outFolder, '/list.txt')
cat(identifier, file=listFile, append=T, fill=T)
# ++
