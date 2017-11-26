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
library("DNAshapeR")

# read arguments
args <- commandArgs(trailingOnly = TRUE)
seqFile <- args[1]
outFolder <- args[2]
palinCore <- F

coreLen <- as.numeric(args[3])  #
keepCore <- F
flank_len <- as.integer(args[4])
coreStart <- flank_len + 1          #if keepcore=T, of course, we need features for all seqs
# some settings about the input file
cc <- c('character', 'numeric')
seqCol <- 1
sigCol <- 2
logarithm <- F
shuffleShape <- F
normalize <- T
interact <- F
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

##############PREPARE FOR SEQUENCE ENCODING
# encode features
## construct mapping ...
map_1mer <- rbind(diag(4), diag(4))
row.names(map_1mer) <- c("A", "C", "G", "T", "a", "c", "g", "t")

if(!palinCore){
  ### pre-alloc ...
  n <- nchar(seqs[1])
  m <- length(seqs)
  feature_1mer <- matrix(data=0, nrow=m, ncol=4*n)
  ### ...
  for (i in 1 : m){
    feature_1mer[i,] <- encode_1mer(seqs[i], map_1mer)
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
  }
  if(keepCore){
    feature_1mer <- feature_1mer[, 1:(n/2*4)]      
  }else{
    feature_1mer <- feature_1mer[, 1:((n/2-coreLen/2)*4)]
  }
}

##############PREPARE THINGS FOR SHAPE ENCODING
# re-assemble outFolder (to avoid problem with tailing '/') --
subs <- unlist(strsplit(outFolder, '/'))
outFolder <- paste0(subs, collapse='/')
# --

# create filename for the temp fasta file +++
tmpFaFile <- paste0(outFolder, '/', identifier, '.fa')
# +++
## ...
m <- length(seqs)

#if fasta file exist, delete the file
if (file.exists(tmpFaFile)){
  cmd <- paste0('rm ', tmpFaFile)
  system(cmd)
}

for (i in 1 : m){    
    #### --- prepare fasta file
    cat(paste0('>', i), file=tmpFaFile, append=T, fill=T)
    cat(seqs[i], file=tmpFaFile, append=T, fill=T)
    #### ---
  }
# predict shape features --
pred = getShape(tmpFaFile)
# --

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
#junkFiles <- paste0(outFolder, '/', identifier, '.fa*')
#cmd <- paste0('/bin/rm -f ', junkFiles)
#system(cmd)
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


if(logarithm){
  y <- log2(signals)
}else{
  y <- signals
}


# write combinations of features to files ++
# before_avg_mgw <- cbind(feature_mgw[,1:(flank_len-4)], t(tail(t(feature_mgw),(flank_len-4))))   
# before_avg_roll <- cbind(feature_roll[,1:(flank_len-3)], t(tail(t(feature_roll),(flank_len-3))))   
# before_avg_prot <- cbind(feature_prot[,1:(flank_len-4)], t(tail(t(feature_prot),(flank_len-4))))   
# before_avg_helt <- cbind(feature_helt[,1:(flank_len-3)], t(tail(t(feature_helt),(flank_len-3))))   

# features_selected <- cbind(rowMeans(before_avg_mgw), rowMeans(before_avg_roll), rowMeans(before_avg_prot), rowMeans(before_avg_helt))\
#DNA flank seq features:
features_selected <- cbind(y,rep(c(1), nrow(feature_1mer)), feature_1mer)
outFile <- paste0(outFolder, '/', identifier, '.10000000000')
write.table(features_selected, outFile, sep=" ", quote=F, row.names=F, col.names=F)  

#DNA flank seq MGW:
#features_selected <- cbind(y,rep(c(1), nrow(features_selected)), feature_mgw)
outFile <- paste0(outFolder, '/', identifier, '.00010000000')
write.table(feature_mgw, outFile, sep=" ", quote=F, row.names=F, col.names=F)  

#DNA flank seq Roll:
#features_selected <- cbind(y,rep(c(1), nrow(features_selected)), feature_roll)
outFile <- paste0(outFolder, '/', identifier, '.00001000000')
write.table(feature_roll, outFile, sep=" ", quote=F, row.names=F, col.names=F)  

#DNA flank seq prot:
#features_selected <- cbind(y,rep(c(1), nrow(features_selected)), feature_roll)
outFile <- paste0(outFolder, '/', identifier, '.00000100000')
write.table(feature_prot, outFile, sep=" ", quote=F, row.names=F, col.names=F)  

#DNA flank seq helt:
#features_selected <- cbind(y,rep(c(1), nrow(features_selected)), feature_roll)
outFile <- paste0(outFolder, '/', identifier, '.00000010000')
write.table(feature_helt, outFile, sep=" ", quote=F, row.names=F, col.names=F)  

#all DNA shape combined.
features_selected <- cbind(y,rep(c(1), nrow(feature_mgw)), feature_mgw, feature_roll, feature_prot, feature_helt)
outFile <- paste0(outFolder, '/', identifier, '.00011110000')
write.table(features_selected, outFile, sep=" ", quote=F, row.names=F, col.names=F)  

listFile <- paste0(outFolder, '/list.txt')
cat(identifier, file=listFile, append=T, fill=T)
# ++
