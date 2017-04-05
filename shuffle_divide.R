# arguments
args <- commandArgs(trailingOnly=T)

tf <- as.character(args[1])
outFolder <- as.character(args[2])
k <- as.numeric(args[3])
featureList <- args[4]

#listFile <- '../encoded/mad-max-myc-2all/list.txt'
#outFolder <- '../input/mad-max-myc-2all'
#k <- 10

# re-assemble outFolder (to avoid problem with tailing '/')
subs <- unlist(strsplit(outFolder, '/'))
outFolder <- paste0(paste0(subs, collapse='/'), '/')

# create output folder, print out warning if folder already existed
dir.create(outFolder)

# recover inFolder
subs <- unlist(strsplit(outFolder, '/'))
inFolder <- paste0(paste0(head(subs, -1), collapse='/'), '/', 'encoded/')

# what combinations of features are to be made ...
cmb <- read.table(featureList, header=F, colClasses=c('character', 'character'))
combinations <- cmb[, 2]
#combinations <- c('10000000000', '10011110000', '11000000000', '11011110000', '10100000000', '10111110000', '11100000000')
# ...

re_shuffle <- TRUE
cb <- combinations[1]
# read data
encoded_data <- as.matrix(read.table(paste0(inFolder, tf, '.', cb), header=F))
# randomly divide data into k subsets
m <- nrow(encoded_data)
if(re_shuffle){ # this makes sure different feature combinations use the same shuffling
shuffle <- sample.int(m, m)
re_shuffle <- FALSE
}
subsize <- floor(m/k)
for(i in 1 : k){
if(i != k){
  testIndice <- c(((i-1) * subsize + 1) : (i*subsize)) 
}else{
  testIndice <- c(((i-1) * subsize + 1) : m)
}
trainSet <- shuffle[-testIndice]
testSet <- shuffle[testIndice]
trainFilename <- paste0(outFolder, tf, '.fold', formatC(i, width=2, flag="0"), '.train')
testFilename <- paste0(outFolder, tf, '.fold', formatC(i, width=2, flag="0"), '.test')
write.table(matrix(trainSet, ncol=1), trainFilename, quote=F, sep="\t", row.names=F, col.names=F)
write.table(matrix(testSet, ncol=1), testFilename, quote=F, sep="\t", row.names=F, col.names=F)  
}

