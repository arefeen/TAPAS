library(matrixStats)
library(stats)

calculatingDenominator <- function(condition1, condition2) {
  sum1 = rowSums(log(condition1))
  sum2 = rowSums(log(condition2))

  return((sum1 + sum2) / (ncol(condition1) + ncol(condition2)))
}

readingData <- function(nameOfTheFile) {
  dataFile <- read.table(nameOfTheFile)

  return (dataFile)
}

sizeFactorFunction <- function(condition1, condition2, denominator) {

  sizeFactor1 <- apply( condition1, 2, function(cnts)
    exp( median( ( log(cnts) - denominator )[ is.finite(denominator)] ) ) )


  sizeFactor2 <- apply( condition2, 2, function(cnts)
    exp( median( ( log(cnts) - denominator )[ is.finite(denominator)] ) ) )

  return (list(sfCond1 = sizeFactor1, sfCond2 = sizeFactor2))
}

gettingMean <- function(givenData, sizeFactor) {
  return (rowMeans( t( t(givenData) / sizeFactor )))
}

log2FoldChangeCalculator <- function() {
  fileName1 <- "sample1_sl.txt"
  fileName2 <- "sample2_sl.txt"

  dataFile1 <- readingData(fileName1)
  dataFile2 <- readingData(fileName2)


  if ((ncol(dataFile1) >= 6) & (ncol(dataFile2) >= 6)) {
        denominatorOfBothSample <- calculatingDenominator(dataFile1[, c(5:ncol(dataFile1))], dataFile2[, c(5:ncol(dataFile2))])
        sizeFactorList <- sizeFactorFunction(dataFile1[, c(5:ncol(dataFile1))], dataFile2[, c(5:ncol(dataFile2))], denominatorOfBothSample)

        if( !any( is.na(sizeFactorList$sfCond1 ) ) || !any(is.na(sizeFactorList$sfCond2))) {
                meanCond1 <- gettingMean(dataFile1[, c(5:ncol(dataFile1))], sizeFactorList$sfCond1)
                meanCond2 <- gettingMean(dataFile2[, c(5:ncol(dataFile2))], sizeFactorList$sfCond2)

                meanCond1 <- pmax(meanCond1, 0.00000000000001)
                meanCond2 <- pmax(meanCond2, 0.00000000000001)
   	}
   }

   log2FoldChange <- log2(meanCond2 / meanCond1)
   write.table(cbind(dataFile1[, c(1:4)], log2FoldChange), file = "diff_result_for_sl.txt", sep="\t", eol="\n", quote = FALSE, row.names=FALSE, col.names=FALSE) 
}
