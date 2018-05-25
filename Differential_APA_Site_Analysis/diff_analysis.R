#install.packages("matrixStats", repos="http://cran.rstudio.com/")
library(matrixStats)
#install.packages("locfit", repos="http://cran.rstudio.com/")
library(locfit)
#install.packages("locfit", repos="http://cran.rstudio.com/")
library(stats)
load("scvBiasCorrectionFits.rda")


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
  
  # second parameter 2 means the function is applied to the column of the matrix
  sizeFactor1 <- apply( condition1, 2, function(cnts)
    exp( median( ( log(cnts) - denominator )[ is.finite(denominator)] ) ) )

 #print(condition1)
  
  sizeFactor2 <- apply( condition2, 2, function(cnts)
    exp( median( ( log(cnts) - denominator )[ is.finite(denominator)] ) ) )
  
  return (list(sfCond1 = sizeFactor1, sfCond2 = sizeFactor2))
}

gettingMean <- function(givenData, sizeFactor) {
  return (rowMeans( t( t(givenData) / sizeFactor )))
}

gettingVariance <- function(givenData, sizeFactor) {
  return (rowVars( t( t(givenData) / sizeFactor )))
}

adjustScvForBias <- function( scv, nsamples ) {
   stopifnot( nsamples > 1 )
   if( nsamples - 1 > length( scvBiasCorrectionFits ) )
      scv
   else
      ifelse( scv > .02,
         pmax( safepredict( scvBiasCorrectionFits[[ nsamples-1 ]], scv ), 1e-8 * scv ),
         scv )   # For scv < .02, our fit is too coarse, but no correction seems necessary anyway
}

safepredict <- function( fit, x ) {
  # A wrapper around predict to avoid the issue that predict.locfit cannot
  # propagate NAs and NaNs properly.
  
  res <- rep.int( NA_real_, length(x) )
  res[ is.finite(x) ] <- predict( fit, x[is.finite(x)] )
  res
}

parametricDispersionFit <- function( means, disps )
{
   coefs <- c( .1, 1 )
   iter <- 0
   while(TRUE) {
      residuals <- disps / ( coefs[1] + coefs[2] / means )
      good <- which( (residuals > 1e-4) & (residuals < 15) )
      fit <- glm( disps[good] ~ I(1/means[good]),
         family=Gamma(link="identity"),na.action=na.exclude, start=coefs )
      oldcoefs <- coefs
      coefs <- coefficients(fit)
      if( !all( coefs > 0 ) )
         stop( "Parametric dispersion fit failed. Try a local fit and/or a pooled estimation. (See '?estimateDispersions')" )
      if( sum( log( coefs / oldcoefs )^2 ) < 1e-6 )
         break
      iter <- iter + 1
      if( iter > 10 ) {
         warning( "Dispersion fit did not converge." )
         break }
   }

   names( coefs ) <- c( "asymptDisp", "extraPois" )
   ans <- function( q )
      coefs[1] + coefs[2] / q
   attr( ans, "coefficients" ) <- coefs
   ans
}

estimateAndFitDispersion <- function(means, variances, sizeFactors) {

  xim <- mean(1/ sizeFactors)
  dispsAll <- ( variances - xim * means ) / means^2

  variances <- variances[means > 0]
  disps <- dispsAll[means > 0]
  means <- means[means > 0]
  
  fit <- do.call( "locfit", c(
                    list(
                    	variances ~ do.call( "lp", c( list( log(means) ), list()) )),
                        family = "gamma") )

  rm (means)
  rm (variances)
  ans <- function( q )
	pmax( ( safepredict( fit, log(q) ) - xim * q ) / q^2, 1e-8 )

  return (list(disps=dispsAll, dispFunc=ans))
}

mainFunction <- function() {
  fileName1 <- "sample1.txt"
  fileName2 <- "sample2.txt"
  
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

  		varCond1 <- gettingVariance(dataFile1[, c(5:ncol(dataFile1))], sizeFactorList$sfCond1)
  		varCond2 <- gettingVariance(dataFile2[, c(5:ncol(dataFile2))], sizeFactorList$sfCond2)
 
		varCond1 <- pmax(varCond1, 0.00000000000001)
		varCond2 <- pmax(varCond2, 0.00000000000001)

  		dispAndFunc1 <- estimateAndFitDispersion(meanCond1, varCond1, sizeFactorList$sfCond1)

		if (is.null(dispAndFunc1$dispFunc)) {
			print('######A')
			finalDisp1 <- dispAndFunc1$disps
  		} else {
			print('@@@@A')
			fittedDispEsts1 <- dispAndFunc1$dispFunc(meanCond1)
			finalDisp1 <- pmax(dispAndFunc1$disps, fittedDispEsts1)
		}

		dispAndFunc2 <- estimateAndFitDispersion(meanCond2, varCond2, sizeFactorList$sfCond2)
	
 		if (is.null(dispAndFunc2$dispFunc)) {
			print('#####B')
			finalDisp2 <- dispAndFunc2$disps
		} else {
			print('@@@@@@B')
			fittedDispEsts2 <- dispAndFunc2$dispFunc(meanCond2)
			finalDisp2 <- pmax(dispAndFunc2$disps, fittedDispEsts2)
		}

		kAs <- rowSums(dataFile1[, c(5:ncol(dataFile1))])
		kBs <- rowSums(dataFile2[, c(5:ncol(dataFile2))])	

		# for null hypothesis two conditions are combined to calculate the mu value
		mus <- rowMeans( cbind(t( t(dataFile1[, c(5:ncol(dataFile1))]) / sizeFactorList$sfCond1 ), t( t(dataFile2[, c(5:ncol(dataFile2))]) / sizeFactorList$sfCond2 ) ) )

		# the common mu value is used to calculate the variance of individual condition
   		fullVarsA <- pmax( mus * sum(sizeFactorList$sfCond1) + dispAndFunc1$dispFunc(mus) * mus^2 * sum(sizeFactorList$sfCond1^2), mus * sum(sizeFactorList$sfCond1) * (1+1e-8) )
		fullVarsB <- pmax( mus * sum(sizeFactorList$sfCond2) + dispAndFunc2$dispFunc(mus) * mus^2 * sum(sizeFactorList$sfCond2^2), mus * sum(sizeFactorList$sfCond2) * (1+1e-8) )

		# this gives the dispersion for the common mu to individual condition
		sumDispsA <- ( fullVarsA - mus * sum( sizeFactorList$sfCond1 ) ) / ( mus * sum( sizeFactorList$sfCond1 ) )^2
		sumDispsB <- ( fullVarsB - mus * sum( sizeFactorList$sfCond1 ) ) / ( mus * sum( sizeFactorList$sfCond2 ) )^2

		pval <- sapply( seq(along=kAs), function(i) {

      			if( kAs[i] == 0 & kBs[i] == 0 ) {
         			return( NA )
			}

      			# probability of all possible counts sums with the same total count:
      			ks <- 0 : ( kAs[i] + kBs[i] )
      				ps <- dnbinom(ks, mu = mus[i] * sum( sizeFactorList$sfCond1 ), size = 1/sumDispsA[i] ) *
            				dnbinom( kAs[i] + kBs[i] - ks, mu = mus[i] * sum( sizeFactorList$sfCond2 ), size = 1/sumDispsB[i] )

      			# probability of observed count sums:
      			pobs <- dnbinom( kAs[i], mu = mus[i] * sum( sizeFactorList$sfCond1 ), size = 1/sumDispsA[i] ) *
              			dnbinom( kBs[i], mu = mus[i] * sum( sizeFactorList$sfCond2 ), size = 1/sumDispsB[i] )

			if( kAs[i] * sum( sizeFactorList$sfCond2 ) < kBs[i] * sum( sizeFactorList$sfCond1 ) ) {
         			numer <- ps[ 1 : (kAs[i]+1) ]
			} else {
         			numer <- ps[ (kAs[i]+1) : length(ps) ]
			}

      			min( 1, 2 * sum(numer) / sum(ps) )
			
		})
	
		foldChange <- meanCond2 / meanCond1
		log2FoldChange <- log2(meanCond2 / meanCond1)
		padj <- p.adjust( pval, method="BH" )
		#padj <- format(round(p.adjust(pval, method="BH"), 3), nsmall=3)

		write.table(cbind(dataFile1[, c(1:4)], foldChange, log2FoldChange, pval, padj), file = "diff_result.txt", sep="\t", eol="\n", quote = FALSE, row.names=FALSE, col.names=FALSE)
	}
  } else {
	foldChange <- dataFile1[, c(5:ncol(dataFile1))] / dataFile2[, c(5:ncol(dataFile2))]
        log2FoldChange <- abs(log2(dataFile1[, c(5:ncol(dataFile1))]) - log2(dataFile2[, c(5:ncol(dataFile2))]))
	write.table(cbind(dataFile1[, c(1:4)], foldChange, log2FoldChange), file = "diff_result.txt", sep="\t", eol="\n", quote = FALSE, row.names=FALSE, col.names=FALSE)
  }
}

