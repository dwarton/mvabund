# Prototype coefplot, based off coefplot from gllvm by Jenni Niku, Francis Hui and Sara Taskinen
# Plots the specified coefficients and confidence intervals for a manyglm object
# Confidence intervals are calculated by assuming normality (to be fixed)

# Julian Byrnes January 2018
# julian.francis.byrnes@gmail.com

#' @title  Plots the coefficients of the covariates of a manyglm object with confidence intervals.
#' @description A way to plot the coefficients of the covariates of a manyglm object. Modifies code from Niku, Hui and Taskinen's coefplot.gllvm
#'
#' @param object A manyglm object
#' @param y.label Whether all the Y variables should be labelled
#' @param which.Xcoef Which X covariates should be included in the plot.
#' @param which.Ys Which Y variables should be included in the plot.
#' @param incl.intercept Whether the intercept coefficient should be included.
#' @param cex.ylab A plotting parameter. The default is 0.5.
#' @param mfrow Plotting parameter
#' @param mar Plotting parameter
#' @param ... Other plotting parameters
#' @return none
#' @export
coefplot.manyglm <- function(
    object,
    y.label = TRUE,
    which.Xcoef = NULL,
    which.Ys = NULL,
    incl.intercept = FALSE,
    cex.ylab = 0.5,
    mfrow = NULL,
    mar = NULL, ...){


  # Sanity checks
  if(is.null(object$coefficients))
    stop("There are no X covariates to plot.")
  if(is.null(object$y))
    stop("There are no Y values to plot.")

  # Check which coefficients to plot
  # Remove intercept if not included.
  if(is.null(which.Xcoef))
    if(!incl.intercept)
      which.Xcoef <- 2:nrow(object$coefficients)
    else
      which.Xcoef <- 1:nrow(object$coefficients)

  if(is.null(which.Ys))
    which.Ys <- 1:ncol(object$coefficients)


  coeffs <- as.data.frame(object$coefficients)[which.Xcoef,which.Ys]

  # Deal with any NAs by throwing an error
  if(anyNA(coeffs))
    stop("NAs not allowed in coefficient matrix.")

  numXvars <- nrow(coeffs)
  Xnames <- rownames(coeffs)
  numYs <- ncol(coeffs)
  Ynames <- colnames(coeffs)

  # Set plot settings
  if(missing(mfrow))
    mfrow <- c(1, numXvars)
  if (is.null(mar))
    original_par <- par(mfrow = mfrow)
  else
    original_par <- par(mfrow = mfrow, mar = mar)

  # For each X covariate, order the Ys from lowest to highest coefficients
  for (i in 1:numXvars) {
    plotCoeff <- sort(as.numeric(coeffs[i,]))
    plotCoefIndex <- order(as.numeric(coeffs[i,]))
    plotSD <- as.data.frame(object$stderr.coefficients)[which.Xcoef,which.Ys]
    plotSD <- plotSD[i,][plotCoefIndex]

    # Calculate (asympt normal) confidence intervals using standard error
    lowerCI <- unlist(plotCoeff - 1.96 * plotSD)
    upperCI <- unlist(plotCoeff + 1.96 * plotSD)

    # Determine whether they lie on one side of the 0 line
    # And colour differently accordingly
    col.seq <- rep("black", numYs)
    col.seq[lowerCI < 0 & upperCI > 0] <- "grey"

    # Plot
    plot(x = unlist(plotCoeff),
         y = 1:numYs,
         yaxt = "n",
         col=col.seq,
         xlab=Xnames[i],
         ylab="",
         xlim= c(min(lowerCI), max(upperCI)),
         pch = "x",
         cex.lab=1.3, ...)
    # More plot stuff
    segments(
      x0 = lowerCI,
      y0 = 1:numYs,
      x1 = upperCI,
      y1 = 1:numYs,
      col=col.seq)
    abline(v=0, lty=1)
    if(y.label)
      axis(2,
        at=1:numYs,
        labels=Ynames[plotCoefIndex],
        las=1,
        cex.axis=cex.ylab)
  }

  # gllvm

  # Checks for X covariates in model, exits if not.
  # Gets the names of the coefficients, length thereof
  # Additionally, gets the names of the species to which they are fitted

  # Default settings for coefplot - mfrow <- c(1, number of X coefficients)
  # & mar current setting
  # Otherwise can manually control the defaults

  # Then calculates confidence intervals and sorts
  # Values with confidence interval NOT intersecting the line are black; the rest are grey

  # reset par after the plot
  par(original_par)
}

coefplot <- function (object, ...) {
   UseMethod("coefplot")
}
