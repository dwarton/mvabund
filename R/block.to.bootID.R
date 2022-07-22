# used in the summary and anova functions of manyglm and manylm to get a matrix of 
# boot id's based off of a block, see ?anova.manyglm for more details
block_to_bootID <- function (block, bootID, nRows, nBoot, resamp) {
  tb = table(block)
  nLevels = length(tb)
  if (any(tb != nRows/nLevels) && resamp != "case") {
    print(tb)
    stop("Sorry, unless you are using case resampling, block needs to be a balanced factor - same number of rows for each level. Try using resamp='case'.")
  } else {
    blockIDs = vector("list",nLevels)
    for(i.level in 1:nLevels)
      blockIDs[[i.level]] = which(block==names(tb)[i.level])
    unlistIDs = unlist(blockIDs) # needed to match each resampled observation with its correct location
  }
  # then each iteration...
  # generate a bootID matrix if required
  if(is.null(bootID)){
    if(resamp=="perm.resid")
        samp <- t( sapply(rep(nLevels,nBoot),sample,replace=FALSE) )
    else
        samp <- matrix(sample(nLevels, nLevels * nBoot, replace=TRUE), ncol=nLevels)
  }
  else {
    samp <- bootID
  }
  bootID <-  matrix(NA,nBoot,nRows)
  for(iBoot in 1:nBoot) {
    # unlistIDs is needed to make sure each unlisted blockID ends up in the right place
    IDsboot = unlist(blockIDs[samp[iBoot, ]])
    if(length(IDsboot)<nRows) #redo it but this time ensure enough blocks
    {
      nIDs = sample(nLevels,nRows,replace=TRUE)
      IDsboot = unlist(blockIDs[nIDs])
    }
    bootID[iBoot, unlistIDs] = unlist(IDsboot[1:nRows])
  }
  bootID = bootID - 1 #to fit the format in C, 0 to nRows.
  if(interactive()) cat(paste("Using block resampling...","\n"))
  bootID
}
