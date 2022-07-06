manyany = function(formula, fn, family="negative.binomial", data=NULL, composition = FALSE, block = NULL, get.what="details", var.power=NA, na.action = "na.exclude", ...)
{
  #MANYANY - applies a function of your choice to each column of YMAT and computes logLik by taxon.
  # FORMULA is the formula to use in the call to FM.
  # FN is a character vector giving the name of the function to be applied to each taxon.  e.g. "glm"
  # a FAMILY argument needs to be specified. It can be a list with different families for different variables (with length matching the number of columns in YMAT)
  # COMPOSITION is a logical switch indicating whether or not to do a compositional analysis (i.e., include a
  #  row effect in the model to account for changes in total abundance across samples).
  # VAR.POWER is needed for tweedie distributions - the power parameter needs to be specified here as well as in the family argument.
  # ... any further arguments required by FN.
  #
  # Examples:
  # require(mvabund)
  # data(spider)
  # abund=spider$abund
  # X=data.frame(spider$x)
  # 
  # a manygam:
  # library(mgcv)
  # ft=manyany(abund~s(soil.dry),"gam",data=X,family="poisson")
  # 
  # a manyglmm:
  # library(lme4)
  # gr = rep(1:2,each=14)
  # ft=manyany(abund~soil.dry+1|gr,"lmer",data=X,family="poisson")
  #
  ## A manyglm:
  # ft=manyany(abund~soil.dry,"glm",data=X,family="poisson")
  ## note this gives the same answer as:  
  # ft=manyglm(mvabund(abund)~X$soil,family="poisson")
  
  #Now compatible with ordinal data, but note that linear predictor is for the lowest
  #category only, whereas fitted is for the observed category.
  
  
  tol=1.e-8 #for truncation of linear predictor
  if(missing(data)) # Only coerce to model frame if not specified by the user.
    mf =model.frame(formula, parent.frame())
  else
    mf = model.frame(formula,data=data)

  # get response and its dimensions
  yMat = model.response(mf)
  yMat = as.matrix(yMat)
  n.rows = dim(yMat)[1]
  n.vars = dim(yMat)[2]
  
  # for clm, need to change family argument and make each column of response a factor
  if(fn=="clm")
  {
    allargs <- match.call(expand.dots = FALSE)
    dots <- allargs$...
    if( "link" %in% names(dots) )
      link <- dots$link
    else
      link="logit"
    if(link=="loglog")
      fam.i = binomial("cloglog") #to avoid errors since "loglog" is not defined
    else
      fam.i = binomial(link) #although not binomial
    fam.i$family = "ordinal"
    fam = vector(mode="list",length=n.vars)
    family=fam
    # convert data to a dataframe of factors:
    yMat = data.frame(yMat, stringsAsFactors=TRUE) #converting to data frame so factor input is read as factors
    for(iVar in 1:n.vars)
      yMat[,iVar] = as.factor(yMat[,iVar])
  }

  # get names for response, or assign if empty
  yNames = dimnames(yMat)
  if(is.null(yNames[[1]]))
    yNames[[1]] = 1:n.rows
  if(length(yNames)==1)
    yNames[[2]] = paste("y",1:n.vars,sep="")
  #    yNames[[2]] = "y" #to avoid issues later.
  
  call=match.call()
  
  # fix this bit later  
  if(composition==TRUE)
  {
    yVec    = as.vector(yMat)
    rows     = factor(rep(1:n.rows,n.vars))
    cols     = factor(rep(1:n.vars,each=n.rows))
    mf    = data.frame(mf[rows,],rows,cols)
    formula = formula(paste(formula[2],"~rows+cols+cols:(",formula[3],")",sep=""))
    n.rows.orig = n.rows #save for later
    n.vars.orig = n.vars #save for later
    n.rows  = length(yVec)
    n.vars  = 1
    names(yVec) = paste( yNames[[2]][cols], ".", yNames[[1]][rows], sep="")
    yMat    = as.matrix(yVec)
    if(inherits(family,"family")==FALSE & length(family)>1)
      stop("when using composition=TRUE, family argument must have length one.")
    if(is.null(block))  #to make sure resampling is by row of original data, not of vectorised data.
       block = rows
    else
       block = block[rows]
  }

  #If family is specified once, turn it into a list of n.vars family objects
  if(inherits(family,"family") || length(family)==1)
  {
    fam = family #temporary store to slot into a big list
    family = vector("list",n.vars)
    for(i.var in 1:n.vars)
      family[[i.var]] = fam
  }
  if(length(family)!=n.vars)
      stop("family argument has length more than one but not equal to the number of columns in yMat (!?)")

  if(length(var.power)==1)
    var.power=rep(var.power,n.vars)
    
  fam = family

  # now ensure each family is a proper family function
  for(i.var in 1:n.vars)
  {
    if (is.character(family[[i.var]])) 
    {
      if (family[[i.var]] == "negbinomial" || family[[i.var]]=="negative.binomial")
      {
        fam[[i.var]] = negative.binomial(10^6)
        fam[[i.var]]$family = family[[i.var]]
      }
      else if (family[[i.var]] == "binomial(link=logit)")
      {
        fam[[i.var]] = binomial()
        fam[[i.var]]$family = family[[i.var]]
      }
      else if (family[[i.var]] == "binomial(link=cloglog)" || family[[i.var]] == "cloglog")
      {
        fam[[i.var]] = binomial("cloglog")
        fam[[i.var]]$family = family[[i.var]]
      }
      else
      {
        fam.fn = get(fam[[i.var]], mode = "function", envir = parent.frame())        
        fam[[i.var]] = fam.fn()
      }  
    }
    if(fn=="clm")
      fam[[i.var]] = family[[i.var]] = fam.i
    if(fam[[i.var]]$family=="binomial")
      warning("The binomial option of manyany currently assumes you have binary (presence/absence) response")
  }

  # find response variable in mf (should be first but better safe than sorry)
  nameOfResponse  = as.character(formula[[2]])
  whichIsResponse = which(names(mf)==nameOfResponse)

    # set up empty objects
  manyfit = vector(mode="list",length=n.vars)
  fits = matrix(NA,n.rows,n.vars)
  etas = matrix(NA,n.rows,n.vars)
  params = manyfit
  logL = rep(NA,n.vars)
  
  # fit model sequentially for each variable
  for(i.var in 1:n.vars)
  {
    # change response to just column iVar
    mf[[1]] = yMat[,i.var]
    
    # refit model via do.call
    manyfit[[i.var]] = do.call(fn, list(formula=formula, family=family[[i.var]], data=mf, na.action=na.action, ...)) #note use of family argument as originally specified
    
    # store logL, or get from dviance if undefined
    logL[i.var]  = logLik(manyfit[[i.var]])
    if(is.na(logL[i.var]))
      logL[i.var] = -0.5*deviance(manyfit[[i.var]]) #just in case logL function is undefined, e.g. tweedie 

    # only get extra stuff if get.what says to... this is skipped by anova.manyany
    if(get.what=="details"||get.what=="models")
    {
            fits[,i.var] = fitted(manyfit[[i.var]])

      etas[,i.var] = switch(fn,
                            "lmer"=manyfit[[i.var]]@eta,
                            "clm"=predict(manyfit[[i.var]],type="linear.predictor",newdata=mf)$eta1,
                            predict(manyfit[[i.var]])
                         )
      #need to then truncate as if on logit scale...
      if(substr(fam[[i.var]]$family,1,3)=="bin" || fam[[i.var]]$family=="ordinal") #truncate linear predictor to more reasonable range
      {
        etas[,i.var] = pmax(etas[,i.var], fam[[i.var]]$linkfun(tol)/2)
        etas[,i.var] = pmin(etas[,i.var], fam[[i.var]]$linkfun(1-tol)/2)
      }
      if(fam[[i.var]]$link=="log"||fam[[i.var]]$link=="mu^0") #truncate linear predictor to more reasonable range
        etas[,i.var] = pmax(etas[,i.var], log(tol)/2)
      if(i.var==1)
      {
        cf = try(coef(manyfit[[i.var]]),silent=TRUE) #don't know if this function is defined
        if(inherits(cf, "try-error"))
        {
          do.coef   = FALSE
          coefs     = NULL
        } 
        else
        {
          coefs     = vector(mode="list",n.vars)
          coefs[[1]] = cf
          names(coefs[[1]])=dimnames(cf)[[1]]
          if(composition==FALSE & n.vars>1) #only name y variable if not compositional model
            names(coefs)=yNames[[2]]
          do.coef   = TRUE
        }
      }
      else
      {
        if(do.coef==TRUE)
          coefs[[i.var]] = coef(manyfit[[i.var]])
      }
      if(fam[[i.var]]$family=="poisson")
        params[[i.var]] = list(q=yMat[,i.var],lambda=fits[,i.var])
      if(substr(fam[[i.var]]$family,1,3)=="bin")
        params[[i.var]] = list(q=yMat[,i.var],prob=fits[,i.var],size=1)
      if(fam[[i.var]]$family=="Tweedie")
        params[[i.var]] = list(q=yMat[,i.var], power=var.power[i.var], mu=fits[,i.var], phi=summary(manyfit[[i.var]])$disp)
      if(fam[[i.var]]$family=="ordinal")
        params[[i.var]] = list(q=yMat[,i.var], mu=predict(manyfit[[i.var]], type="cum.prob"), muAll=predict(manyfit[[i.var]],type="cum.prob",newdata=mf[-whichIsResponse])$cprob2)
      if(grepl("egative",fam[[i.var]]$family) || fam[[i.var]]$family == "negbinomial")
      {
        if(any(names(manyfit[[i.var]])=="theta"))
          theta=manyfit[[i.var]]$theta
        else
          {
            if(any(names(manyfit[[i.var]])=="phi"))
              theta = 1/manyfit[[i.var]]$phi
            else # otherwise it must be fixed and tied up in the family argument 
              theta = 1/(fam[[i.var]]$var(1)-1)
          }
        params[[i.var]] = list(q=yMat[,i.var],mu=fits[,i.var],size=theta)
      }
      if(fam[[i.var]]$family=="gaussian")
      {
        s.ft=summary(manyfit[[i.var]])
        if(any(names(s.ft)=="sigma"))
          sd=s.ft$sigma
        else
          sd=s.ft$scale
        params[[i.var]] = list(q=yMat[,i.var],mean=fits[,i.var],sd=sd)
      }
    } #end get.what if statement
  } #end i.var loop
  object=list(logL=logL, get.what=get.what)

  #now format predictions and get residuals, if required
  if(get.what=="details"||get.what=="models")
  {
    if(composition==TRUE) #reshape to original data size if required
    {
      fits   = matrix(fits,n.rows.orig,n.vars.orig)
      etas   = matrix(etas,n.rows.orig,n.vars.orig)
    }    
    else #only name y variables for logL and params if not compositional model
    {
      if(n.vars>1)
      {
        names(logL) = yNames[[2]]
        names(params) = yNames[[2]]
      }
    }
    attributes(logL)$df = attributes(logLik(manyfit[[i.var]]))$df
    attributes(logL)$nobs = n.rows
    class(logL) = "logLik"
    resids = residuals.manyany(list(params=params, family=fam, composition=composition, fitted.values=fits, get.what=get.what))
    dimnames(resids) = yNames
    dimnames(fits)   = yNames
    dimnames(etas)   = yNames
    mf[[1]] = yMat #DW, 3/2/22 change: return full response
    object=list(logL=logL,fitted.values=fits,residuals=resids,linear.predictor=etas,family=fam, coefficients = coefs, call=call,params=params,model=mf, terms = terms(manyfit[[i.var]]), formula=formula, block=block, composition=composition, get.what=get.what)
#    object=list(logL=logL,fitted.values=fits,residuals=resids,linear.predictor=etas,family=fam, coefficients = coefs, call=call,params=params,model=model.frame(manyfit[[i.var]]), terms = terms(manyfit[[i.var]]), formula=formula, block=block, composition=composition, get.what=get.what)
  }
  if(get.what=="models") #also return the model fits, if requested
  {
    object$fits = manyfit
    names(object$fits) = yNames[[2]]
  }
  class(object)=c("manyany", class(manyfit[[i.var]]) )
  return(object)
} #end manyany function


print.manyany <- function(x, digits = max(3L, getOption("digits") - 3L), ...)
{
  n.vars = dim(x$fitted)[2]
  cat("\nCall:  ", paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n\n", sep = "")
  cat("Number of rows:\t   ", dim(x$fitted)[1], "\t Number of columns:\t   ", n.vars)
  cat("\n")
  cat("Number of parameters in model:\t   ", length(unlist(x)))
  cat("\n\n")
  cat("Residual Deviance:\t   ", format(signif(-2*sum(x$logL), digits)))
  cat("\n\n")
  cat("Dunn-Smyth Residuals:")
  cat("\n")
  cat(summary(qnorm(as.vector(x$residuals))))
  cat("\n\n-")
}
  
  
logLik.manyany <- function(object, ...)
{
  object$logL
}

coef.manyany <- function(object, ...)
{
  object$coefficients
}

residuals.manyany<- function(object, ...)
{
  if(object$get.what!="details" & object$get.what!="models")
    stop("To compute residuals, set get.what='details' in your manyany call")
  tol=1.e-6
  params = object$params
  n.rows = length(params[[1]]$q)
  n.vars = length(params)
  if(length(object$family)==1)
    family = rep(object$family,n.vars)
  else
    family=object$family
  resids=matrix(NA,n.rows,n.vars)
  dimnames(resids)[[1]] = names(params[[1]]$yMat)
  dimnames(resids)[[2]] = names(params)
  for(i.var in 1:n.vars)
  {
    if(family[[i.var]]$family=="ordinal")
    {
      u = runif(n.rows)
      resids[,i.var] = u*params[[i.var]]$mu$cprob1 + (1-u)*params[[i.var]]$mu$cprob2      
    }
    else
    {
      param.minus = params[[i.var]]
      param.minus$q = params[[i.var]]$q - tol
      if(grepl("egative",family[[i.var]]$family) || family[[i.var]]$family == "negbinomial")
        pfn = "pnbinom"
      if(family[[i.var]]$family=="poisson")
        pfn = "ppois"
      if(substr(family[[i.var]]$family,1,3)=="bin")
        pfn = "pbinom"
      if(family[[i.var]]$family=="gaussian")
      {
        pfn = "pnorm"
        param.minus$q = params[[i.var]]$q
      } 
      if(family[[i.var]]$family=="Tweedie")
        pfn = "ptweedie"
      u = runif(n.rows)
      #to avoid any values identically 1:
      pMinus = pmin(do.call(pfn, param.minus), 1-tol)
      resids[,i.var] = u*do.call(pfn, params[[i.var]]) + (1-u)*pMinus
    }
  }
  if(object$composition==TRUE) #reshape to original data size if required
    resids = matrix(resids, dim(object$fitted)[1], dim(object$fitted)[2])
  resids=qnorm(resids)
  return(resids)
}


plot.manyany=function(x, ...)
{
  object = x
  if(object$get.what!="details" & object$get.what!="models")
    stop("To plot your fit, set get.what='details' in your manyany call")
  
  # DW, 1/3/19: removed qnorm from next call, now done in resids function 
  Dunn.Smyth.Residuals=residuals.manyany(object)

  Fitted.values = object$linear
  
  # add colours if not already there, rainbow sorted by total abundance...
  if(hasArg("col")==F)
  {
    n.rows = dim(object$fitt)[1]
    n.vars = dim(object$fitt)[2]
    sumRank = rank(apply(object$linear,2,sum),ties.method="first")
    col=rep(rainbow(n.vars+1)[sumRank],each=n.rows)
    plot(Dunn.Smyth.Residuals~Fitted.values, ann=F, col=col, ...)
  }
  else
    plot(Dunn.Smyth.Residuals~Fitted.values, ann=F, ...)
  
  #add xlab and ylab if not already provided...
  args=match.call()
  mx=match("xlab",names(args),0L)
  if(mx==0)
    xlab="Linear Predictor"
  else
    xlab=deparse(args[mx])

  my=match("ylab",names(args),0L)
  if(my==0)
    ylab="Dunn-Smyth Residuals"
  else
    ylab=deparse(args[my])
  
  title(xlab=xlab,ylab=ylab)

} 



