################ get_polys for getting orthogonal polynomials ###################

get_polys = function( X, X.des.train=NULL)
{
  # get_polys will take a matrix of env vars (or trait vars), and standardise the quant ones
  # as well as return orthogonal poly's. Importantly, if training matrices are given as input,
  # these will be used in matrix construction.
  
  
  if(is.null(dim(X)))
    X = data.frame(X, stringsAsFactors=TRUE)
  n.sites = dim(X)[1]
  n.params  = dim(X)[2]
  if(is.null(X.des.train))
    n.train.sites = n.sites
  else
    n.train.sites <- dim(X.des.train$X)[1]
  if(is.null(X.des.train$var.type))
    var.type = rep("quantitative",n.params)
  else
    var.type = X.des.train$var.type
  for (i in 1:n.params)
  {
    
    # test if binary quantitative, if so, change to a factor to avoid poly error.  But only if training data
    if(is.null(X.des.train$var.type) & is.factor(X[,i])==FALSE)
    {
      testBinary = try(duplicated(X[,i],nmax=2), silent=TRUE)
      if(inherits(testBinary,"logical"))
      {
        X[,i] = factor(X[,i])
        warning(paste0("Binary variable '", names(X)[i], "' found and changed to a factor"))
      }
    }
    
    if(is.factor(X[,i]))
    {
      n.levels    = length(levels(X[,i]))
      if(n.levels==2)
      {
        var.type[i]="binary" #treat as quantitative but don't find its square
        #change variable name to indicate what it is doing
        dimnames(X)[[2]][i]=paste(names(X)[i],levels(X[,i])[2],sep="")
        #change entry in X to numerical
        X[,i]  = as.numeric(as.factor(X[,i]))*2-3
      }
      else
      {
        var.type[i]="factor"
        contrasts(X[,i])[contrasts(X[,i])==0] = -1
      }
    }
  }
  
  # to return standardised values of quant vars, with coeff, where needed:
  is.quant = which(var.type=="quantitative")
  n.quant = length(is.quant)
  if( n.quant>0 )
  {
    X.squ = X[,0]
    degs = c()
    names.X.squ = c()
    poly.coefs = as.list( rep( NA, n.quant ) )
    for(i.quant in is.quant)
    {
      poly.i = poly( X[,i.quant], degree=2, coefs=X.des.train$coefs[[i.quant]] )
      X.squ = cbind( X.squ, poly.i )
      degs  = c( degs, attr(poly.i, "degree") )
      poly.coefs[[ i.quant ]] = attr(poly.i, "coefs")   
      names.X.squ = c( names.X.squ, dimnames(X)[[2]][i.quant], paste( dimnames(X)[[2]][i.quant], ".squ", sep="") )
    }
    X.squ = X.squ * sqrt(n.train.sites)
    dimnames(X.squ)[[2]] = names.X.squ
    X[,var.type=="quantitative"] = X.squ[,degs==1]
    #get rid of the linear terms:
    X.squ = data.frame(X.squ[,degs==2])
  }
  else
  {
    X.squ=NULL
    poly.coefs=NULL
  }
  # to return orthogonal poly values of quant vars (with coeff):
  return( list( X=X, X.squ=X.squ, var.type=var.type, coefs=poly.coefs ) )
}
