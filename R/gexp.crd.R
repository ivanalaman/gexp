gexp.crd <- function(mu         = mu,
                     sigma      = sigma,
                     r          = r,
                     ef         = ef,
                     contrasts  = contrasts,
                     rd         = rd,
                     newfactors = newfactors,
                     randomized = randomized)
{
  if(is.null(ef)) stop("You must specify at least a factor") 

  sigma <- as.matrix(sigma)

  if(is.null(newfactors)){
    aux_factor <- lapply(ef,
                         function(x) as.matrix(x))

    names(aux_factor) <- paste('X',1:length(ef),sep='')

    aux_factor1 <- as.list(tolower(names(aux_factor)))
    aux_factor2 <- lapply(aux_factor,
                          function(x) 1:dim(x)[1])

    factors <- mapply(function(x,y) paste(x,y,sep=''),
                      aux_factor1,
                      aux_factor2,
                      SIMPLIFY = FALSE)

    names(factors) <- names(aux_factor)
  } else {
    if(!is.list(newfactors)){
      stop('This argument must be a list. See examples!')
    }
    factors <- newfactors
  }

  factors$r <- 1:r

  dados <- expand.grid(factors,
                       KEEP.OUT.ATTRS = FALSE)

  aux_X1 <- paste('~',
                  paste(names(dados)[-length(names(dados))],
                        collapse='+')) 

  
  if(is.null(contrasts)){
    contrasts <- lapply(factors[1:length(ef)],function(x)diag(length(x)))
  }

  names(contrasts) <- names(dados)[1:length(ef)]

  X  <- model.matrix(eval(parse(text=aux_X1)),
                     dados,
                     contrasts.arg=contrasts)
  Z <- NULL

  e <- mvtnorm::rmvnorm(n = dim(X)[1],  
                        sigma = sigma)

  if(length(mu)!=0 & length(mu) == 1){
    betas <- as.matrix(c(mu, unlist(ef))) 
  } else if(length(mu)!=0 & length(mu) > 1){
    betas <- rbind(mu, do.call('rbind',ef)) 
  } else {
    betas <- as.matrix(unlist(ef))
  }

  yl <- X%*%betas + e

  colnames(yl) <- paste('Y',1:dim(yl)[2],sep='')

  Y <- round(yl, rd)

  dados <- cbind(dados,Y)

  if(randomized){

    dados <- dados[sample(dim(dados)[1]),]

  }

  res <- list(X   = X,
              Z   = Z,
              Y   = Y,
              dfm = dados)

}
