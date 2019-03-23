gexp.rcbd <- function(mu        = mu,
                      err       = err,
                      r         = r,
                      fl        = fl,
                      blkl      = blkl,
                      fe        = fe,
                      blke      = blke,
                      contrasts = contrasts,
                      round     = round,
                      random    = random, ...)  
{
  if(is.null(fe))  fe <- list(f1 = rep(1,3))

  if(is.null(fl)){
    aux_factor <- lapply(fe,
                         function(x) as.matrix(x))

    names(aux_factor) <- paste('X', 1:length(fe), sep='')

    aux_factor1 <- as.list(tolower(names(aux_factor)))
    aux_factor2 <- lapply(aux_factor,
                          function(x) 1:dim(x)[1])

    factors <- mapply(function(x, y) paste(x, y, sep=''),
                      aux_factor1,
                      aux_factor2,
                      SIMPLIFY = FALSE)

    names(factors) <- names(aux_factor)
  } else {
    if(!is.list(fl)){
      stop('This argument must be a list. See examples!')
    }
    factors <- fl
  }

  factors$r <- 1:r

  ifelse(is.null(blkl),
         factors$Block <- 1:dim(as.matrix(blke))[1],
         factors[[names(blkl)]] <- unlist(blkl))

  dados <- expand.grid(factors,
                       KEEP.OUT.ATTRS = FALSE)

  aux_lf <- names(dados) 

  dados[[aux_lf[length(aux_lf)]]] <- factor(dados[[aux_lf[length(aux_lf)]]])

  lf <- aux_lf[-c((length(aux_lf)-1):length(aux_lf))] 
  aux_X1 <- paste('~ ',
                  aux_lf[length(aux_lf)],
                  '+',
                  paste(lf,
                        collapse='+')) 

  if(is.null(contrasts)){
    contrasts <- lapply(factors[aux_lf!='r'], function(x)diag(length(x)))
  } else{
    if((length(fe)+1) != length(contrasts))stop('You must be include all the contrasts!')
  }

  names(contrasts) <- names(dados)[aux_lf!='r']

  X  <- model.matrix(eval(parse(text=aux_X1)),
                     dados,
                     contrasts.arg=contrasts) 

  Z <- NULL

  if(is.null(err)){
   
    e <- mvtnorm::rmvnorm(n = dim(X)[1],  
                          sigma = diag(ncol(as.matrix(fe[[1]]))))
 
  } else {
    if(!is.matrix(err)) stop("This argument must be a matrix n x 1 univariate or n x p multivariate!")

    e <- err
  }

  if(length(mu)!=0 & length(mu) == 1){

    betas <- as.matrix(c(mu, blke, unlist(fe)))

  } else if(length(mu)!=0 & length(mu) > 1){

    betas <- rbind(mu, blke, do.call('rbind', fe))

  } else if(is.null(mu) & length(fe) > 1 & all(unlist(lapply(fl, is.ordered))==TRUE)){#Todos os fatores são quantitativos e só há interesse em contrastes polinomiais

    aux_betas <- lapply(fe[-1],
                        function(x)x[-1])
    aux_betas1 <- c(fe[1], aux_betas)
    aux_betas2 <- lapply(aux_betas1, as.matrix)
    aux_betas3 <- do.call('rbind', aux_betas2)
    betas <- as.matrix(c(blke, aux_betas3))

  } else {

    aux_betas <- lapply(fe, as.matrix)
    aux_betas2 <- do.call('rbind', aux_betas)
    betas <- as.matrix(c(blke, aux_betas2))
    
  }

  #} else {
  #  betas <- as.matrix(c(blke, unlist(fe)))
  #}

  yl <- X%*%betas + e

  colnames(yl) <- paste('Y', 1:dim(yl)[2], sep='') 

  Y <- round(yl, round)

  # J.C.Faria
  if(is.null(mu)){
    dados <- lapply(dados, 
                    function(x) if(is.ordered(factor(x))) as.numeric(as.character(x)) else x)

    dados <- as.data.frame(dados)                
  }                  

  dados <- cbind(dados, Y)

  if(random){

    dados <- dados[sample(dim(dados)[1]),]

  }

  res <- list(X   = X,
              Z   = Z,
              Y   = Y,
              dfm = dados)

}
