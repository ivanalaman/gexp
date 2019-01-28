gexp.rcbd <- function(mu        = mu,
                      err       = err,
                      r         = r,
                      factorsl  = factorsl,
                      blocksl   = blocksl,
                      ef        = ef,
                      eb        = eb,
                      contrasts = contrasts,
                      round     = round,
                      random    = random)  
{
  if(is.null(ef)) stop("You must specify at least a factor") 

  if(is.null(factorsl)){ 
    aux_factor <- lapply(ef,
                         function(x) as.matrix(x))

    names(aux_factor) <- paste('X', 1:length(ef), sep='')

    aux_factor1 <- as.list(tolower(names(aux_factor)))
    aux_factor2 <- lapply(aux_factor,
                          function(x) 1:dim(x)[1])

    factors <- mapply(function(x, y) paste(x, y, sep=''),
                      aux_factor1,
                      aux_factor2,
                      SIMPLIFY = FALSE)

    names(factors) <- names(aux_factor)
  } else {
    if(!is.list(factorsl)){
      stop('This argument must be a list. See examples!')
    }
    factors <- factorsl 
  }

  factors$r <- 1:r

  ifelse(is.null(blocksl),
         factors$Block <- 1:dim(as.matrix(eb))[1],
         factors[[names(blocksl)]] <- unlist(blocksl))

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
    if((length(ef)+1) != length(contrasts))stop('You must be include all the contrasts!')
  }

  names(contrasts) <- names(dados)[aux_lf!='r']

  X  <- model.matrix(eval(parse(text=aux_X1)),
                     dados,
                     contrasts.arg=contrasts) 

  Z <- NULL

  if(is.null(err)){
   
    e <- mvtnorm::rmvnorm(n = dim(X)[1],  
                          sigma = diag(ncol(as.matrix(ef[[1]])))) 
 
  } else {
    if(!is.matrix(err)) stop("This argument must be a matrix n x 1 univariate or n x p multivariate!")

    e <- err
  }

  if(length(mu)!=0 & length(mu) == 1){

    betas <- as.matrix(c(mu, eb, unlist(ef))) 

  } else if(length(mu)!=0 & length(mu) > 1){

    betas <- rbind(mu, eb, do.call('rbind', ef)) 

  } else if(is.null(mu) & length(ef) > 1 & all(unlist(lapply(factorsl,is.ordered))==TRUE)){#Todos os fatores são quantitativos e só há interesse em contrastes polinomiais

    aux_betas <- lapply(ef[-1],
                        function(x)x[-1])
    aux_betas1 <- c(ef[1],aux_betas)
    aux_betas2 <- lapply(aux_betas1,as.matrix)
    aux_betas3 <- do.call('rbind',aux_betas2)
    betas <- as.matrix(c(eb,aux_betas3))

  } else {

    aux_betas <- lapply(ef,as.matrix)
    aux_betas2 <- do.call('rbind',aux_betas)
    betas <- as.matrix(c(eb,aux_betas2))
    
  }

  #} else {
  #  betas <- as.matrix(c(eb, unlist(ef)))
  #}

  yl <- X%*%betas + e

  colnames(yl) <- paste('Y', 1:dim(yl)[2], sep='') 

  Y <- round(yl, round)

  dados <- cbind(dados, Y)

  if(random){

    dados <- dados[sample(dim(dados)[1]),]

  }

  res <- list(X   = X,
              Z   = Z,
              Y   = Y,
              dfm = dados)

}
