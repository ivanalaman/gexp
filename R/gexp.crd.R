gexp.crd <- function(mu        = mu,
                     err       = err,
                     r         = r,
                     factorsl  = factorsl,
                     ef        = ef,
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

  dados <- expand.grid(factors,
                       KEEP.OUT.ATTRS = FALSE)

  aux_X1 <- paste('~',
                  paste(names(dados)[-length(names(dados))],
                        collapse='+')) 


  if(is.null(contrasts)){
    contrasts <- lapply(factors[1:length(ef)], function(x)diag(length(x)))
  }

  names(contrasts) <- names(dados)[1:length(ef)]

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

    betas <- as.matrix(c(mu, unlist(ef))) 

  } else if(length(mu)!=0 & length(mu) > 1){

    betas <- rbind(mu, do.call('rbind', ef)) 

  } else if(is.null(mu) & length(ef) > 1 & all(unlist(lapply(factorsl,is.ordered))==TRUE)){#Todos os fatores são quantitativos e só há interesse em contrastes polinomiais

    aux_betas <- lapply(ef[-1],
                        function(x)x[-1])
    aux_betas1 <- c(ef[1],aux_betas)
    aux_betas2 <- lapply(aux_betas1,as.matrix)
    betas <- do.call('rbind',aux_betas2)

   } else {

    aux_betas <- lapply(ef,as.matrix)
    betas <- do.call('rbind',aux_betas)
  
  }

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
