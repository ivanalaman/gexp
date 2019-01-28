gexp.lsd <- function(mu        = mu, 
                     err       = err,
                     factorsl  = factorsl,
                     rowsl     = rowsl,
                     colsl     = colsl,
                     ef        = ef,  
                     erow      = erow, 
                     ecol      = ecol, 
                     nrand     = nrand,    
                     contrasts = contrasts,
                     round     = round, 
                     random    = random)
{

  if(length(ef) != 1){
    stop('Use only one factor!')
  }

  f1 <- ef[[1]]
  n <- length(f1)

  if(is.null(factorsl)){
    levelss <- paste('x', 1:n, sep='')
  } else {
    if(length(factorsl)!=1){
      stop('Use only one factor!')
    }

    levelss <- unlist(factorsl)
  }

  ifelse(is.null(rowsl),
    levelsrows <- factor(rep(1:n, rep(n, n))),
    levelsrows <- factor(rep(unlist(rowsl), rep(n, n))))

  ifelse(is.null(colsl),
         levelscols <- factor(rep(1:n, n)),
         levelscols <- factor(rep(unlist(colsl), n)))

  aux_dados  <- latin(n,
                      levelss = levelss,
                      nrand = nrand)

  dados = data.frame(Row    = levelsrows,
                     Column = levelscols,
                     X1     = c(aux_dados))

  if(!is.null(rowsl)){
    names(dados) <- gsub('Row',
                         names(rowsl),
                         names(dados))
  }
  
  if(!is.null(colsl)){
    names(dados) <- gsub('Column',
                         names(colsl),
                         names(dados))
  }
   
  if(!is.null(factorsl)){
    names(dados) <- gsub('X1',
                         names(factorsl),
                         names(dados))
  }
  
  factors <- lapply(dados, levels)

  aux_X1 <- paste('~ ',
                  paste(names(dados),
                        collapse='+'))

  if(is.null(contrasts)){
    contrasts <- lapply(factors, function(x)diag(length(x)))
  } else{
    if((length(ef)+length(erow)+length(ecol)) != length(contrasts))stop('You must be include all the contrasts!')
  }

  names(contrasts) <- names(dados)

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

    betas <- as.matrix(c(mu, erow, ecol, f1)) 

  } else if(length(mu)!=0 & length(mu) > 1){

    betas <- rbind(mu, erow, ecol, f1) 

 } else if(is.null(mu) & length(ef) > 1 & all(unlist(lapply(factorsl,is.ordered))==TRUE)){#Todos os fatores são quantitativos e só há interesse em contrastes polinomiais

    aux_betas <- lapply(ef[-1],
                        function(x)x[-1])
    aux_betas1 <- c(ef[1],aux_betas)
    aux_betas2 <- lapply(aux_betas1,as.matrix)
    aux_betas3 <- do.call('rbind',aux_betas2)
    betas <- as.matrix(c(erow,ecol,aux_betas3))

  } else {

    aux_betas <- lapply(ef,as.matrix)
    aux_betas2 <- do.call('rbind',aux_betas)
    betas <- as.matrix(c(erow,ecol,aux_betas2))
    
  }
   
    #   } else {
    #     betas <- as.matrix(c(erow, ecol, f1))
    #   }

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
