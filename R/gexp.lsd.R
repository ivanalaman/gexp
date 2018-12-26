gexp.lsd <- function(mu           = mu, 
                     error        = error,
                     labelfactors = labelfactors,
                     labellsdrows = labellsdrows,
                     labellsdcols = labellsdcols,
                     ef           = ef,  
                     erow         = erow, 
                     ecol         = ecol, 
                     nrand        = nrand,    
                     contrasts    = contrasts,
                     rd           = rd, 
                     randomized   = randomized)
{

  if(length(ef) != 1){
    stop('Use only one factor!')
  }

  f1 <- ef[[1]]
  n <- length(f1)

  if(is.null(labelfactors)){
    levelss <- paste('x',1:n,sep='')
  } else {
    if(length(labelfactors)!=1){
      stop('Use only one factor!')
    }

    levelss <- unlist(labelfactors)
  }

  ifelse(is.null(labellsdrows),
    levelsrows <- factor(rep(1:n,rep(n,n))),
    levelsrows <- factor(rep(unlist(labellsdrows),rep(n,n))))

  ifelse(is.null(labellsdcols),
         levelscols <- factor(rep(1:n,n)),
         levelscols <- factor(rep(unlist(labellsdcols),n)))

  aux_dados  <- latin(n,
                      levelss = levelss,
                      nrand = nrand)

  dados = data.frame(Row    = levelsrows,
                     Column = levelscols,
                     X1     = c(aux_dados))

  if(!is.null(labellsdrows)){
    names(dados) <- gsub('Row',
                         names(labellsdrows),
                         names(dados))
  }
  
  if(!is.null(labellsdcols)){
    names(dados) <- gsub('Column',
                         names(labellsdcols),
                         names(dados))
  }
   
  if(!is.null(labelfactors)){
    names(dados) <- gsub('X1',
                         names(labelfactors),
                         names(dados))
  }
  
  factors <- lapply(dados,levels)

  aux_X1 <- paste('~ ',
                  paste(names(dados),
                        collapse='+'))

  if(is.null(contrasts)){
    contrasts <- lapply(factors,function(x)diag(length(x)))
  } else{
    if((length(ef)+length(erow)+length(ecol)) != length(contrasts))stop('You must be include all the contrasts!')
  }

  names(contrasts) <- names(dados)

  X  <- model.matrix(eval(parse(text=aux_X1)),
                     dados,
                     contrasts.arg=contrasts) 

  Z <- NULL
  
  if(is.null(error)){
    e <- mvtnorm::rmvnorm(n = dim(X)[1],  
                          sigma = as.matrix(1))
  } else {
    if(!is.matrix(error)) stop("This argument must be a matrix n x 1 univariate or n x p multivariate!")

    e <- error
  } 
 
  if(length(mu)!=0 & length(mu) == 1){
    betas <- as.matrix(c(mu, erow, ecol, f1)) 
  } else if(length(mu)!=0 & length(mu) > 1){
    betas <- rbind(mu, erow, ecol, f1) 
  } else {
    betas <- as.matrix(c(erow,ecol,f1))
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
