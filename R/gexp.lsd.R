gexp.lsd <- function(mu         = mu, 
                     sigma      = sigma,
                     newfactors = newfactors,  
                     ef         = ef,  
                     erow       = erow, 
                     ecol       = ecol, 
                     nrand      = nrand,    
                     contrasts  = contrasts,
                     rd         = rd, 
                     randomized = randomized)
{

  if(length(ef) !=1){
    stop('Use only one factor!')
  }

  sigma <- as.matrix(sigma)

  f1 <- ef[[1]]
  n <- length(f1)

  if(is.null(newfactors)){
    levelss <- paste('x',1:n,sep='')
  } else {

    if(length(newfactors)!=1){
      stop('Use only one factor!')
    }

    levelss <- unlist(newfactors)

  }

  aux_dados  <- latin(n,
                      levelss = levelss,
                      nrand = nrand)

  dados = data.frame(Row    = factor(rep(1:n,rep(n,n))),
                     Column = factor(rep(1:n,n)),
                     X1   = c(aux_dados))

  if(!is.null(newfactors)){
    names(dados) <- gsub('X1',names(newfactors),names(dados))
  }

  factors <- lapply(dados,levels)
  #factors <- names(dados)

  aux_X1 <- paste('~ Row + Column +',
                  names(factors[3]))

  #   aux_X2 <- paste('list(',
  #                   paste(factors,
  #                         paste('= contrasts(dados$',
  #                               factors,
  #                               sep=''),
  #                         ',',
  #                         'contrasts = FALSE)',
  #                         collapse=','),
  #                   ')')
  #    X = model.matrix(eval(parse(text=aux_X1)),
  #                    dados,
  #                    contrasts.arg = eval(parse(text=aux_X2)))   
  #    

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

  e <- mvtnorm::rmvnorm(n = dim(X)[1],  
                        sigma = sigma)

  #   ifelse(length(mu) == 1,
  #          betas <- as.matrix(c(mu, erow, ecol, f1)),
  #          betas <- rbind(mu, erow, ecol, f1))

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
