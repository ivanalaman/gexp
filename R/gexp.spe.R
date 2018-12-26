gexp.spe <- function(mu           = mu,
                     error        = error,
                     errorp       = errorp,
                     r            = r, 
                     labelfactors = labelfactors,  
                     labelblocks  = labelblocks,
                     labellsdrows = labellsdrows,
                     labellsdcols = labellsdcols,
                     ef           = ef, 
                     einter       = einter, 
                     eb           = eb, 
                     erow         = erow, 
                     ecol         = ecol, 
                     contrasts    = contrasts,
                     nrand        = nrand,    
                     rd           = rd,   
                     randomized   = randomized)  
{
  if(is.null(ef)) stop("You must specify at least a factor") 
  if(is.null(labelfactors)){   
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
    if(!is.list(labelfactors)){
      stop('This argument must be a list. See examples!')
    }
    factors <- labelfactors   
  }

  if(is.null(eb) & is.null(erow) & is.null(ecol)){ #é um DIC

    factors$r <- factor(1:r)

    dados <- expand.grid(factors,
                         KEEP.OUT.ATTRS = FALSE)

    aux_X1 <- paste('~',
                    paste(names(dados)[-length(names(dados))],
                          collapse='*')) 

    if(is.null(contrasts)){
      contrasts <- lapply(factors[1:length(ef)],function(x)diag(length(x)))
    }

    names(contrasts) <- names(dados)[1:length(ef)]

    X  <- model.matrix(eval(parse(text=aux_X1)),
                       dados,
                       contrasts.arg=contrasts) 


    parceprincipal <- names(factors)[1] 

    aux_Z1 <- paste('~ 0 + ',
                    parceprincipal,
                    ':r',
                    sep='')

    auxfactor <- c(parceprincipal,'r')
    aux_Z2 <-  paste('list(',
                     paste(auxfactor,
                           paste('= contrasts(dados$',
                                 auxfactor,
                                 sep=''),
                           ',',
                           'contrasts = FALSE)',
                           collapse=','),
                     ')') 

    Z = model.matrix(eval(parse(text=aux_Z1)),
                     dados,
                     contrasts.arg = eval(parse(text=aux_Z2)))             

    if(is.null(error)){
      e <- mvtnorm::rmvnorm(n = dim(X)[1],  
                            sigma = as.matrix(1))
    } else {
      if(!is.matrix(error)) stop("This argument must be a matrix n x 1 univariate or n x p multivariate!")

      e <- error
    }   

    if(is.null(errorp)){
      e_plot <- mvtnorm::rmvnorm(n = dim(Z)[2],  
                            sigma = as.matrix(1))
    } else {
      if(!is.matrix(errorp)) stop("This argument must be a matrix n x 1 univariate or n x p multivariate!")

      e_plot <- errorp
    }   

    #+ Help to interaction effects!
    auxl <- lapply(ef,length)
    niv <- unlist(auxl)
    n <- length(niv)
    resl <- list()
    for(i in 1:(n-1))resl[[i]] <- combn(niv,i+1)
    resl1 <- lapply(resl,function(x)apply(x,2,prod))
    resl2 <- sum(unlist(resl1))

    if(length(einter)!=resl2){
      tcltk::tkmessageBox(message=paste('The number of effects this argument is ',resl2),
                          icon = 'error')
    }
    #End

    if(length(mu)!=0 & length(mu) == 1){
      betas <- as.matrix(c(mu, unlist(ef), einter)) 
    } else if(length(mu)!=0 & length(mu) > 1){
      betas <- rbind(mu, do.call('rbind',ef), einter) 
    } else {
      betas <- as.matrix(c(unlist(ef), einter))
    }

    yl <- X%*%betas + Z%*%e_plot + e

    colnames(yl) <- paste('Y',1:dim(yl)[2],sep='')

    Y <- round(yl, rd)

  } else if(!is.null(eb) & is.null(erow) & is.null(ecol)){#é um DBC

    factors$r <- 1:r  

    ifelse(is.null(labelblocks),
           factors$Block <- 1:dim(as.matrix(eb))[1],
           factors[[names(labelblocks)]] <- unlist(labelblocks))

    dados <- expand.grid(factors,
                         KEEP.OUT.ATTRS = FALSE)

    aux_lf <- names(dados) 

    dados[[aux_lf[length(aux_lf)]]] <- factor(dados[[aux_lf[length(aux_lf)]]])

    lf <- aux_lf[-c((length(aux_lf)-1):length(aux_lf))] 
    aux_X1 <- paste('~ ',
                    aux_lf[length(aux_lf)],
                    '+',
                    paste(lf,
                          collapse='*')) 

    if(is.null(contrasts)){
      contrasts <- lapply(factors[aux_lf!='r'],function(x)diag(length(x)))
    } else{
      if((length(ef)+length(eb)) != length(contrasts))stop('You must be include all the contrasts!')
    }

    names(contrasts) <- names(dados)[aux_lf!='r']

    X  <- model.matrix(eval(parse(text=aux_X1)),
                       dados,
                       contrasts.arg=contrasts) 

    aux_Z1 <-  paste('~ 0 + ',
                     paste(lf[1],
                           aux_lf[length(aux_lf)],
                           sep=':'))

    aux_Z2 <- paste('list(',
                    paste(c(lf[1],
                            aux_lf[length(aux_lf)]),
                          paste('= contrasts(dados$',
                                c(lf[1],
                                  aux_lf[length(aux_lf)]),
                                sep=''),
                          ',',
                          'contrasts = FALSE)',
                          collapse=','),
                    ')')

    Z <- model.matrix(eval(parse(text=aux_Z1)), 
                      dados, 
                      contrasts.arg = eval(parse(text=aux_Z2)))

    if(is.null(error)){
      e <- mvtnorm::rmvnorm(n = dim(X)[1],  
                            sigma = as.matrix(1))
    } else {
      if(!is.matrix(error)) stop("This argument must be a matrix n x 1 univariate or n x p multivariate!")

      e <- error
    }   

    if(is.null(errorp)){
      e_plot <- mvtnorm::rmvnorm(n = dim(Z)[2],  
                                 sigma = as.matrix(1))
    } else {
      if(!is.matrix(errorp)) stop("This argument must be a matrix n x 1 univariate or n x p multivariate!")

      e_plot <- errorp
    }   

    #+ Help to interaction effects!
    auxl <- lapply(ef,length)
    niv <- unlist(auxl)
    n <- length(niv)
    resl <- list()
    for(i in 1:(n-1))resl[[i]] <- combn(niv,i+1)
    resl1 <- lapply(resl,function(x)apply(x,2,prod))
    resl2 <- sum(unlist(resl1))

    if(length(einter)!=resl2){
      tcltk::tkmessageBox(message=paste('The number of effects this argument is ',resl2),
                          icon = 'error')
    }
    #End
    
    if(length(mu)!=0 & length(mu) == 1){
      betas <- as.matrix(c(mu, eb, unlist(ef), einter)) 
    } else if(length(mu)!=0 & length(mu) > 1){
      betas <- rbind(mu, eb, do.call('rbind',ef), einter) 
    } else {
      betas <- as.matrix(c(eb,unlist(ef), einter))
    }

    yl <- X%*%betas + Z%*%e_plot + e

    colnames(yl) <- paste('Y',1:dim(yl)[2],sep='')

    Y <- round(yl, rd)

  } else { #é um DQL

    if(dim(as.matrix(ef[[1]]))[1] <= 2) stop("The minimum effect this factors is three.")

    aux_trats1 <- suppressWarnings(do.call('interaction',factors))
    aux_trats2 <- sort(levels(aux_trats1))

    nsubplot <- length(factors[[2]])

    aux_trats3 <- matrix(aux_trats2,
                         ncol=nsubplot,
                         byrow=TRUE)

    aux_trats4 <- apply(aux_trats3,1,function(x)paste(x,collapse=' '))
    n <- length(factors[[1]])# deve ser do mesmo comprimento do fator no primeiro slot da lista (parcela)

    aux_trats5 <- latin(n, 
                        levelss = aux_trats4,
                        nrand = nrand)

    aux_trats7 <- as.matrix(c(aux_trats5))
    aux_trats8 <- apply(aux_trats7, 1, function(x)unlist(strsplit(x,' ')))
    aux_trats9 <- as.matrix(c(aux_trats8))
    aux_trats10 <- strsplit(as.character(aux_trats9[,1]),'[.]')
    trats <- do.call('rbind',aux_trats10)
    colnames(trats) <- names(factors)

    aux_column <- rep(1:n, rep(nsubplot,n))

    ifelse(is.null(labellsdrows),
           levelsrows <- factor(rep(1:n, rep(length(aux_column),n))),
           levelsrows <- factor(rep(unlist(labellsdrows),rep(length(aux_column),n))))

    ifelse(is.null(labellsdcols),
           levelscols <- factor(rep(rep(1:n, rep(nsubplot,n)),n)),
           levelscols <- factor(rep(rep(unlist(labellsdcols),rep(nsubplot,n)),n))) 
     
    dados  <- data.frame(Row    = levelsrows,
                         Column = levelscols,
                         trats)

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

    aux_X1 <- paste('~ ',
                    paste(names(dados)[1:2],
                          collapse = '+'),
                    '+',
                    paste(names(dados)[-c(1:2)],
                          collapse='*')) 

    aufactors <- lapply(dados,levels)  

    if(is.null(contrasts)){
      contrasts <- lapply(aufactors,
                          function(x)diag(length(x)))
    } else {
      if((length(ef)+length(erow)+length(ecol)) != length(contrasts))stop('You must be include all the contrasts!')
    }

    names(contrasts) <- names(dados)

    X  <- model.matrix(eval(parse(text=aux_X1)),
                       dados,
                       contrasts.arg=contrasts) 

    parceprincipal <- names(factors)[1]

    aux_Z1 <- paste('~ 0 + ',
                    parceprincipal,
                    ':',
                    paste(names(dados)[1:2],
                          collapse=':'),
                    sep='')

    auxfactor <- c(parceprincipal,names(dados)[1:2])
    aux_Z2 <-  paste('list(',
                     paste(auxfactor,
                           paste('= contrasts(dados$',
                                 auxfactor,
                                 sep=''),
                           ',',
                           'contrasts = FALSE)',
                           collapse=','),
                     ')') 

    Z = model.matrix(eval(parse(text=aux_Z1)),
                     dados,
                     contrasts.arg = eval(parse(text=aux_Z2)))             
    if(is.null(error)){
      e <- mvtnorm::rmvnorm(n = dim(X)[1],  
                            sigma = as.matrix(1))
    } else {
      if(!is.matrix(error)) stop("This argument must be a matrix n x 1 univariate or n x p multivariate!")

      e <- error
    }   

    if(is.null(errorp)){
      e_plot <- mvtnorm::rmvnorm(n = dim(Z)[2],  
                                 sigma = as.matrix(1))
    } else {
      if(!is.matrix(errorp)) stop("This argument must be a matrix n x 1 univariate or n x p multivariate!")

      e_plot <- errorp
    }   

    #+ Help to interaction effects!
    auxl <- lapply(ef,length)
    niv <- unlist(auxl)
    n <- length(niv)
    resl <- list()
    for(i in 1:(n-1))resl[[i]] <- combn(niv,i+1)
    resl1 <- lapply(resl,function(x)apply(x,2,prod))
    resl2 <- sum(unlist(resl1))

    if(length(einter)!=resl2){
      tcltk::tkmessageBox(message=paste('The number of effects this argument is ',resl2),
                          icon = 'error')
    }
   
    if(length(mu)!=0 & length(mu) == 1){
      betas <- as.matrix(c(mu, erow, ecol, unlist(ef), einter)) 
    } else if(length(mu)!=0 & length(mu) > 1){
      betas <- rbind(mu, erow, ecol, do.call('rbind',ef), einter) 
    } else {
      betas <- as.matrix(c(erow,ecol,do.call('rbind',ef), einter))
    }

    yl <- X%*%betas + Z%*%e_plot + e

    colnames(yl) <- paste('Y',1:dim(yl)[2],sep='')

    Y <- round(yl, rd)
  }

  dados <- cbind(dados,Y)

  if(randomized){

    dados <- dados[sample(dim(dados)[1]),]

  }

  res <- list(X   = X,
              Z   = Z,
              Y   = Y,
              dfm = dados)

}
