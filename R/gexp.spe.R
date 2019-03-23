gexp.spe <- function(mu        = mu,
                     err       = err,
                     errp      = errp,
                     r         = r, 
                     fl        = fl,
                     blkl      = blkl,
                     rowl      = rowl,
                     coll      = coll,
                     fe        = fe,
                     inte      = inte,
                     blke      = blke,
                     rowe      = rowe,
                     cole      = cole,
                     contrasts = contrasts,
                     nrand     = nrand,    
                     round     = round,   
                     random    = random, ...) 
{
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

  if(is.null(blke) & is.null(rowe) & is.null(cole)){ #é um DIC

    factors$r <- factor(1:r)

    dados <- expand.grid(factors,
                         KEEP.OUT.ATTRS = FALSE)

    aux_X1 <- paste('~',
                    paste(names(dados)[-length(names(dados))],
                          collapse='*')) 

    if(is.null(contrasts)){
      contrasts <- lapply(factors[1:length(fe)], function(x)diag(length(x)))
    }

    names(contrasts) <- names(dados)[1:length(fe)]

    X  <- model.matrix(eval(parse(text=aux_X1)),
                       dados,
                       contrasts.arg=contrasts) 


    parceprincipal <- names(factors)[1] 

    aux_Z1 <- paste('~ 0 + ',
                    parceprincipal,
                    ':r',
                    sep='')

    auxfactor <- c(parceprincipal, 'r')
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

    if(is.null(err)){

      e <- mvtnorm::rmvnorm(n = dim(X)[1],  
                            sigma = diag(ncol(as.matrix(fe[[1]]))))

    } else {
      if(!is.matrix(err)) stop("This argument must be a matrix n x 1 univariate or n x p multivariate!")

      e <- err
    }   

    if(is.null(errp)){
      e_plot <- mvtnorm::rmvnorm(n = dim(Z)[2],  
                                 sigma = diag(ncol(as.matrix(fe[[1]]))))
    } else {
      if(!is.matrix(errp)) stop("This argument must be a matrix n x 1 univariate or n x p multivariate!")

      e_plot <- errp
    }   

    #+ Help to interaction effects!
    auxl <- lapply(fe, length)
    niv <- unlist(auxl)
    n <- length(niv)
    resl <- list()
    for(i in 1:(n-1))resl[[i]] <- combn(niv, i+1)
    resl1 <- lapply(resl, function(x)apply(x, 2, prod))
    resl2 <- sum(unlist(resl1))

    if(length(inte)!=resl2){
      tcltk::tkmessageBox(message=paste('The number of effects this argument is ', resl2),
                          icon = 'error')
    }
    #End

    if(length(mu)!=0 & length(mu) == 1){

      betas <- as.matrix(c(mu, unlist(fe), inte))

    } else if(length(mu)!=0 & length(mu) > 1){

      betas <- rbind(mu, do.call('rbind', fe), inte)

    } else if(is.null(mu) & length(fe) > 1 & all(unlist(lapply(fl, is.ordered))==TRUE)){#Todos os fatores são quantitativos e só há interesse em contrastes polinomiais

      aux_betas <- lapply(fe[-1],
                          function(x)x[-1])
      aux_betas1 <- c(fe[1], aux_betas)
      aux_betas2 <- lapply(aux_betas1, as.matrix)
      aux_betas3 <- do.call('rbind', aux_betas2)
      betas <- as.matrix(c(aux_betas3, inte))

    } else {

      aux_betas <- lapply(fe, as.matrix)
      aux_betas2 <- do.call('rbind', aux_betas)
      betas <- as.matrix(c(aux_betas2, inte))

    }

    #     } else {
    # 
    #       betas <- as.matrix(c(unlist(fe), inte))
    # 
    #     }
    # 
    yl <- X%*%betas + Z%*%e_plot + e

    colnames(yl) <- paste('Y', 1:dim(yl)[2], sep='')

    Y <- round(yl, round)

  } else if(!is.null(blke) & is.null(rowe) & is.null(cole)){#é um DBC

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
                          collapse='*')) 

    if(is.null(contrasts)){
      contrasts <- lapply(factors[aux_lf!='r'], function(x)diag(length(x)))
    } else{
      if((length(fe)+length(blke)) != length(contrasts))stop('You must be include all the contrasts!')
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

    if(is.null(err)){

      e <- mvtnorm::rmvnorm(n = dim(X)[1],  
                            sigma = diag(ncol(as.matrix(fe[[1]]))))

    } else {
      if(!is.matrix(err)) stop("This argument must be a matrix n x 1 univariate or n x p multivariate!")

      e <- err

    }   

    if(is.null(errp)){

      e_plot <- mvtnorm::rmvnorm(n = dim(Z)[2],  
                                 sigma = diag(ncol(as.matrix(fe[[1]]))))

    } else {
      if(!is.matrix(errp)) stop("This argument must be a matrix n x 1 univariate or n x p multivariate!")

      e_plot <- errp

    }   

    #+ Help to interaction effects!
    auxl <- lapply(fe, length)
    niv <- unlist(auxl)
    n <- length(niv)
    resl <- list()
    for(i in 1:(n-1))resl[[i]] <- combn(niv, i+1)
    resl1 <- lapply(resl, function(x)apply(x, 2, prod))
    resl2 <- sum(unlist(resl1))

    if(length(inte)!=resl2){
      tcltk::tkmessageBox(message=paste('The number of effects this argument is ', resl2),
                          icon = 'error')
    }
    #End

    if(length(mu)!=0 & length(mu) == 1){

      betas <- as.matrix(c(mu, blke, unlist(fe), inte))

    } else if(length(mu)!=0 & length(mu) > 1){

      betas <- rbind(mu, blke, do.call('rbind', fe), inte)

    } else if(is.null(mu) & length(fe) > 1 & all(unlist(lapply(fl, is.ordered))==TRUE)){#Todos os fatores são quantitativos e só há interesse em contrastes polinomiais

      aux_betas <- lapply(fe[-1],
                          function(x)x[-1])
      aux_betas1 <- c(fe[1], aux_betas)
      aux_betas2 <- lapply(aux_betas1, as.matrix)
      aux_betas3 <- do.call('rbind', aux_betas2)
      betas <- as.matrix(c(blke, aux_betas3, inte))

    } else {

      aux_betas <- lapply(fe, as.matrix)
      aux_betas2 <- do.call('rbind', aux_betas)
      betas <- as.matrix(c(blke, aux_betas2, inte))

    }

    #     } else {
    # 
    #       betas <- as.matrix(c(blke, unlist(fe), inte))
    # 
    #     }
    # 
    yl <- X%*%betas + Z%*%e_plot + e

    colnames(yl) <- paste('Y', 1:dim(yl)[2], sep='')

    Y <- round(yl, round)

  } else { #é um DQL

    if(dim(as.matrix(fe[[1]]))[1] <= 2) stop("The minimum effect this factors is three.")

    aux_trats1 <- suppressWarnings(do.call('interaction', factors))
    aux_trats2 <- sort(levels(aux_trats1))

    nsubplot <- length(factors[[2]])

    aux_trats3 <- matrix(aux_trats2,
                         ncol=nsubplot,
                         byrow=TRUE)

    aux_trats4 <- apply(aux_trats3, 1, function(x)paste(x, collapse=' '))
    n <- length(factors[[1]])# deve ser do mesmo comprimento do fator no primeiro slot da lista (parcela)

    aux_trats5 <- latin(n, 
                        levelss = aux_trats4,
                        nrand = nrand)

    aux_trats7 <- as.matrix(c(aux_trats5))
    aux_trats8 <- apply(aux_trats7, 1, function(x)unlist(strsplit(x, ' ')))
    aux_trats9 <- as.matrix(c(aux_trats8))
    aux_trats10 <- strsplit(as.character(aux_trats9[, 1]), '[.]')
    trats <- do.call('rbind', aux_trats10)
    colnames(trats) <- names(factors)

    aux_column <- rep(1:n, rep(nsubplot, n))

    ifelse(is.null(rowl),
           levelsrows <- factor(rep(1:n, rep(length(aux_column), n))),
           levelsrows <- factor(rep(unlist(rowl), rep(length(aux_column), n))))

    ifelse(is.null(coll),
           levelscols <- factor(rep(rep(1:n, rep(nsubplot, n)), n)),
           levelscols <- factor(rep(rep(unlist(coll), rep(nsubplot, n)), n)))

    dados  <- data.frame(Row    = levelsrows,
                         Column = levelscols,
                         trats)

    if(!is.null(rowl)){
      names(dados) <- gsub('Row',
                           names(rowl),
                           names(dados))
    }

    if(!is.null(coll)){
      names(dados) <- gsub('Column',
                           names(coll),
                           names(dados))
    }

    aux_X1 <- paste('~ ',
                    paste(names(dados)[1:2],
                          collapse = '+'),
                    '+',
                    paste(names(dados)[-c(1:2)],
                          collapse='*')) 

    aufactors <- lapply(dados, levels)  

    if(is.null(contrasts)){
      contrasts <- lapply(aufactors,
                          function(x)diag(length(x)))
    } else {
      if((length(fe)+length(rowe)+length(cole)) != length(contrasts))stop('You must be include all the contrasts!')
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

    auxfactor <- c(parceprincipal, names(dados)[1:2])
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
    if(is.null(err)){

      e <- mvtnorm::rmvnorm(n = dim(X)[1],  
                            sigma = diag(ncol(as.matrix(fe[[1]]))))

    } else {
      if(!is.matrix(err)) stop("This argument must be a matrix n x 1 univariate or n x p multivariate!")

      e <- err

    }   

    if(is.null(errp)){

      e_plot <- mvtnorm::rmvnorm(n = dim(Z)[2],  
                                 sigma = diag(ncol(as.matrix(fe[[1]]))))

    } else {
      if(!is.matrix(errp)) stop("This argument must be a matrix n x 1 univariate or n x p multivariate!")

      e_plot <- errp

    }   

    #+ Help to interaction effects!
    auxl <- lapply(fe, length)
    niv <- unlist(auxl)
    n <- length(niv)
    resl <- list()
    for(i in 1:(n-1))resl[[i]] <- combn(niv, i+1)
    resl1 <- lapply(resl, function(x)apply(x, 2, prod))
    resl2 <- sum(unlist(resl1))

    if(length(inte)!=resl2){
      tcltk::tkmessageBox(message=paste('The number of effects this argument is ', resl2),
                          icon = 'error')
    }

    if(length(mu)!=0 & length(mu) == 1){

      betas <- as.matrix(c(mu, rowe, cole, unlist(fe), inte))

    } else if(length(mu)!=0 & length(mu) > 1){

      betas <- rbind(mu, rowe, cole, do.call('rbind', fe), inte)

    } else if(is.null(mu) & length(fe) > 1 & all(unlist(lapply(fl, is.ordered))==TRUE)){#Todos os fatores são quantitativos e só há interesse em contrastes polinomiais

      aux_betas <- lapply(fe[-1],
                          function(x)x[-1])
      aux_betas1 <- c(fe[1], aux_betas)
      aux_betas2 <- lapply(aux_betas1, as.matrix)
      aux_betas3 <- do.call('rbind', aux_betas2)
      betas <- as.matrix(c(rowe, cole, aux_betas3, inte))

    } else {

      aux_betas <- lapply(fe, as.matrix)
      aux_betas2 <- do.call('rbind', aux_betas)
      betas <- as.matrix(c(rowe, cole, aux_betas2, inte))

    }

    #     } else {
    #       betas <- as.matrix(c(rowe, cole, do.call('rbind', fe), inte))
    #     }
    # 
    yl <- X%*%betas + Z%*%e_plot + e

    colnames(yl) <- paste('Y', 1:dim(yl)[2], sep='')

    Y <- round(yl, round)
  }

  # J.C.Faria
  if(is.null(mu)){
    dados <- lapply(dados, 
                    function(x) if(is.ordered(factor(x))) as.numeric(as.character(x)) else x)

    dados <- as.data.frame(dados)                
  }                  

  dados <- cbind(dados, Y)

  if(random){

    dados <- dados[sample(dim(dados)[1]), ]

  }

  res <- list(X   = X,
              Z   = Z,
              Y   = Y,
              dfm = dados)

}
