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
  #-----------------------------------------------#
  #++++++++ Help to interaction effects! +++++++++#
  #-----------------------------------------------# 
  if(length(mu)==1){
    auxl <- lapply(fe,
                   length)
  } else {
    aux <- lapply(fe,function(x)apply(x,2,length))
    auxl <- lapply(aux,unique)
  }

  niv <- unlist(auxl) 

  n <- length(niv)
  resl <- list()
  for(i in 1:(n-1))resl[[i]] <- combn(niv, i+1)
  resl1 <- lapply(resl,
                  function(x) apply(x,
                                    2,
                                    prod))
  resl2 <- sum(unlist(resl1))

  if(length(mu)==1){
    if(is.null(inte)){
      inte <- rep(1,
                  resl2)
    }else{
      if(length(inte)!=resl2)    # J.C.Faria
        stop(paste("The length of the 'inte' argument must be: ", 
                   resl2,
                   sep='')) 
    }
  } else {
    if(is.null(inte)){
      inte <- matrix(rep(1,
                         resl2*length(mu)),
                     ncol=length(mu)) 
    } else {
      if(dim(inte)[1]!=resl2)    # J.C.Faria
        stop(paste("The length of the 'inte' argument must be: ", 
                   resl2,
                   sep=''))  
    }
  }

  #-----------------------------------------------#
  #++++++++ Help to interaction effects! +++++++++#
  #-----------------------------------------------#  

  if(is.null(fl)){
    aux_factor <- lapply(fe,
                         function(x) as.matrix(x))

    names(aux_factor) <- paste('X',
                               1:length(fe),
                               sep='')

    aux_factor1 <- as.list(tolower(names(aux_factor)))
    aux_factor2 <- lapply(aux_factor,
                          function(x) 1:dim(x)[1])

    factors <- mapply(function(x, y) paste(x, y, sep=''),
                      aux_factor1,
                      aux_factor2,
                      SIMPLIFY=FALSE)

    names(factors) <- names(aux_factor)

    quanti <- FALSE
    quali <- TRUE 

  } else {
    if(!is.list(fl)){
      stop('This argument must be a list. See examples!')
    }

    quanti <- all(lapply(fl,function(x)is.numeric(x)) == TRUE)
    quali  <- all(lapply(fl,function(x)is.numeric(x)) != TRUE)

    #Se não for nem quanti nem quali é pq é um híbrido!Neste caso, vamos encontrar as posições.              
    posquanti <- which(unlist(lapply(fl,is.numeric)) == TRUE)#em qual posição estão os quanti

    if(quanti){
      factors <- lapply(fl,as.ordered)
    } else if(!quanti & !quali){#híbrido
      factors <- fl
      factors[posquanti] <- lapply(fl[posquanti],as.ordered)
    } else {
      factors <- fl
    }  
  }

  if(is.null(blke) & is.null(rowe) & is.null(cole)){ #é um DIC

    factors$r <- factor(1:r)

    dados <- expand.grid(factors,
                         KEEP.OUT.ATTRS=FALSE)

    aux_X1 <- paste('~',
                    paste(names(dados)[-length(names(dados))],
                          collapse='*')) 

    #     if(is.null(contrasts)){
    #       contrasts <- lapply(factors[1:length(fe)],
    #                           function(x)diag(length(x)))
    #     }

    #names(contrasts) <- names(dados)[1:length(fe)]
    auxfactors <- factors[names(factors)!='r']
    if(quali){#Só qualitativos
      contrast <- lapply(auxfactors,
                         function(x)diag(length(x))) 
    } else if(quanti){#Só quantitativos
      contrast <- lapply(auxfactors,
                         function(x)contr.poly(length(x)))  
    } else {#híbrido
      contrast <- lapply(auxfactors,
                         function(x)diag(length(x)))
      contrast[posquanti] <- lapply(auxfactors[posquanti],
                                    function(x)contr.poly(length(x)))   
    }  

    if(!is.null(contrasts)){
      contrast[names(contrasts)] <- contrasts
      contrasts <- contrast
    }else{
      contrasts <- contrast
    }

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
                           'contrasts=FALSE)',
                           collapse=','),
                     ')') 

    Z <- model.matrix(eval(parse(text=aux_Z1)),
                      dados,
                      contrasts.arg=eval(parse(text=aux_Z2)))

    if(is.null(err)){

      e <- mvtnorm::rmvnorm(n=dim(X)[1],
                            sigma=diag(ncol(as.matrix(fe[[1]]))))

    } else {
      if(!is.matrix(err))
        stop("This argument must be a matrix n x 1 univariate or n x p multivariate!")

      e <- err
    }   

    if(is.null(errp)){
      e_plot <- mvtnorm::rmvnorm(n=dim(Z)[2],
                                 sigma=diag(ncol(as.matrix(fe[[1]]))))
    } else {
      if(!is.matrix(errp))
        stop("This argument must be a matrix n x 1 univariate or n x p multivariate!")

      e_plot <- errp
    }   

    #+ Help to interaction effects!
    #     auxl <- lapply(fe, length)
    #     niv <- unlist(auxl)
    #     n <- length(niv)
    #     resl <- list()
    #     for(i in 1:(n-1))resl[[i]] <- combn(niv, i+1)
    #     resl1 <- lapply(resl,
    #                     function(x) apply(x,
    #                                       2,
    #                                       prod))
    #     resl2 <- sum(unlist(resl1))
    # 
    #     if(length(inte)!=resl2){
    #       tcltk::tkmessageBox(message=paste('The number of effects this argument is ', resl2),
    #                           icon='error')
    #     }
    #End

    #     if(length(mu)!=0 & length(mu) == 1){
    # 
    #       betas <- as.matrix(c(mu,
    #                            unlist(fe),
    #                            inte))
    # 
    #     } else if(length(mu)!=0 & length(mu) > 1){
    # 
    #       betas <- rbind(mu,
    #                      do.call('rbind',
    #                              fe),
    #                      inte)
    # 
    #     } else if(is.null(mu) & length(fe) > 1 & all(unlist(lapply(fl, is.ordered)) == TRUE)) { # Todos os fatores são quantitativos e só há interesse em contrastes polinomiais
    # 
    #       aux_betas <- lapply(fe[-1],
    #                           function(x)x[-1])
    #       aux_betas1 <- c(fe[1],
    #                       aux_betas)
    #       aux_betas2 <- lapply(aux_betas1,
    #                            as.matrix)
    #       aux_betas3 <- do.call('rbind',
    #                             aux_betas2)
    #       betas <- as.matrix(c(aux_betas3,
    #                            inte))
    #     } else {
    # 
    #       aux_betas <- lapply(fe,
    #                           as.matrix)
    #       aux_betas2 <- do.call('rbind',
    #                             aux_betas)
    #       betas <- as.matrix(c(aux_betas2,
    #                            inte))
    #     }
    # 
    #     } else {
    # 
    #       betas <- as.matrix(c(unlist(fe), inte))
    # 
    #     }
    # 
    if(length(mu) == 1){#univariado
      betas <- as.matrix(c(mu, 
                           unlist(fe),
                           inte))
    } else {#multivariado
      betas <- rbind(mu,
                     do.call('rbind', 
                             fe),
                     inte)
    } 

    yl <- X%*%betas + Z%*%e_plot + e

    colnames(yl) <- paste('Y',
                          1:dim(yl)[2],
                          sep='')

    Y <- round(yl,
               round)

  } else if(!is.null(blke) & is.null(rowe) & is.null(cole)){#é um DBC

    factors$r <- 1:r  

    ifelse(is.null(blkl),
           factors$Block <- 1:dim(as.matrix(blke))[1],
           factors[[names(blkl)]] <- unlist(blkl))

    dados <- expand.grid(factors,
                         KEEP.OUT.ATTRS=FALSE)

    aux_lf <- names(dados) 
    lb <- aux_lf[length(aux_lf)]#label block

    dados[[lb]] <- factor(dados[[lb]])

    lf <- aux_lf[-c((length(aux_lf)-1):length(aux_lf))] 
    aux_X1 <- paste('~ ',
                    lb,
                    '+',
                    paste(lf,
                          collapse='*')) 

    #     if(is.null(contrasts)){
    #       contrasts <- lapply(factors[aux_lf!='r'],
    #                           function(x)diag(length(x)))
    #     } else{
    #       if((length(fe)+length(blke)) != length(contrasts))
    #         stop('You must be include all the contrasts!')
    #     }
    # 
    #names(contrasts) <- names(dados)[aux_lf!='r']
    auxfactors <- factors[names(factors)!="r"]
    if(quali){#Só qualitativos
      contrast <- lapply(auxfactors,
                         function(x)diag(length(x))) 
    } else if(quanti){#Só fatores quantitativos. O bloco é qualitativo é entra aqui!
      contrast <- lapply(auxfactors,
                         function(x)contr.poly(length(x)))
      contrast[[lb]] <- diag(length(blke))                                              
    } else {#híbrido
      contrast <- lapply(auxfactors,
                         function(x)diag(length(x)))
      contrast[posquanti] <- lapply(auxfactors[posquanti],
                                    function(x)contr.poly(length(x))) 
    }  
    #     contrasts <- lapply(factors[1:length(fe)],
    #                         function(x)diag(length(x)))    
    #}
    names(contrast) <- names(dados)[aux_lf!='r']

    if(!is.null(contrasts)){
      contrast[names(contrasts)] <- contrasts
      contrasts <- contrast
    }else{
      contrasts <- contrast
    }


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
                          'contrasts=FALSE)',
                          collapse=','),
                    ')')

    Z <- model.matrix(eval(parse(text=aux_Z1)), 
                      dados, 
                      contrasts.arg=eval(parse(text=aux_Z2)))

    if(is.null(err)){

      e <- mvtnorm::rmvnorm(n=dim(X)[1],
                            sigma=diag(ncol(as.matrix(fe[[1]]))))

    } else {
      if(!is.matrix(err))
        stop("This argument must be a matrix n x 1 univariate or n x p multivariate!")

      e <- err
    }   

    if(is.null(errp)){

      e_plot <- mvtnorm::rmvnorm(n=dim(Z)[2],
                                 sigma=diag(ncol(as.matrix(fe[[1]]))))

    } else {
      if(!is.matrix(errp))
        stop("This argument must be a matrix n x 1 univariate or n x p multivariate!")

      e_plot <- errp
    }   

    #+ Help to interaction effects!
    #     auxl <- lapply(fe,
    #                    length)
    #     niv <- unlist(auxl)
    #     n <- length(niv)
    #     resl <- list()
    #     for(i in 1:(n-1))resl[[i]] <- combn(niv, i+1)
    #     resl1 <- lapply(resl,
    #                     function(x) apply(x,
    #                                       2,
    #                                       prod))
    #     resl2 <- sum(unlist(resl1))
    # 
    #     if(length(inte)!=resl2){
    #       tcltk::tkmessageBox(message=paste('The number of effects this argument is ', resl2),
    #                           icon='error')
    #     }
    #End

    #     if(length(mu)!=0 & length(mu) == 1){
    # 
    #       betas <- as.matrix(c(mu,
    #                            blke,
    #                            unlist(fe),
    #                            inte))
    # 
    #     } else if(length(mu)!=0 & length(mu) > 1){
    # 
    #       betas <- rbind(mu,
    #                      blke,
    #                      do.call('rbind',
    #                              fe),
    #                      inte)
    # 
    #     } else if(is.null(mu) & length(fe) > 1 & all(unlist(lapply(fl, is.ordered)) == TRUE)) { # Todos os fatores são quantitativos e só há interesse em contrastes polinomiais
    # 
    #       aux_betas <- lapply(fe[-1],
    #                           function(x)x[-1])
    #       aux_betas1 <- c(fe[1],
    #                       aux_betas)
    #       aux_betas2 <- lapply(aux_betas1,
    #                            as.matrix)
    #       aux_betas3 <- do.call('rbind',
    #                             aux_betas2)
    #       betas <- as.matrix(c(blke,
    #                            aux_betas3,
    #                            inte))
    #     } else {
    #       aux_betas <- lapply(fe,
    #                           as.matrix)
    #       aux_betas2 <- do.call('rbind',
    #                             aux_betas)
    #       betas <- as.matrix(c(aux_betas2[1,],
    #                            blke,
    #                            aux_betas2[-1,],
    #                            inte))
    #     }

    #     } else {
    # 
    #       betas <- as.matrix(c(blke, unlist(fe), inte))
    # 
    #     }
    # 
    if(length(mu) == 1){#univariado
      betas <- as.matrix(c(mu,
                           blke, 
                           unlist(fe), 
                           inte))
    } else {#multivariado
      betas <- rbind(mu,
                     blke,
                     do.call('rbind', 
                             fe),
                     inte)
    }  

    yl <- X%*%betas + Z%*%e_plot + e

    colnames(yl) <- paste('Y',
                          1:dim(yl)[2],
                          sep='')

    Y <- round(yl,
               round)

  } else { #é um DQL

    if(dim(as.matrix(fe[[1]]))[1] <= 2)
      stop("The minimum effect this factors is three.")

  aux_trats1 <- suppressWarnings(do.call('interaction', factors))
  aux_trats2 <- sort(levels(aux_trats1))

  nsubplot <- length(factors[[2]])

  aux_trats3 <- matrix(aux_trats2,
                       ncol=nsubplot,
                       byrow=TRUE)

  aux_trats4 <- apply(aux_trats3,
                      1,
                      function(x) paste(x, collapse=' '))

  n <- length(factors[[1]]) # deve ser do mesmo comprimento do fator no primeiro slot da lista (parcela)

  aux_trats5 <- latin(n, 
                      levelss=aux_trats4,
                      nrand=nrand)

  aux_trats7 <- as.matrix(c(aux_trats5))
  aux_trats8 <- apply(aux_trats7,
                      1,
                      function(x) unlist(strsplit(x, ' ')))
  aux_trats9 <- as.matrix(c(aux_trats8))
  aux_trats10 <- strsplit(as.character(aux_trats9[, 1]), '[.]')
  trats <- do.call('rbind',
                   aux_trats10)
  colnames(trats) <- names(factors)

  aux_column <- rep(1:n, rep(nsubplot, n))

  ifelse(is.null(rowl),
         levelsrows <- factor(rep(1:n, rep(length(aux_column), n))),
         levelsrows <- factor(rep(unlist(rowl), rep(length(aux_column), n))))

  ifelse(is.null(coll),
         levelscols <- factor(rep(rep(1:n, rep(nsubplot, n)), n)),
         levelscols <- factor(rep(rep(unlist(coll), rep(nsubplot, n)), n)))

  dados <- data.frame(Row=levelsrows,
                      Column=levelscols,
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

  if(!quali){
    dados[,-c(1:2)][[posquanti]] <- as.ordered(dados[,-c(1:2)][[posquanti]])
  }

  aux_X1 <- paste('~ ',
                  paste(names(dados)[1:2],
                        collapse='+'),
                  '+',
                  paste(names(dados)[-c(1:2)],
                        collapse='*')) 

  #   aufactors <- lapply(dados,
  #                       levels)
  # 
  #   if(is.null(contrasts)){
  #     contrasts <- lapply(aufactors,
  #                         function(x)diag(length(x)))
  #   } else {
  #     if((length(fe)+2) != length(contrasts))  # J.C.Faria
  #       stop('You must be include all the contrasts!')
  #   }
  # 
  #names(contrasts) <- names(dados)
  auxfactors <- lapply(dados,
                       levels)

  if(quali){#Só qualitativos
    contrast <- lapply(auxfactors,
                       function(x)diag(length(x))) 
  } else if(quanti){#Só fatores quantitativos. O bloco é qualitativo é entra aqui!
    contrast <- lapply(auxfactors,
                       function(x)contr.poly(length(x)))
    contrast[[names(dados)[1]]] <- diag(length(cole))
    contrast[[names(dados)[2]]] <- diag(length(rowe))
  } else {#híbrido
    contrast <- lapply(auxfactors,
                       function(x)diag(length(x)))

    contrast[-c(1:2)][posquanti] <- lapply(auxfactors[-c(1:2)][posquanti],
                                           function(x)contr.poly(length(x))) 
  }  

  if(!is.null(contrasts)){
    contrast[names(contrasts)] <- contrasts
    contrasts <- contrast
  }else{
    contrasts <- contrast
  }

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
                         'contrasts=FALSE)',
                         collapse=','),
                   ')') 

  Z <- model.matrix(eval(parse(text=aux_Z1)),
                    dados,
                    contrasts.arg=eval(parse(text=aux_Z2)))
  if(is.null(err)){

    e <- mvtnorm::rmvnorm(n=dim(X)[1],
                          sigma=diag(ncol(as.matrix(fe[[1]]))))

  } else {
    if(!is.matrix(err))
      stop("This argument must be a matrix n x 1 univariate or n x p multivariate!")

    e <- err
  }   

  if(is.null(errp)){

    e_plot <- mvtnorm::rmvnorm(n=dim(Z)[2],
                               sigma=diag(ncol(as.matrix(fe[[1]]))))

  } else {
    if(!is.matrix(errp))
      stop("This argument must be a matrix n x 1 univariate or n x p multivariate!")

    e_plot <- errp
  }   

  #+ Help to interaction effects!
  #     auxl <- lapply(fe,
  #                    length)
  #     niv <- unlist(auxl)
  #     n <- length(niv)
  #     resl <- list()
  #     for(i in 1:(n-1))resl[[i]] <- combn(niv, i+1)
  #     resl1 <- lapply(resl,
  #                     function(x) apply(x,
  #                                       2,
  #                                       prod))
  #     resl2 <- sum(unlist(resl1))
  # 
  #     if(length(inte)!=resl2){
  #       tcltk::tkmessageBox(message=paste('The number of effects this argument is ', resl2),
  #                           icon='error')
  #     }

  #   if(length(mu)!=0 & length(mu) == 1){
  # 
  #     betas <- as.matrix(c(mu,
  #                          rowe,
  #                          cole,
  #                          unlist(fe),
  #                          inte))
  # 
  #   } else if(length(mu)!=0 & length(mu) > 1){
  # 
  #     betas <- rbind(mu,
  #                    rowe,
  #                    cole,
  #                    do.call('rbind',
  #                            fe),
  #                    inte)
  # 
  #   } else if(is.null(mu) & length(fe) > 1 & all(unlist(lapply(fl, is.ordered)) == TRUE)) { # Todos os fatores são quantitativos e só há interesse em contrastes polinomiais
  # 
  #     aux_betas <- lapply(fe[-1],
  #                         function(x)x[-1])
  #     aux_betas1 <- c(fe[1],
  #                     aux_betas)
  #     aux_betas2 <- lapply(aux_betas1,
  #                          as.matrix)
  #     aux_betas3 <- do.call('rbind',
  #                           aux_betas2)
  #     betas <- as.matrix(c(rowe,
  #                          cole,
  #                          aux_betas3,
  #                          inte))
  #   } else {#Temos um híbrido
  #     aux_betas <- lapply(fe,
  #                         as.matrix)
  #     
  #     pos <- which(lapply(fl,is.ordered)==TRUE)
  #     aux_betas2 <- aux_betas[[pos]]#VERIFICAR PARA K fatores. Só está sendo considerado 2.  
  #     aux_betas3 <- aux_betas2[-1,]
  #     aux_betas4 <- matrix(c(aux_betas[[-pos]],aux_betas3),ncol=1)
  #     
  #aux_betas3 <- do.call('rbind',
  #                      aux_betas) 
  #                           
  #     betas <- as.matrix(c(aux_betas2[1,1],
  #                          rowe,
  #                          cole,
  #                          aux_betas4,
  #                          inte))
  #   }

  #     } else {
  #       betas <- as.matrix(c(rowe, cole, do.call('rbind', fe), inte))
  #     }
  # 
  if(length(mu) == 1){#univariado
    betas <- as.matrix(c(mu, 
                         rowe, 
                         cole, 
                         unlist(fe),
                         inte))
  } else {#multivariado
    betas <- rbind(mu,
                   rowe,
                   cole,
                   do.call('rbind', fe),
                   inte)
  }  

  yl <- X%*%betas + Z%*%e_plot + e

  colnames(yl) <- paste('Y',
                        1:dim(yl)[2],
                        sep='')

  Y <- round(yl, round)
  }

  # J.C.Faria
  if(!quali){
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
