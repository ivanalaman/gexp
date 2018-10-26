gerexp.fe <- function(mu = mu,
                      sigma = sigma, 
                      r = r,   
                      ef = ef, 
                      einter = einter, 
                      eb = eb, 
                      erow = erow, 
                      ecol = ecol, 
                      nrand = nrand,    
                      rd = rd,   
                      randomized = randomized) 
{
  if(is.null(ef)) stop("You must specify at least a factor") 

  sigma <- as.matrix(sigma)

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

  if(is.null(eb) & is.null(erow) & is.null(ecol)){ #é um DIC

    factors$r <- 1:r

    dados <- expand.grid(factors,
                         KEEP.OUT.ATTRS = FALSE)

    aux_X1 <- paste('~',
                    paste(names(dados)[-length(names(dados))],
                          collapse='*')) 

    aux_X2 <- paste('list(',
                    paste(names(factors)[-length(names(dados))],
                          paste('= contrasts(dados$',
                                names(factors)[-length(names(dados))],
                                sep=''),
                          ',',
                          'contrasts = FALSE)',
                          collapse=','),
                    ')')

    X  <- model.matrix(eval(parse(text=aux_X1)),
                       dados,
                       contrasts.arg = eval(parse(text=aux_X2)))             

    e <- mvtnorm::rmvnorm(n = dim(X)[1],  
                          sigma = sigma)

    ifelse(length(mu) == 1,
           betas <- as.matrix(c(mu, unlist(ef), einter)),
           betas <- rbind(mu, do.call('rbind',ef), einter)) 

    yl <- X%*%betas + e

    colnames(yl) <- paste('Y',1:dim(yl)[2],sep='')

    Y <- round(yl, rd)

  } else if(!is.null(eb) & is.null(erow) & is.null(ecol)){#é um DBC

    factors$Block <- 1:dim(as.matrix(eb))[1]

    dados <- expand.grid(factors,
                         KEEP.OUT.ATTRS = FALSE)

    dados$Block <- factor(dados$Block)

    aux_X1 <- paste('~ Block +',
                    paste(names(dados)[-length(names(dados))],
                          collapse='*')) 

    aux_X2 <- paste('list(',
                    paste(names(factors),
                          paste('= contrasts(dados$',
                                names(factors),
                                sep=''),
                          ',',
                          'contrasts = FALSE)',
                          collapse=','),
                    ')')

    X = model.matrix(eval(parse(text=aux_X1)),
                     dados,
                     contrasts.arg = eval(parse(text=aux_X2)))             

    e <- mvtnorm::rmvnorm(n = dim(X)[1],  
                          sigma = sigma)

    ifelse(length(mu) == 1,
           betas <- as.matrix(c(mu, eb, unlist(ef), einter)),
           betas <- rbind(mu, eb, do.call('rbind',ef), einter))

    yl <- X%*%betas + e

    colnames(yl) <- paste('Y',1:dim(yl)[2],sep='') 

    Y <- round(yl, rd)

  } else { #é um DQL

    aux_trats1 <- do.call('interaction',factors)
    aux_trats2 <- levels(aux_trats1)

    n <- length(aux_trats2)
    aux_trats3 <- latin(n, nrand = nrand)
    aux_trats4 <- mgsub(LETTERS[1:n],aux_trats2,aux_trats3)
    aux_trats5 <- as.matrix(c(aux_trats4))
    aux_trats6 <- strsplit(as.character(aux_trats5[,1]),'[.]')
    trats <- do.call('rbind',aux_trats6)

    dados  <- data.frame(Row    = factor(rep(1:n,rep(n,n))),
                         Column = factor(rep(1:n,n)),
                         trats = trats)
    names(dados) <- c('Row','Column',paste('X',
                                           1:dim(trats)[2],
                                           sep=''))

    aux_X1 <- paste('~ Row + Column +',
                    paste(names(dados)[-c(1:2)],
                          collapse='*')) 

    aux_X2 <- paste('list(',
                    paste(names(dados),
                          paste('= contrasts(dados$',
                                names(dados),
                                sep=''),
                          ',',
                          'contrasts = FALSE)',
                          collapse=','),
                    ')')

    X = model.matrix(eval(parse(text=aux_X1)),
                     dados,
                     contrasts.arg = eval(parse(text=aux_X2)))             

    e <- mvtnorm::rmvnorm(n = dim(X)[1],  
                          sigma = sigma)

    ifelse(length(mu) == 1,
           betas <- as.matrix(c(mu, erow, ecol, unlist(ef), einter)),
           betas <- rbind(mu, erow, ecol, do.call('rbind',ef), einter))

    yl <- X%*%betas + e

    colnames(yl) <- paste('Y',1:dim(yl)[2],sep='')

    Y <- round(yl, rd)

  }

  Z <- NULL

  dados <- cbind(dados,Y)

  if(randomized){

    dados <- dados[sample(dim(dados)[1]),]

  }

  res <- list(X   = X,
              Z   = Z,
              Y   = Y,
              dfm = dados)

}   
