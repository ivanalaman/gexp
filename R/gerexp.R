gerexp <- function(mu    = NULL, #é um escalar numérico ou um vetor se tiver mais de uma variável resposta
                   s2    = 1L,   #é um escalar numérico ou uma matrix de variância-covariância para dados multivariados
                   s2sp  = 1L,   #é um escalar numérico ou uma matrix de variância-covariância para dados multivariados no caso de type ser igual a SPLIT
                   r     = 2L,   #é um escalar que se refere ao nº de repetições
                   f     = NULL, #é uma lista em que cada dimensão se refere ao efeito de um fator que pode ser um vetor numérico ou uma matriz (caso multivariado)
                   inter = NULL, #é um vetor numérico ou uma matriz (caso multivariado) que se refere aos efeitos da interação
                   b     = NULL, #é um vetor numérico ou uma matriz (caso multivariado) que se refere aos efeitos dos blocos
                   erow  = NULL, #é um vetor numérico ou uma matriz (caso multivariado) que se refere aos efeitos das linhas
                   ecol  = NULL, #é um vetor numérico ou uma matriz (caso multivariado) que se refere aos efeitos das colunas
                   type  = c('DIC','DBC','DQL','FAT','SPLIT'),
                   rd    = 2L,   # é um escalar numérico
                   ...) 
{

  if(is.null(f)) stop("You must specify at least a factor") 

  s2 <- as.matrix(s2)

  aux_factor <- lapply(f,
                       function(x) as.matrix(x))

  names(aux_factor) <- paste('X',1:length(f),sep='')

  aux_factor1 <- as.list(tolower(names(aux_factor)))
  aux_factor2 <- lapply(aux_factor,
                        function(x) 1:dim(x)[1])

  factors <- mapply(function(x,y) paste(x,y,sep=''),
                    aux_factor1,
                    aux_factor2,
                    SIMPLIFY = FALSE)

  names(factors) <- names(aux_factor)

  switch(match.arg(type),
         DIC = {        


           if(any(!is.null(inter) | !is.null(b) | !is.null(erow) | !is.null(ecol))) stop("Only f must be specified as effect of treatment!")

           factors$r <- 1:r

           dados <- expand.grid(factors,
                                KEEP.OUT.ATTRS = FALSE)

           aux_X1 <- paste('~',
                           paste(names(dados)[-length(names(dados))],
                                 collapse='+')) 

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
           #print(X)

           e <- mvtnorm::rmvnorm(n = dim(X)[1],  
                                 sigma = s2)

           ifelse(length(mu) == 1,
                  betas <- as.matrix(c(mu, unlist(f))),
                  betas <- rbind(mu, do.call('rbind',f)))

           yl <- X%*%betas + e

           colnames(yl) <- paste('Y',1:dim(yl)[2],sep='')

           Y <- round(yl, rd)

         },
         DBC = {

           if(is.null(b)) stop("You must specify the block effect with b argument!")

           factors$Block <- 1:dim(as.matrix(b))[1]

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
                                 sigma = s2)

           ifelse(length(mu) == 1,
                  betas <- as.matrix(c(mu, b, unlist(f))),
                  betas <- rbind(mu, b, do.call('rbind',f)))

           yl <- X%*%betas + e

           colnames(yl) <- paste('Y',1:dim(yl)[2],sep='') 

           Y <- round(yl, rd)

         },
         DQL = {

           if(length(f)!=1) stop("Use only one factor this design!")

           f1 <- f[[1]]
           n <- length(f1)

           aux_dados  <- latin(n)

           dados = data.frame(Row    = factor(rep(1:n,rep(n,n))),
                              Column = factor(rep(1:n,n)),
                              X1   = c(aux_dados))

           X = model.matrix(~ Row + Column + X1, 
                            dados,
                            contrasts.arg = list(Row=contrasts(dados$Row,
                                                               contrasts = FALSE),
                                                 Column = contrasts(dados$Column,
                                                                    contrasts=FALSE),
                                                 X1 = contrasts(dados$X1,
                                                                contrasts=FALSE)))

           e <- mvtnorm::rmvnorm(n = dim(X)[1],  
                                 sigma = s2)

           ifelse(length(mu) == 1,
                  betas <- as.matrix(c(mu, erow, ecol, f1)),
                  betas <- rbind(mu, erow, ecol, f1))

           yl <- X%*%betas + e

           colnames(yl) <- paste('Y',1:dim(yl)[2],sep='')  

           Y <- round(yl, rd)

         },
         FAT = {

           if(is.null(b) & is.null(erow) & is.null(ecol)){ #é um DIC

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
                                   sigma = s2)

             ifelse(length(mu) == 1,
                    betas <- as.matrix(c(mu, unlist(f), inter)),
                    betas <- rbind(mu, do.call('rbind',f), inter)) 

             yl <- X%*%betas + e

             colnames(yl) <- paste('Y',1:dim(yl)[2],sep='')

             Y <- round(yl, rd)

           } else if(!is.null(b) & is.null(erow) & is.null(ecol)){#é um DBC

             factors$Block <- 1:dim(as.matrix(b))[1]

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
                                   sigma = s2)

             ifelse(length(mu) == 1,
                    betas <- as.matrix(c(mu, b, unlist(f), inter)),
                    betas <- rbind(mu, b, do.call('rbind',f), inter))

             yl <- X%*%betas + e

             colnames(yl) <- paste('Y',1:dim(yl)[2],sep='') 

             Y <- round(yl, rd)

           } else { #é um DQL

             aux_trats1 <- do.call('interaction',factors)
             aux_trats2 <- levels(aux_trats1)

             n <- length(aux_trats2)
             aux_trats3 <- latin(n)
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
                                   sigma = s2)

             ifelse(length(mu) == 1,
                    betas <- as.matrix(c(mu, erow, ecol, unlist(f), inter)),
                    betas <- rbind(mu, erow, ecol, do.call('rbind',f), inter))

             yl <- X%*%betas + e

             colnames(yl) <- paste('Y',1:dim(yl)[2],sep='')

             Y <- round(yl, rd)
             
           }
         },
         SPLIT = {

           s2sp <- as.matrix(s2sp)

           if(is.null(b) & is.null(erow) & is.null(ecol)){ #é um DIC

             factors$r <- factor(1:r)

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
             Z <- model.matrix(~ 0 + X1:r, 
                               dados, 
                               contrasts.arg = list(X1 = contrasts(dados$X1, contrasts = FALSE),
                                                    r = contrasts(dados$r, contrasts = FALSE)))

             e <- mvtnorm::rmvnorm(n = dim(X)[1],  
                                   sigma = s2)

             e_plot <- mvtnorm::rmvnorm(n = dim(Z)[2],  
                                        sigma = s2sp)

             ifelse(length(mu) == 1,
                    betas <- as.matrix(c(mu, unlist(f), inter)),
                    betas <- rbind(mu, do.call('rbind',f), inter))

             yl <- X%*%betas + Z%*%e_plot + e

             colnames(yl) <- paste('Y',1:dim(yl)[2],sep='')

             Y <- round(yl, rd)

           } else if(!is.null(b) & is.null(erow) & is.null(ecol)){#é um DBC

             factors$Block <- 1:dim(as.matrix(b))[1]

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

             Z <- model.matrix(~ 0 + X1:Block, 
                               dados, 
                               contrasts.arg = list(X1 = contrasts(dados$X1, contrasts = FALSE),
                                                    Block = contrasts(dados$Block, contrasts = FALSE)))

             e <- mvtnorm::rmvnorm(n = dim(X)[1],  
                                   sigma = s2)

             e_plot <- mvtnorm::rmvnorm(n = dim(Z)[2],  
                                        sigma = s2sp)

             ifelse(length(mu) == 1,
                    betas <- as.matrix(c(mu, b, unlist(f), inter)),
                    betas <- rbind(mu, b, do.call('rbind',f), inter))

             yl <- X%*%betas + Z%*%e_plot + e

             colnames(yl) <- paste('Y',1:dim(yl)[2],sep='')

             Y <- round(yl, rd)

         } else { #é um DQL

           if(dim(as.matrix(f[[1]]))[1] <= 2) stop("The minimum effect this factors is three.")
           
           aux_trats1 <- suppressWarnings(do.call('interaction',factors))
           aux_trats2 <- sort(levels(aux_trats1))
          
           nsubplot <- length(factors[[2]])

           aux_trats3 <- matrix(aux_trats2,
                                ncol=nsubplot,
                                byrow=TRUE)

           aux_trats4 <- apply(aux_trats3,1,function(x)paste(x,collapse=' '))
           n <- length(factors[[1]])# deve ser do mesmo comprimento do fator no primeiro slot da lista (parcela)

           aux_trats5 <- latin(n)

           aux_trats6 <- mgsub(LETTERS[1:n],aux_trats4,aux_trats5)
           aux_trats7 <- as.matrix(c(aux_trats6))
           aux_trats8 <- apply(aux_trats7, 1, function(x)unlist(strsplit(x,' ')))
           aux_trats9 <- as.matrix(c(aux_trats8))
           aux_trats10 <- strsplit(as.character(aux_trats9[,1]),'[.]')
           trats <- do.call('rbind',aux_trats10)

           aux_column <- rep(1:n, rep(nsubplot,n))
           roww <- rep(1:n, rep(length(aux_column),n))
           column <- rep(aux_column,n)

           dados  <- data.frame(Row    = factor(roww),
                                Column = factor(column),
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

           Z <- model.matrix(~ 0 + X1:Row:Column, 
                             dados, 
                             contrasts.arg = list(X1 = contrasts(dados$X1, contrasts = FALSE),
                                                  Row = contrasts(dados$Row, contrasts = FALSE),
                                                  Column = contrasts(dados$Column, contrasts = FALSE))) 

           e <- mvtnorm::rmvnorm(n = dim(X)[1],  
                                 sigma = s2)

           e_plot <- mvtnorm::rmvnorm(n = dim(Z)[2],  
                                      sigma = s2sp)

           ifelse(length(mu) == 1,
                  betas <- as.matrix(c(mu, erow, ecol, unlist(f), inter)),
                  betas <- rbind(mu, erow, ecol, do.call('rbind',f), inter))

           yl <- X%*%betas + Z%*%e_plot + e

           colnames(yl) <- paste('Y',1:dim(yl)[2],sep='')

           Y <- round(yl, rd)

           }
        }
  )

  dados <- cbind(dados,Y)

  return(dados)
}
