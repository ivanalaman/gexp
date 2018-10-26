gerexp.crd <- function(mu = mu,
                       sigma = sigma,
                       r = r,
                       ef = ef,
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
  Z <- NULL

  e <- mvtnorm::rmvnorm(n = dim(X)[1],  
                        sigma = sigma)

  ifelse(length(mu) == 1,
         betas <- as.matrix(c(mu, unlist(ef))),
         betas <- rbind(mu, do.call('rbind',ef)))

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
