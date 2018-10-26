gerexp.lsd <- function(mu = mu, 
                       sigma = sigma, 
                       ef = ef,  
                       erow = erow, 
                       ecol = ecol, 
                       nrand = nrand,    
                       rd = rd,   
                       randomized = randomized)
{

  sigma <- as.matrix(sigma)

  f1 <- ef[[1]]
  n <- length(f1)

  aux_dados  <- latin(n, nrand = nrand)

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

  Z <- NULL

  e <- mvtnorm::rmvnorm(n = dim(X)[1],  
                        sigma = sigma)

  ifelse(length(mu) == 1,
         betas <- as.matrix(c(mu, erow, ecol, f1)),
         betas <- rbind(mu, erow, ecol, f1))

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
