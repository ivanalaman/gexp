gexp.lsd <- function(mu        = mu, 
                     err       = err,
                     fl        = fl,
                     rowl      = rowl,
                     coll      = coll,
                     fe        = fe,
                     rowe      = rowe,
                     cole      = cole,
                     nrand     = nrand,    
                     contrasts = contrasts,
                     round     = round, 
                     random    = random, ...)
{

  if(length(fe) != 1){
    stop('Use only one factor!')
  }

  f1 <- fe[[1]]
  n <- length(f1)

  if(is.null(fl)){
    levelss <- paste('x', 1:n, sep='')
  } else {
    if(length(fl)!=1){
      stop('Use only one factor!')
    }

    levelss <- unlist(fl)
  }

  ifelse(is.null(rowl),
    levelsrows <- factor(rep(1:n, rep(n, n))),
    levelsrows <- factor(rep(unlist(rowl), rep(n, n))))

  ifelse(is.null(coll),
         levelscols <- factor(rep(1:n, n)),
         levelscols <- factor(rep(unlist(coll), n)))

  aux_dados  <- latin(n,
                      levelss = levelss,
                      nrand = nrand)

  dados = data.frame(Row    = levelsrows,
                     Column = levelscols,
                     X1     = c(aux_dados))

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
   
  if(!is.null(fl)){
    names(dados) <- gsub('X1',
                         names(fl),
                         names(dados))
  }
  
  factors <- lapply(dados, levels)

  aux_X1 <- paste('~ ',
                  paste(names(dados),
                        collapse='+'))

  if(is.null(contrasts)){
    contrasts <- lapply(factors, function(x)diag(length(x)))
  } else{
    if((length(fe)+length(rowe)+length(cole)) != length(contrasts))stop('You must be include all the contrasts!')
  }

  names(contrasts) <- names(dados)

  X  <- model.matrix(eval(parse(text=aux_X1)),
                     dados,
                     contrasts.arg=contrasts) 

  Z <- NULL
  
  if(is.null(err)){

   e <- mvtnorm::rmvnorm(n = dim(X)[1],  
                          sigma = diag(ncol(as.matrix(fe[[1]]))))

  } else {
    if(!is.matrix(err)) stop("This argument must be a matrix n x 1 univariate or n x p multivariate!")

    e <- err
  } 
 
  if(length(mu)!=0 & length(mu) == 1){

    betas <- as.matrix(c(mu, rowe, cole, f1))

  } else if(length(mu)!=0 & length(mu) > 1){

    betas <- rbind(mu, rowe, cole, f1)

 } else if(is.null(mu) & length(fe) > 1 & all(unlist(lapply(fl, is.ordered))==TRUE)){#Todos os fatores são quantitativos e só há interesse em contrastes polinomiais

    aux_betas <- lapply(fe[-1],
                        function(x)x[-1])
    aux_betas1 <- c(fe[1], aux_betas)
    aux_betas2 <- lapply(aux_betas1, as.matrix)
    aux_betas3 <- do.call('rbind', aux_betas2)
    betas <- as.matrix(c(rowe, cole, aux_betas3))

  } else {

    aux_betas <- lapply(fe, as.matrix)
    aux_betas2 <- do.call('rbind', aux_betas)
    betas <- as.matrix(c(rowe, cole, aux_betas2))
    
  }
   
    #   } else {
    #     betas <- as.matrix(c(rowe, cole, f1))
    #   }

  yl <- X%*%betas + e

  colnames(yl) <- paste('Y', 1:dim(yl)[2], sep='')  

  Y <- round(yl, round)

  # J.C.Faria
  if(is.null(mu)){
    dados <- lapply(dados, 
                    function(x) if(is.ordered(factor(x))) as.numeric(as.character(x)) else x)

    dados <- as.data.frame(dados)                
  }                  

  dados <- cbind(dados, Y)

  if(random){

    dados <- dados[sample(dim(dados)[1]),]

  }

  res <- list(X   = X,
              Z   = Z,
              Y   = Y,
              dfm = dados)

} 
