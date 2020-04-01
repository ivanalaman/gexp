gexp.spe_lsd <- function(x,
                         ...)
{
  ifelse(is.null(x$fe),
         fe <- list(f1 = rep(1, 
                             3),
                    f2 = rep(1,
                             2)),
         fe <- x$fe)

  if(x$r != 1){
    x$r <- 1
    warning('Internaly replicates was set to one (r=1)!')
  }

  intee <- makeInteraction(mu = x$mu,
                           fe = fe,
                           inte = x$inte)#Preparando o argumento inte dos efeitos da interação 

  treatments <- makeTreatments(fl = x$fl,
                               fe = fe,
                               quali = x$qualiquanti$quali,
                               quanti = x$qualiquanti$quanti,
                               posquanti = x$qualiquanti$posquanti)

  contrast <- makeContrasts(factors = treatments,
                            quali = x$qualiquanti$quali,
                            quanti = x$qualiquanti$quanti,
                            posquanti = x$qualiquanti$posquanti)

  n <- length(treatments[[1]])

  ifelse(is.null(x$rowe),
         rowe <- rep(1,n),
         rowe <- x$rowe)

  ifelse(is.null(x$cole),
         cole <- rowe,
         cole <- x$cole)

  rowcolumn <- list()

  ifelse(is.null(x$rowl),
         {
           rowcolumn$Row <- factor(1:dim(as.matrix(rowe))[1])   
           contrast[["Row"]] <- diag(dim(as.matrix(rowe))[1])  
         },
         {
           rowcolumn[[names(x$rowl)]] <-factor(unlist(x$rowl))
           contrast[[names(x$rowl)]] <- diag(dim(as.matrix(rowe))[1])    
         })

  ifelse(is.null(x$coll),
         {
           rowcolumn$Column <- factor(1:dim(as.matrix(cole))[1])  
           contrast[["Column"]] <- diag(dim(as.matrix(cole))[1])
         },
         {
           rowcolumn[[names(x$coll)]] <- factor(unlist(x$coll))
           contrast[[names(x$coll)]] <- diag(dim(as.matrix(cole))[1])
         })

  if(!is.null(x$contrasts)){
    contrast[names(x$contrasts)] <- x$contrasts
    contrasts <- contrast
  }else{
    contrasts <- contrast
  }

  combfactors <- suppressWarnings(do.call('interaction', 
                                          treatments)) 

  scombfactors <- sort(levels(combfactors))

  nsubplot <- do.call('prod',
                      lapply(treatments[-1],
                             length))

  mcombfac <- matrix(scombfactors,
                     ncol = nsubplot,
                     byrow = TRUE)

  levelsfactors <- apply(mcombfac,
                         1,
                         function(x) paste(x, 
                                           collapse = ' '))

  n <- length(treatments[[1]]) # deve ser do mesmo comprimento do fator no primeiro slot da lista (parcela)

  sorttreatment <- latin(n, 
                         levelss = levelsfactors,
                         nrand = 0)

  msortfac <- as.matrix(c(sorttreatment))

  aux_trats8 <- apply(msortfac,
                      1,
                      function(x) unlist(strsplit(x, 
                                                  ' ')))
  aux_trats9 <- as.matrix(c(aux_trats8))

  aux_trats10 <- strsplit(as.character(aux_trats9[, 1]), 
                          '[.]')

  trats <- do.call('rbind',
                   aux_trats10)

  colnames(trats) <- names(treatments)

  trats <- as.data.frame(trats, stringsAsFactors=TRUE)

  if(!x$qualiquanti$quali){
    trats[,x$qualiquanti$posquanti] <- as.ordered(trats[,x$qualiquanti$posquanti])
  }

  rowcolumn[[1]] <- rep(rowcolumn[[1]],
                        rep(length(scombfactors),
                            n))
  rowcolumn[[2]] <- rep(rep(rowcolumn[[2]],
                            rep(nsubplot,n)),
                        n)

  combrowcolumn <- data.frame(rowcolumn)

  dados <- data.frame(combrowcolumn,
                      trats)

  cformula <- paste('~',
                    paste(names(rowcolumn),
                          collapse = '+'),
                    '+',
                    paste(names(treatments),
                          collapse = '*'))                         

  XB <- makeXBeta(cformula, 
                  dados, 
                  mu = x$mu, 
                  fe = fe,
                  blke = x$blke, 
                  rowe = rowe,
                  cole = cole, 
                  inte = intee, 
                  contrasts = contrasts)

  facplott <- colnames(trats)[1]

  cformulaplot <- paste('~ 0 + ',
                        facplott,
                        ':',
                        paste(names(combrowcolumn),
                              collapse = ':'),
                        sep = '')

  parceprincipal <- c(facplott, 
                      names(combrowcolumn))

  zformula <-  paste('list(',
                     paste(parceprincipal,
                           paste('= contrasts(dados$',
                                 parceprincipal,
                                 sep = ''),
                           ',',
                           'contrasts = FALSE)',
                           collapse = ','),
                     ')') 

  Z <- model.matrix(eval(parse(text = cformulaplot)),
                    dados,
                    contrasts.arg = eval(parse(text = zformula)))

  if(is.null(x$err)){
    e <- mvtnorm::rmvnorm(n = dim(XB$XB)[1],
                          sigma = diag(length(x$mu)))
  }else{
    if(!is.matrix(x$err))
      stop("This argument must be a matrix n x 1 univariate or n x p multivariate!")

    e <- x$err
  }

  if(is.null(x$errp)){
    e_plot <- mvtnorm::rmvnorm(n = dim(Z)[2],
                               sigma = diag(length(x$mu)))
  } else {
    if(!is.matrix(x$errp))
      stop("This argument must be a matrix n x 1 univariate or n x p multivariate!")

    e_plot <- x$errp
  }   

  yl <- XB$XB + Z%*%e_plot + e

  colnames(yl) <- paste('Y',
                        1:dim(yl)[2],
                        sep = '')

  Y <- round(yl,
             x$round)
  #J.C.Faria
  if(!x$qualiquanti$quali){
    dados <- lapply(dados, 
                    function(x) if(is.ordered(factor(x))) as.numeric(as.character(x)) else x)

    dados <- as.data.frame(dados)                
  }                  

  dados <- cbind(dados,
                 Y)

  res <- list(X = XB$X,
              Z = Z,
              Y = Y,
              dfm = dados)

  class(res) <- c(paste('gexp',
                        class(x),
                        sep = '.'),
                  'gexp',
                  'list')

  return(res)
}
