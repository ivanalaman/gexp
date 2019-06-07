gexp.simple_crd <- function(x,
                            ...)
{
  ifelse(is.null(x$fe),
         fe <- list(f1 = rep(1, 
                             3)),
         fe <- x$fe)

  factors <- list(r = 1:x$r)

  treatments <- makeTreatments(fl = x$fl,
                               fe = fe,
                               quali = x$qualiquanti$quali,
                               quanti = x$qualiquanti$quanti,
                               posquanti = x$qualiquanti$posquanti)

  contrast <- makeContrasts(factors = treatments,
                            quali = x$qualiquanti$quali,
                            quanti = x$qualiquanti$quanti,
                            posquanti = x$qualiquanti$posquanti)

  if(!is.null(x$contrasts)){
    contrast[names(x$contrasts)] <- x$contrasts
    contrasts <- contrast
  }else{
    contrasts <- contrast
  }

  cformula <- paste('~',
                    paste(names(treatments),
                          collapse = '+'))                         

  factors <- c(factors,
               treatments)

  dados <- expand.grid(factors,
                       KEEP.OUT.ATTRS = FALSE)

  XB <- makeXBeta(cformula, 
                  dados, 
                  mu = x$mu, 
                  fe = fe,
                  blke = x$blke, 
                  rowe = x$rowe,
                  cole = x$cole, 
                  inte = x$inte, 
                  contrasts = contrasts)

  Z <- NULL

  if(is.null(x$err)){
    e <- mvtnorm::rmvnorm(n = dim(XB$XB)[1],
                          sigma = diag(length(x$mu)))
  }else{
    if(!is.matrix(x$err))
      stop("This argument must be a matrix n x 1 univariate or n x p multivariate!")

    e <- x$err
  }

  yl <- XB$XB + e

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
