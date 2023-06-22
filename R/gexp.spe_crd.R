gexp.spe_crd <- function(x,
                         ...)
{
  ifelse(is.null(x$fe),
         fe <- list(f1 = rep(1, 
                             3),
                    f2 = rep(1,
                             2)),
         fe <- x$fe)

  factors <- list(r = factor(1:x$r))

  intee <- makeInteraction(mu = x$mu,
                           fe = fe,
                           inte = x$inte)#Preparando o argumento inte dos efeitos da interacao

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
                          collapse = '*'))                         

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
                  inte = intee, 
                  contrasts = contrasts)

  facplott <- names(treatments)[1]

  cformulaplot <- paste('~ 0 + ',
                        facplott,
                        ':r',
                        sep = '')

  plott <- c(facplott, 'r')
  zformula <-  paste('list(',
                     paste(plott,
                           paste('= contrasts(dados$',
                                 plott,
                                 sep = ''),
                           ',',
                           'contrasts=FALSE)',
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

  res <- list(X   = XB$X,
              Z   = Z,
              Y   = Y,
              dfm = dados)

  class(res) <- c(paste('gexp',
                        class(x),
                        sep = '.'),
                  'gexp',
                  'list')
  
  return(res)  
}
