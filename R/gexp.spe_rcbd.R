gexp.spe_rcbd <- function(x,
                          ...)
{
  ifelse(is.null(x$fe),
         fe <- list(f1 = rep(1, 
                             3),
                    f2 = rep(1,
                             2)),
         fe <- x$fe)

  ifelse(is.null(x$blke),
         blke <- rep(1, 
                     3),
         blke <- x$blke)

  factors <- list(r = 1:x$r)
  contrast <- list()

  ifelse(is.null(x$blkl),
         {
           factors$Block <- factor(1:dim(as.matrix(blke))[1])
           contrast[["Block"]] <- diag(dim(as.matrix(blke))[1]) 
         },
         {
           factors[[names(x$blkl)]] <- factor(unlist(x$blkl))
           contrast[[names(x$blkl)]] <- diag(dim(as.matrix(blke))[1]) 
         })                                

  intee <- makeInteraction(mu = x$mu,
                           fe = fe,
                           inte = x$inte)#Preparando o argumento inte dos efeitos da interação

  treatments <- makeTreatments(fl = x$fl,
                               fe = fe,
                               quali = x$qualiquanti$quali,
                               quanti = x$qualiquanti$quanti,
                               posquanti = x$qualiquanti$posquanti)

  contrasttreatments <- makeContrasts(factors = treatments,
                                      quali = x$qualiquanti$quali,
                                      quanti = x$qualiquanti$quanti,
                                      posquanti = x$qualiquanti$posquanti)

  contrast <- c(contrast,
                contrasttreatments)  

  if(!is.null(x$contrasts)){
    contrast[names(x$contrasts)] <- x$contrasts
    contrasts <- contrast
  }else{
    contrasts <- contrast
  }

  factors <- c(factors,
               treatments)

  cformula <- paste('~',
                    names(factors)[2],
                    '+',
                    paste(names(treatments),
                          collapse = '*'))                         

  dados <- expand.grid(factors,
                       KEEP.OUT.ATTRS = FALSE)#montando o data.frame

  XB <- makeXBeta(cformula, 
                  dados, 
                  mu = x$mu, 
                  fe = fe,
                  blke = blke, 
                  rowe = x$rowe,
                  cole = x$cole, 
                  inte = intee, 
                  contrasts = contrasts)

  facplott <- names(treatments)[1]

  cformulaplot <- paste('~ 0 + ',
                        paste(facplott,
                              names(factors)[2],
                              sep = ':'))

  plott <- c(facplott, 
             names(factors)[2])

  zformula <-  paste('list(',
                     paste(plott,
                           paste('= contrasts(dados$',
                                 plott,
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
