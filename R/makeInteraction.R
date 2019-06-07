makeInteraction <- function(mu,
                            fe,
                            inte)
{
  #-----------------------------------------------#
  #++++++++ Help to interaction effects! +++++++++#
  #-----------------------------------------------#
  if(length(mu) == 1){
    auxl <- lapply(fe,
                   length)
  } else {
    aux <- lapply(fe,
                  function(x)apply(x,
                                   2,
                                   length))
    auxl <- lapply(aux,
                   unique)
  }

  niv <- unlist(auxl) 

  n <- length(niv)

  resl <- list()

  for(i in 1:(n-1)){
    resl[[i]] <- combn(niv, 
                       i+1)
  }

  resl1 <- lapply(resl,
                  function(x) apply(x,
                                    2,
                                    prod))

  resl2 <- sum(unlist(resl1))

  if(length(mu) == 1){
    if(is.null(inte)){
      inte <- rep(1,
                  resl2)
    } else {
      if(length(inte) != resl2)    # J.C.Faria
        stop(paste("The length of the 'inte' argument must be: ", 
                   resl2,
                   sep = '')) 
    }
  } else {
    if(is.null(inte)){
      inte <- matrix(rep(1,
                         resl2*length(mu)),
                     ncol = length(mu)) 
    } else {
      if(dim(inte)[1] != resl2)    # J.C.Faria
        stop(paste("The length of the 'inte' argument must be: ", 
                   resl2,
                   sep = ''))  
    }
  }
  return(inte)
}
